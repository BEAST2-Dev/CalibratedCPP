package calibratedcpp.lphy.tree;

import calibratedcpp.lphy.prior.Calibration;
import calibratedcpp.lphy.prior.CalibrationArray;
import lphy.base.distribution.DistributionConstants;
import lphy.base.evolution.birthdeath.BirthDeathConstants;
import lphy.base.evolution.tree.TaxaConditionedTreeGenerator;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;

import java.util.*;

import static calibratedcpp.lphy.tree.CPPUtils.*;

/**
 * Shared calibrated-CPP tree-building logic (clade decomposition, weighted slot
 * assignment, coalescing) for a coalescent point process whose node ages come
 * from some birth-death law. Concrete subclasses only need to supply that law
 * via {@link #logCDF} (the node-age log-CDF); the truncated age sampling, stem-age
 * sampling and the uncalibrated fallback are all derived from it here, mirroring how
 * the BEAST side derives everything from {@code calculateLogNodeAgeCDF}
 * (see {@code CalibratedCoalescentPointProcess}).
 */
public abstract class AbstractCalibratedCPPTree extends TaxaConditionedTreeGenerator implements GenerativeDistribution<TimeTree> {

    Value<Number> rho;
    Value<CalibrationArray> calibrations;
    Value<String[]> otherNames;
    Value<Number> stemAge;
    Value<Number> rootAge;
    double conditionAge;
    boolean rootConditioned = false;

    List<String> nameList;
    Set<String> usedNames;

    public static final String calibrationsName = "calibrations";
    public static final String stemAgeName = "stemAge";
    public static final String otherTaxaNames = "otherNames";
    public static final String rootAgeName = "rootAge";

    protected AbstractCalibratedCPPTree(Value<Integer> n, Value<Number> rho, Value<CalibrationArray> calibrations,
                                         Value<String[]> otherNames, Value<Number> stemAge, Value<Number> rootAge) {
        super(n, null, null);

        // No conditioning is required: with neither calibrations nor rootAge, the uncalibrated path
        // samples the stem (fixed n) or the origin (random n) from the process itself.
        this.rho = rho;
        this.calibrations = calibrations;
        this.otherNames = otherNames;
        this.stemAge = stemAge;
        this.rootAge = rootAge;
    }

    // ****** node age density hooks — the only thing concrete subclasses implement ******

    /** Recompute any resolved rate values from the current parameter bindings. Called once per sample(). */
    protected abstract void resolveRates();

    protected abstract double logCDF(double t);

    protected double cdf(double t) { return Math.exp(logCDF(t)); }

    /**
     * Draw nSims node ages from the law truncated to [lowerTime, upperTime] by
     * inverse-transform sampling of {@link #cdf}. This generic version works for any
     * model that can compute {@link #logCDF}; the constant-rate subclasses override it
     * with a tuned closed-form sampler (see {@code CPPUtils.sampleTimes}).
     */
    protected double[] sampleAges(double lowerTime, double upperTime, int nSims) {
        double qLo = cdf(lowerTime);
        double qHi = Double.isInfinite(upperTime) ? 1.0 : cdf(upperTime);
        double[] out = new double[nSims];
        for (int i = 0; i < nSims; i++) {
            out[i] = invertCDF(qLo + random.nextDouble() * (qHi - qLo));
        }
        return out;
    }

    /**
     * Sample a stem/origin age older than greaterThan for nTaxa taxa. The origin has
     * CDF Q(t)^nTaxa (the max of nTaxa i.i.d. node ages), left-truncated at greaterThan.
     * Derived from {@link #cdf}; constant-rate subclasses override with a closed form.
     */
    protected double sampleStemAge(double greaterThan, int nTaxa) {
        double qaN = Math.pow(cdf(greaterThan), nTaxa);
        double p = qaN + random.nextDouble() * (1.0 - qaN);
        return invertCDF(Math.pow(p, 1.0 / nTaxa));
    }

    /**
     * Numerical inverse of the monotone {@link #cdf}: solve Q(t) = p for t &gt;= 0 via an
     * exponential bracket followed by bisection. Override with a closed form when the model
     * has one. (Bisection is the robust starting point; swap for Brent once it works.)
     */
    protected double invertCDF(double p) {
        if (p <= 0.0) return 0.0;
        double hi = 1.0;
        while (cdf(hi) < p && hi < 1e12) hi *= 2.0;
        double lo = 0.0;
        for (int i = 0; i < 100; i++) {
            double mid = 0.5 * (lo + hi);
            if (cdf(mid) < p) lo = mid; else hi = mid;
        }
        return 0.5 * (lo + hi);
    }

    /** Build nTips leaf nodes: supplied taxa names padded with indices, shuffled with the framework RNG. */
    private List<TimeTreeNode> buildShuffledLeaves(String[] taxaNames, int nTips) {
        List<String> names = new ArrayList<>(nTips);
        if (taxaNames != null) {
            for (int i = 0; i < Math.min(taxaNames.length, nTips); i++) names.add(taxaNames[i]);
        }
        for (int i = names.size(); i < nTips; i++) names.add(String.valueOf(i));
        for (int i = names.size() - 1; i > 0; i--) {          // Fisher-Yates on the seeded RNG
            int j = random.nextInt(i + 1);
            String tmp = names.get(i); names.set(i, names.get(j)); names.set(j, tmp);
        }
        List<TimeTreeNode> nodes = new ArrayList<>(nTips);
        for (String name : names) {
            TimeTreeNode leaf = new TimeTreeNode(0.0);
            leaf.setId(name);
            nodes.add(leaf);
        }
        return nodes;
    }

    /**
     * Fixed n, fixed root: n-2 internal depths from q/Q(rootAge) plus the root age placed uniformly
     * among them (so it is the maximum depth). Root-conditioned, no stem. Node ages come from this
     * generator's own law via {@link #sampleAges}, so it works for any subclass.
     */
    private RandomVariable<TimeTree> sampleUncalibratedTree(String[] taxaNames, double rootAge) {
        int n = getN().value();
        conditionAge = rootAge;
        List<TimeTreeNode> nodes = buildShuffledLeaves(taxaNames, n);

        // index 0 is the root reference; the remaining n-1 are the between-tip depths: n-2 from
        // q/Q(rootAge), plus rootAge placed uniformly so the maximum depth equals the root age.
        List<Double> times = new ArrayList<>(n);
        times.add(rootAge);
        if (n > 2) {
            for (double age : sampleAges(0, rootAge, n - 2)) times.add(age);
        }
        times.add(random.nextInt(times.size()) + 1, rootAge);

        coalesce(nodes, times);
        TimeTree tree = new TimeTree();
        tree.setRoot(nodes.getFirst(), true);
        return new RandomVariable<>("CPPTree", tree, this);
    }

    /**
     * Fixed n, <em>random</em> stem: the stem is the max of n i.i.d. node ages — density
     * n·q(t)·Q(t)^(n-1), sampled via {@link #sampleStemAge} — and the n-1 internal depths are then
     * drawn from q/Q(stem) ({@link #sampleAges} truncated to (0, stem)). The root is the deepest of
     * those depths; the stem sits above it.
     */
    private RandomVariable<TimeTree> sampleRandomStemTree(String[] taxaNames, int n) {
        double stem = sampleStemAge(0, n);           // origin ~ n·q·Q^(n-1): max of n i.i.d. node ages
        conditionAge = stem;
        List<TimeTreeNode> nodes = buildShuffledLeaves(taxaNames, n);

        TimeTree tree = new TimeTree();
        if (n == 1) {
            tree.setRoot(nodes.getFirst(), true);
            tree.getRoot().setRootStem(stem);
            return new RandomVariable<>("CPPTree", tree, this);
        }

        // index 0 is the stem reference (>= every depth); the n-1 depths follow
        List<Double> times = new ArrayList<>(n);
        times.add(stem);
        for (double depth : sampleAges(0, stem, n - 1)) times.add(depth);
        coalesce(nodes, times);

        tree.setRoot(nodes.getFirst(), true);
        tree.getRoot().setRootStem(stem);
        return new RandomVariable<>("CPPTree", tree, this);
    }

    /**
     * Random number of tips, origin-conditioned. Between-tip depths are drawn i.i.d. from the law
     * over [0, &infin;) via {@link #sampleAges}; the first draw exceeding the origin caps the tree
     * and is discarded, so the tip count is one more than the number of sub-origin depths (geometric
     * with success probability Q(origin)). The root is the deepest sub-origin node; the origin sits
     * above it as a stem. Law-generic version of the old constant-rate {@code CPPTree} n==0 branch.
     */
    private RandomVariable<TimeTree> sampleRandomNTree(String[] taxaNames, Double originAge) {
        // origin: use the supplied age, else draw one unconditionally from the law
        double origin = (originAge != null) ? originAge : sampleAges(0, Double.POSITIVE_INFINITY, 1)[0];
        conditionAge = origin;

        // grow between-tip depths until one exceeds the origin (that draw caps the tree; discard it)
        List<Double> depths = new ArrayList<>();
        while (true) {
            double h = sampleAges(0, Double.POSITIVE_INFINITY, 1)[0];
            if (h > origin) break;
            depths.add(h);
            if (depths.size() > 1_000_000) {
                throw new RuntimeException("Random-N CPP did not terminate after 1e6 tips; the process " +
                        "is (near-)critical for these parameters. Provide n, or a smaller origin.");
            }
        }
        int nTips = depths.size() + 1;
        List<TimeTreeNode> nodes = buildShuffledLeaves(taxaNames, nTips);

        TimeTree tree = new TimeTree();
        if (nTips == 1) {                       // first draw exceeded the origin: a lone tip on a stem
            tree.setRoot(nodes.getFirst(), true);
            tree.getRoot().setRootStem(origin);
            return new RandomVariable<>("CPPTree", tree, this);
        }

        // coalesce: index 0 is the origin reference (>= every depth, so never the min); depths follow
        List<Double> times = new ArrayList<>(nTips);
        times.add(origin);
        times.addAll(depths);                   // origin + (nTips - 1) depths == nTips entries
        coalesce(nodes, times);

        tree.setRoot(nodes.getFirst(), true);
        tree.getRoot().setRootStem(origin);     // root age = deepest depth <= origin; stem up to origin
        return new RandomVariable<>("CPPTree", tree, this);
    }

    /** Construct a subclade generator of the same concrete type, reusing this generator's rate parameters. */
    protected abstract AbstractCalibratedCPPTree newSubClade(int nTaxa, CalibrationArray subCalibrations);

    @Override
    public RandomVariable<TimeTree> sample() {
        resolveRates();

        // if no calibrations, build an uncalibrated CPP. A null n means "grow a random number
        // of tips from the process" (origin-conditioned); a fixed n conditions on the root age.
        if (getCalibrations() == null) {
            String[] names = getOtherNames() == null ? null : getOtherNames().value();
            if (getN() == null) {                        // random number of tips, origin-conditioned
                return sampleRandomNTree(names,
                        getRootAge() == null ? null : getRootAge().value().doubleValue());
            }
            if (getRootAge() == null) {                  // fixed n, stem sampled from n·q·Q^(n-1)
                return sampleRandomStemTree(names, getN().value());
            }
            return sampleUncalibratedTree(names, getRootAge().value().doubleValue());  // fixed n, fixed root
        }

        // obtain pass in parameters
        double samplingProb = getSamplingProb().value().doubleValue();
        int n = getN().value();
        CalibrationArray calibrationArray = getCalibrations().value();
        Calibration[] calibrations = calibrationArray.getCalibrationArray();

        // initialise params
        double rootAgeValue = 0;
        TimeTree tree = new TimeTree();

        List<String> backUpNames = new ArrayList<>();

        // step1: get valid clade calibrations
        List<Calibration> cladeCalibrations = new ArrayList<>(Arrays.stream(calibrations).toList());

        // sort it with decreasing order
        cladeCalibrations.sort((c1, c2) -> Double.compare(c2.getAge(), c1.getAge()));

        // if root calibration is already in clade calibration
        if (cladeCalibrations.getFirst().getTaxa().length == n) {
            rootConditioned = true;
            rootAgeValue = cladeCalibrations.getFirst().getAge();
            // if only one root calibration, then return cpp
            if (cladeCalibrations.size() == 1) {
                return sampleUncalibratedTree(cladeCalibrations.getFirst().getTaxa(),
                        cladeCalibrations.getFirst().getAge());
            } else {
                // else remove the root calibration from cladeCalibrations
                backUpNames.addAll(Arrays.asList(cladeCalibrations.getFirst().getTaxa()));
                cladeCalibrations.remove(cladeCalibrations.getFirst());
            }
        }

        // step2: get all maximal calibration
        List<Calibration> maximalCalibrations = getMaximalCalibrations(cladeCalibrations);

        // map the taxa names for calibration clades
        int index = 0;
        int cladeSizes = 0;
        String[][] taxaNames = new String[maximalCalibrations.size()][];
        for (Calibration entry : maximalCalibrations) {
            taxaNames[index] = entry.getTaxa();
            cladeSizes += taxaNames[index].length;
            index++;
        }

        // calculate the number of nodes
        // m = non-clade tips + clade roots
        int m = n - cladeSizes + maximalCalibrations.size();

        /* initialise the lists
            A : holding the indices of inactive nodes (wait for assign)
            l : holding the indices of active nodes (has assigned)
            times : holding the times of internal nodes, the first element is root or stem age
            nodeList : holding all nodes
         */
        List<Integer> A = new ArrayList<>(m);
        for (int i = 0; i < m; i++) {
            A.add(i);
        }
        int[] l = new int[m];
        List<Double> times = new ArrayList<>(Collections.nCopies(m, 0.0));
        List<TimeTreeNode> nodeList = new ArrayList<>((Collections.nCopies(m, null)));

        // step3: calculate condition age (root or stem age)
        // if rootConditioned, then condition on root
        // if !rootConditioned, then use stem age or sample one
        if (rootConditioned) {
            int ind;
            ind = random.nextInt(m - 1) + 1; // [1, m-1]
            if (m == 2) {
                ind = 1;
            }

            times.set(ind, rootAgeValue);
            conditionAge = rootAgeValue;
        } else {
            if (getStemAge() != null) {
                conditionAge = getStemAge().value().doubleValue();
            } else {
                int idx = 0;
                while (Double.isNaN(conditionAge) || Double.isInfinite(conditionAge) || conditionAge == 0.0) {
                    conditionAge = sampleStemAge(maximalCalibrations.getFirst().getAge(), n);
                    idx++;
                    if (idx > 200) {
                        throw new RuntimeException("The stem age cannot be sampled because of the bad parameter combination. Please provide a stemAge.");
                    }
                }
            }
        }

        // step4: build clades for each maximalCalibrations
        nameList = new ArrayList<>(n);
        usedNames = new HashSet<>();

        for (Calibration maximalCalibration : maximalCalibrations) {
            String[] currentNames = maximalCalibration.getTaxa();
            String[] uniqueNames = new String[currentNames.length];
            for (int j = 0; j < currentNames.length; j++) {
                String newName = makeUnique(currentNames[j], usedNames);

                nameList.add(newName);
                uniqueNames[j] = newName;
                backUpNames.remove(currentNames[j]);
            }
            maximalCalibration.setTaxa(uniqueNames);
        }

        // loop through all maximalCalibrations
        for (int i = 0; i < maximalCalibrations.size(); i++) {
            // step1: get subclades
            List<Calibration> subClades = getNestedClades(maximalCalibrations.get(i), cladeCalibrations);
            // step2: get sampled element
            // calculate weights
            double w = cdf(conditionAge) - cdf(maximalCalibrations.get(i).getAge());
            // calculate score s for each node
            int[] s = calculateScore(A, m, times);
            // calculate weight for each node
            double[] weights = getWeights(s, w);

            if (A.size() == 1) {
                l[i] = A.getFirst();
            } else {
                // sample one element from A with probability weights
                l[i] = A.get(sampleIndex(weights));
            }

            // step3: construct subtrees
            double cladeAge = maximalCalibrations.get(i).getAge();
            String[] subcladeTaxa = maximalCalibrations.get(i).getTaxa();

            Calibration[] clades = new Calibration[subClades.size() + 1];
            clades[0] = maximalCalibrations.get(i);
            for (int j = 0; j < subClades.size(); j++) {
                Calibration cal = new Calibration(subClades.get(j).getTaxa());
                cal.setAge(subClades.get(j).getAge());
                clades[j + 1] = cal;
            }

            // simulate a tree for these clades, only offer calibrations
            AbstractCalibratedCPPTree subTreeGen = newSubClade(subcladeTaxa.length, new CalibrationArray(clades));

            // put clade into nodeList
            TimeTree subTree = subTreeGen.sample().value();
            nodeList.set(l[i], subTree.getRoot());

            // sample times at l[i] and l[i]+1 conditioned to be older than the clade age,
            if (l[i] > 0 && times.get(l[i]) == 0) {
                times.set(l[i], sampleAges(cladeAge, conditionAge, 1)[0]);
            }

            if (l[i] < m - 1 && times.get(l[i] + 1) == 0) {
                times.set(l[i] + 1, sampleAges(cladeAge, conditionAge, 1)[0]);
            }

            // remove corresponding node in A
            A.remove(Integer.valueOf(l[i]));
        }

        // step5: organise times
        // after calibrations, sample times for remaining unassigned nodes
        List<Integer> zeroIndices = new ArrayList<>();
        for (int i = 1; i < times.size(); i++) {
            if (times.get(i) == 0) zeroIndices.add(i);
        }
        double[] sampledTimes = sampleAges(0, conditionAge, zeroIndices.size());
        for (int j = 0; j < zeroIndices.size(); j++) {
            times.set(zeroIndices.get(j), sampledTimes[j]);
        }

        // set the first node to be the max, make it the root
        if (times.size() > 2) {
            double max = times.get(1);
            for (int i = 1; i < times.size(); i++) {
                if (times.get(i) > max) {
                    max = times.get(i);
                }
            }
            times.set(0, max); // set this the largest
        } else if (times.size() == 2) {
            times.set(0, times.get(1));
        } else {
            throw new RuntimeException("Unreachable");
        }

        if (rootConditioned) {
            if (Math.abs(times.getFirst() - rootAgeValue) > 1e-8) {
                // shouldn't have this thrown theoretically
                throw new RuntimeException("The max age is not root age when root conditioned");
            }
        }

        // or if not root conditioned, set it to stem age
        if (!rootConditioned) {
            times.set(0, conditionAge);
        }

        // step6: fill in nodelist
        // get non-clade taxa
        List<String> outGroupTaxa = new ArrayList<>();
        if (getOtherNames() != null) {
            String[] otherNamesArr = getOtherNames().value();
            for (String otherName : otherNamesArr) {
                String newName = makeUnique(otherName, usedNames);
                nameList.add(newName);
                outGroupTaxa.add(newName);
            }
        } else {
            if (!backUpNames.isEmpty()) {
                for (String name : backUpNames) {
                    nameList.add(name);
                    outGroupTaxa.add(name);
                }
            }
        }

        // fit other names in if there are non-assigned names
        int nameListSize = nameList.size();
        for (int i = 0; i < n - nameListSize; i++) {
            String indexName = String.valueOf(i);
            String newName = makeUnique(indexName, usedNames);
            nameList.add(newName);
            outGroupTaxa.add(newName);
        }

        // get random order for non-clade taxa
        Collections.shuffle(outGroupTaxa);

        // Assign remaining uncalibrated taxa names to available node positions
        int ind = 0;
        for (int i = 0; i < nodeList.size() && ind < outGroupTaxa.size(); i++) {
            if (nodeList.get(i) == null) {
                TimeTreeNode tip = new TimeTreeNode(0);
                tip.setId(outGroupTaxa.get(ind));
                nodeList.set(i, tip);
                ind++;
            }
        }

        // step7: coalesce
        // combine sub-CPPs into the final tree
        coalesce(nodeList, times);

        tree.setRoot(nodeList.getFirst(), true);
        if (!rootConditioned && conditionAge != nodeList.getFirst().getAge()) {
            tree.getRoot().setRootStem(times.getFirst());
        }

        return new RandomVariable<>("CPPTree", tree, this);
    }

    // public for unit test
    public static void coalesce(List<TimeTreeNode> nodeList, List<Double> times) {
        while (nodeList.size() > 1) {
            // start from the youngest node
            int j = indexOfMin(times);
            if (times.size() == 2) {
                j = 1; // make it being the second one if there are only 2 nodes
            }

            // build relationship
            TimeTreeNode child_left = nodeList.get(j - 1);
            TimeTreeNode child_right = nodeList.get(j);

            TimeTreeNode parent = new TimeTreeNode(times.get(j));
            parent.addChild(child_left);
            parent.addChild(child_right);

            child_left.setParent(parent);
            child_right.setParent(parent);

            // adjust indices
            nodeList.set(j - 1, parent);
            nodeList.remove(j);

            // remove the time and age of the second node
            times.remove(j);
        }
    }

    /*
       Functions
    */
    private static int[] calculateScore(List<Integer> A, int m, List<Double> times) {
        int[] s = new int[A.size()];
        for (int j = 0; j < A.size(); j++) {
            int nodeIndex = A.get(j);
            int count = 0;
            // Check if i < m and i+1 is within bounds
            if (nodeIndex < m - 1 && times.get(nodeIndex + 1) == 0) count++;
            // Check if nodeIndex is within bounds and times[i] == 0
            if (nodeIndex > 0 && times.get(nodeIndex) == 0) count++;
            s[j] = count;
        }
        return s;
    }

    private static double[] getWeights(int[] s, double w) {
        double sumOfWeights = 0;
        double[] weights = new double[s.length];

        for (int i = 0; i < s.length; i++) {
            weights[i] = Math.pow(w, s[i]);
            sumOfWeights += weights[i];
        }

        // normalise weights
        for (int i = 0; i < s.length; i++) {
            weights[i] /= sumOfWeights;
        }

        return weights;
    }

    @Override
    public Map<String, Value> getParams() {
        Map<String, Value> map = super.getParams();
        map.put(BirthDeathConstants.rhoParamName, rho);
        map.put(DistributionConstants.nParamName, n);
        if (calibrations != null) map.put(calibrationsName, calibrations);
        if (rootAge != null) map.put(rootAgeName, rootAge);
        if (stemAge != null) map.put(stemAgeName, stemAge);
        if (otherNames != null) map.put(otherTaxaNames, otherNames);
        return map;
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(BirthDeathConstants.rhoParamName)) rho = value;
        else if (paramName.equals(DistributionConstants.nParamName)) n = value;
        else if (paramName.equals(calibrationsName)) calibrations = value;
        else if (paramName.equals(otherTaxaNames)) otherNames = value;
        else if (paramName.equals(stemAgeName)) stemAge = value;
        else if (paramName.equals(rootAgeName)) rootAge = value;
        else throw new IllegalArgumentException("Unknown parameter name: " + paramName);
    }

    public Value<Integer> getN() {
        return getParams().get(DistributionConstants.nParamName);
    }

    public Value<Number> getSamplingProb() {
        return getParams().get(BirthDeathConstants.rhoParamName);
    }

    public Value<CalibrationArray> getCalibrations() {
        return getParams().get(calibrationsName);
    }

    public Value<Number> getStemAge() {
        return getParams().get(stemAgeName);
    }

    public Value<Number> getRootAge() {
        return getParams().get(rootAgeName);
    }

    public Value<String[]> getOtherNames() {
        return getParams().get(otherTaxaNames);
    }

    public Value<Double> getOrigin() {
        return new Value<>("", conditionAge);
    }

    public boolean getRootCondition() {
        return rootConditioned;
    }
}
