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
            double t = invertCDF(qLo + random.nextDouble() * (qHi - qLo));
            // Clamp to [lowerTime, upperTime]: near CDF saturation the numerical inverse can overshoot
            // the nearly-flat upper bound, which would otherwise let a depth exceed the origin.
            out[i] = t < lowerTime ? lowerTime : (t > upperTime ? upperTime : t);
        }
        return out;
    }

    /**
     * Sample a stem/origin age older than greaterThan for nTaxa taxa. The origin has
     * CDF Q(t)^nTaxa (the max of nTaxa i.i.d. node ages), left-truncated at greaterThan.
     * Derived from {@link #cdf}; constant-rate subclasses override with a closed form.
     *
     * <p>Only well defined when Q saturates (lim Q(t) = 1, a supercritical process). For a critical or
     * subcritical process lim Q(t) &lt; 1, so a node age — and hence the max of n — is infinite with
     * positive probability and no finite origin exists; the inverse CDF returns +&infin; there. Rather
     * than let that propagate to NaN branch lengths, fail with guidance to condition on an explicit origin.
     */
    protected double sampleStemAge(double greaterThan, int nTaxa) {
        double qaN = Math.pow(cdf(greaterThan), nTaxa);
        double p = qaN + random.nextDouble() * (1.0 - qaN);
        double stem = invertCDF(Math.pow(p, 1.0 / nTaxa));
        if (!Double.isFinite(stem)) {
            throw new RuntimeException("Cannot sample a finite origin: the node-age CDF does not saturate "
                    + "(lim Q(t) < 1 — a critical or subcritical process, where lambda <= mu). Condition the "
                    + "tree by supplying stemAge or rootAge.");
        }
        return stem;
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
        Set<String> used = new HashSet<>();
        List<String> names = new ArrayList<>(nTips);
        if (taxaNames != null) {
            for (int i = 0; i < Math.min(taxaNames.length, nTips); i++) names.add(makeUnique(taxaNames[i], used));
        }
        for (int i = names.size(); i < nTips; i++) names.add(makeUnique(String.valueOf(i), used));
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
     * Assemble a CPP tree from an origin age and a list of between-tip depths (each &lt;= origin).
     * Leaves are coalesced by repeated min-depth merge, so the root is the deepest depth. When
     * {@code stem} is true the origin sits above the root as a stem; otherwise the origin coincides
     * with the maximum depth (root-conditioned, no stem). Shared by all three uncalibrated cases —
     * they differ only in how (origin, depths) are produced.
     */
    private RandomVariable<TimeTree> buildFromDepths(String[] taxaNames, double origin, List<Double> depths, boolean stem) {
        int nTips = depths.size() + 1;
        conditionAge = origin;
        List<TimeTreeNode> nodes = buildShuffledLeaves(taxaNames, nTips);

        TimeTree tree = new TimeTree();
        if (nTips > 1) {
            // index 0 is the origin reference (>= every depth, so never the min); the depths follow
            List<Double> times = new ArrayList<>(nTips);
            times.add(origin);
            times.addAll(depths);
            coalesce(nodes, times);
        }
        tree.setRoot(nodes.getFirst(), true);
        if (stem) tree.getRoot().setRootStem(origin);
        return new RandomVariable<>("CPPTree", tree, this);
    }

    /**
     * Fixed n, fixed root: n-2 internal depths from q/Q(rootAge) plus the root age placed uniformly
     * among them (so the root is the maximum depth). Root-conditioned, no stem.
     */
    private RandomVariable<TimeTree> sampleFixedRootTree(String[] taxaNames, int n, double rootAge) {
        List<Double> depths = new ArrayList<>(Math.max(0, n - 1));
        if (n >= 2) {
            if (n > 2) for (double age : sampleAges(0, rootAge, n - 2)) depths.add(age);
            depths.add(random.nextInt(depths.size() + 1), rootAge);   // place the root uniformly
        }
        return buildFromDepths(taxaNames, rootAge, depths, false);
    }

    /**
     * Fixed n, <em>random</em> stem: the stem is the max of n i.i.d. node ages — density
     * n·q(t)·Q(t)^(n-1), sampled via {@link #sampleStemAge} — and the n-1 internal depths are drawn
     * from q/Q(stem). The root is the deepest of those depths; the stem sits above it.
     */
    private RandomVariable<TimeTree> sampleRandomStemTree(String[] taxaNames, int n) {
        double stem = sampleStemAge(0, n);           // origin ~ n·q·Q^(n-1): max of n i.i.d. node ages
        List<Double> depths = new ArrayList<>(Math.max(0, n - 1));
        for (double depth : sampleAges(0, stem, n - 1)) depths.add(depth);
        return buildFromDepths(taxaNames, stem, depths, true);
    }

    /**
     * Fixed n, fixed stem: condition on the supplied stem age rather than sampling it. The n-1 internal
     * depths are drawn from the law truncated to [0, stem]; the stem sits above the root. Unlike the
     * sampled-stem case this is well defined for a subcritical process, whose max-of-n origin is infinite
     * with positive probability (Q saturates below 1) and would otherwise yield an infinite root age.
     */
    private RandomVariable<TimeTree> sampleFixedStemTree(String[] taxaNames, int n, double stem) {
        List<Double> depths = new ArrayList<>(Math.max(0, n - 1));
        for (double depth : sampleAges(0, stem, n - 1)) depths.add(depth);
        return buildFromDepths(taxaNames, stem, depths, true);
    }

    /**
     * Fixed stem, random number of tips. Between-tip depths are drawn i.i.d. from the law over
     * [0, &infin;) via {@link #sampleAges}; the first draw exceeding the stem caps the tree and is
     * discarded, so the tip count is one more than the number of sub-stem depths (geometric with success
     * probability 1 - Q(stem)). No node sits at the stem; it forms a stem above the root. Reached only
     * when no taxa were supplied, so names are always indices.
     */
    private RandomVariable<TimeTree> sampleRandomNStemTree(String[] taxaNames, double stem) {
        List<Double> depths = new ArrayList<>();
        while (true) {
            double h = sampleAges(0, Double.POSITIVE_INFINITY, 1)[0];
            if (h > stem) break;                      // this draw caps the tree; discard it
            depths.add(h);
            checkRunawayTipCount(depths.size());
        }
        return buildFromDepths(taxaNames, stem, depths, true);
    }

    /**
     * Fixed root, random number of tips. The root is the crown of two independent sub-clades, so depths
     * accumulate through two geometric runs: the first age exceeding the root is where the root itself is
     * placed (becoming the deepest depth), and the second age exceeding the root ends the tree and is
     * discarded. Root-conditioned, no stem.
     */
    private RandomVariable<TimeTree> sampleRandomNRootTree(String[] taxaNames, double root) {
        List<Double> depths = new ArrayList<>();
        int exceedances = 0;
        while (exceedances < 2) {
            double h = sampleAges(0, Double.POSITIVE_INFINITY, 1)[0];
            if (h > root) {
                if (++exceedances == 1) depths.add(root);   // place the root at the first exceedance
            } else {
                depths.add(h);
            }
            checkRunawayTipCount(depths.size());
        }
        return buildFromDepths(taxaNames, root, depths, false);
    }

    private static void checkRunawayTipCount(int size) {
        if (size > 1_000_000) {
            throw new RuntimeException("Random-N CPP did not terminate after 1e6 tips; the process " +
                    "is (near-)critical for these parameters. Provide n, or a smaller origin.");
        }
    }

    /** Construct a subclade generator of the same concrete type, reusing this generator's rate parameters. */
    protected abstract AbstractCalibratedCPPTree newSubClade(int nTaxa, CalibrationArray subCalibrations);

    @Override
    public RandomVariable<TimeTree> sample() {
        resolveRates();

        // if no calibrations, build an uncalibrated CPP. The tip count is, in priority order, the
        // explicit n, else the number of supplied taxa, else random (grown from the process). A fixed
        // count with no rootAge samples the stem; with a rootAge it conditions on that root.
        if (getCalibrations() == null) {
            String[] names = getOtherNames() == null ? null : getOtherNames().value();
            Integer nTips = (getN() != null) ? getN().value() : (names != null ? names.length : null);
            if (names != null && nTips != null && names.length > nTips) {
                throw new IllegalArgumentException("taxa has " + names.length + " names but n = " + nTips
                        + "; taxa must contain at most n elements.");
            }
            if (nTips == null) {                         // random number of tips (no n, no taxa)
                // A fixed STEM caps a single geometric run: the first node age above the stem ends the tree
                // and is discarded (no node sits at the stem; it forms a stem above the root). A fixed ROOT
                // is the crown of two independent sub-clades: the first age above the root is where the root
                // is placed, then a second run accumulates until the second age above the root ends the tree.
                // With neither, the origin is sampled from the process and treated as a stem.
                if (getStemAge() != null) {
                    return sampleRandomNStemTree(names, getStemAge().value().doubleValue());
                } else if (getRootAge() != null) {
                    return sampleRandomNRootTree(names, getRootAge().value().doubleValue());
                }
                return sampleRandomNStemTree(names, sampleAges(0, Double.POSITIVE_INFINITY, 1)[0]);
            }
            if (getRootAge() != null) {                  // fixed n, fixed root
                return sampleFixedRootTree(names, nTips, getRootAge().value().doubleValue());
            }
            if (getStemAge() != null) {                  // fixed n, fixed stem (conditions on the origin)
                return sampleFixedStemTree(names, nTips, getStemAge().value().doubleValue());
            }
            return sampleRandomStemTree(names, nTips);   // fixed n, stem sampled from n·q·Q^(n-1)
        }

        // obtain pass in parameters
        if (getN() == null) {
            throw new IllegalArgumentException("n must be provided when calibrations are given "
                    + "(a random number of tips is only supported for uncalibrated trees).");
        }
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
                return sampleFixedRootTree(cladeCalibrations.getFirst().getTaxa(), n,
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
            // Youngest merge point among the between-tip depths (indices >= 1). Index 0 is the
            // origin/stem reference and is never a merge point; searching from 1 (rather than a
            // global argmin) keeps that invariant even in degenerate regimes where the origin is
            // tied with, or numerically below, the depths.
            int j = 1;
            double min = times.get(1);
            for (int i = 2; i < times.size(); i++) {
                if (times.get(i) < min) {
                    min = times.get(i);
                    j = i;
                }
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
        // n is optional (omit -> random number of tips): keep it out of the params graph when null,
        // otherwise LPhy traverses a null node. Also drop any null n that super.getParams() added.
        if (n != null) map.put(DistributionConstants.nParamName, n);
        else map.remove(DistributionConstants.nParamName);
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
