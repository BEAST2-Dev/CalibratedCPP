package calibrationprior;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import calibration.CalibrationForest;
import calibration.CalibrationNode;
import org.hipparchus.linear.*;
import org.hipparchus.special.Erf;
import org.hipparchus.special.Gamma;

import java.io.PrintStream;
import java.util.*;
import java.util.Collections;

/**
 * @author Marcus Overwater
 */

@Description("A topology consistent prior on clade calibration ages.")
public class CalibrationPrior extends Distribution {

    public Input<TreeInterface> treeInput =
            new Input<>("tree", "Tree to calibrate", Input.Validate.REQUIRED);

    public Input<List<CalibrationCladePrior>> cladesInput =
            new Input<>("calibration", "List of calibration clades (alternative to calibrationForest).", new ArrayList<>());

    public Input<CalibrationForest> calibrationForestInput =
            new Input<>("calibrationForest",
                    "Constraint tree with [&lower=X,upper=Y] annotations. "
                    + "Alternative to providing an explicit list of calibration clades.");

    private List<CalibrationNode> calibrationNodes = new ArrayList<>();

    @Override
    public void initAndValidate() {
        TreeInterface tree = treeInput.get();
        if (tree == null) throw new IllegalArgumentException("Tree is null");

        List<CalibrationCladePrior> explicitClades = cladesInput.get();
        CalibrationForest forest = calibrationForestInput.get();

        if (!explicitClades.isEmpty() && forest != null) {
            throw new IllegalArgumentException(
                    "Provide calibration clades via 'calibration' or 'calibrationForest', not both.");
        }

        List<CalibrationCladePrior> cladePriors;
        if (forest != null) {
            // A CalibrationForest is the shared currency; CalibrationPrior needs the bounded
            // clades. A forest built from bare taxon sets carries no bounds — flag it rather
            // than silently contributing nothing to the density.
            cladePriors = forest.getCalibrationCladePriors();
            if (cladePriors.isEmpty() && !forest.getAllNodes().isEmpty()) {
                throw new IllegalArgumentException(
                        "CalibrationPrior requires clades with age bounds, but the supplied "
                        + "'calibrationForest' contains only bound-free clades. Add "
                        + "[&lower=X,upper=Y] annotations (or lowerAge/upperAge) to the "
                        + "calibrated clades.");
            }
        } else {
            cladePriors = explicitClades;
        }

        if (cladePriors.isEmpty()) {
            // No calibrations — nothing to fit; calculateLogP() will return 0
            calibrationNodes = Collections.emptyList();
            return;
        }

        // === Step 2: build inclusion hierarchy ===
        CalibrationForest calibrationForest = CalibrationForest.fromPriors(cladePriors);
        calibrationNodes = calibrationForest.getAllNodes();

        for (CalibrationNode node : calibrationNodes) {
            if (node.parent != null) {
                node.getCalibrationCladePrior().isOverlapEdge = node.parent.getCalibrationCladePrior().getLower() < node.getCalibrationCladePrior().getUpper();
            }
        }

        // === Step 3: compute log-moments ===
        for (CalibrationNode n : calibrationNodes) computeLogTargets(n);

        // === Step 3b: validate overlap edges ===
        for (CalibrationNode node : calibrationNodes) {
            if (node.parent == null || !node.getCalibrationCladePrior().isOverlapEdge) continue;

            CalibrationCladePrior child = node.getCalibrationCladePrior();
            CalibrationCladePrior parent = node.parent.getCalibrationCladePrior();

            if (child.sigma2 < parent.sigma2) {
                double childLogWidth  = Math.log(child.getUpper())  - Math.log(child.getLower());
                double parentLogWidth = Math.log(parent.getUpper()) - Math.log(parent.getLower());
                Log.warning(String.format(
                    "WARNING: Calibration clade '%s' (bounds [%.4g, %.4g], log-width %.4g) is more " +
                    "precisely specified than its overlapping ancestor clade '%s' " +
                    "(bounds [%.4g, %.4g], log-width %.4g). " +
                    "For overlapping calibrations the child's calibration interval should span " +
                    "at least as wide a range on the log scale as the parent's. " +
                    "The Beta prior edge variance will be clamped to a small positive value.",
                    node.prior.getID(),
                    child.getLower(), child.getUpper(), childLogWidth,
                    node.parent.prior.getID(),
                    parent.getLower(), parent.getUpper(), parentLogWidth));
            }
        }

        // === Step 4: partition into overlap components and fit ===
        List<CalibrationComponent> comps = CalibrationComponent.partition(calibrationForest);
        for (CalibrationComponent comp : comps) fitComponent(comp);
        // done
    }

    private void fitComponent(CalibrationComponent comp) {
        CalibrationNode root = comp.getRoot();
        List<CalibrationNode> members = comp.getMembers();

        // Map each node to a column index in the least-squares system
        Map<CalibrationNode, Integer> colIndex = new LinkedHashMap<>();
        int idx = 0;
        for (CalibrationNode n : members) {
            if (n != root) colIndex.put(n, idx++);
        }

        int K = colIndex.size();
        if (K == 0) return; // only root in component, nothing to fit

        RealMatrix A = new Array2DRowRealMatrix(K, K);
        double[] bMu = new double[K];
        double[] bVar = new double[K];

        // Build linear system: log(child) - log(root) = sum of edges along path
        for (CalibrationNode child : colIndex.keySet()) {
            int row = colIndex.get(child);

            List<CalibrationNode> path = pathFromRoot(root, child);
            for (CalibrationNode step : path) {
                if (step == root) continue;
                Integer col = colIndex.get(step);
                if (col != null) {
                    A.addToEntry(row, col, 1.0);
                }
            }

            bMu[row] = child.getCalibrationCladePrior().mu - root.getCalibrationCladePrior().mu;
            bVar[row] = child.getCalibrationCladePrior().sigma2 - root.getCalibrationCladePrior().sigma2;
        }

        // Solve least-squares for mean and variance edges
        RealVector mHat = solveQR(A, new ArrayRealVector(bMu));
        RealVector vHat = nnls(A, new ArrayRealVector(bVar));
        for (int i = 0; i < K; i++) {
            if (vHat.getEntry(i) < 1e-8) vHat.setEntry(i, 1e-8);
        }

        for (int i = 0; i < K; i++) {
            if (vHat.getEntry(i) <= 0) vHat.setEntry(i, 1e-8);
        }

        // Assign edges and Beta parameters
        for (Map.Entry<CalibrationNode, Integer> e : colIndex.entrySet()) {
            CalibrationNode child = e.getKey();
            int j = e.getValue();

            child.getCalibrationCladePrior().mEdge = mHat.getEntry(j);
            child.getCalibrationCladePrior().vEdge = vHat.getEntry(j);

            double[] ab = invertLogMomentsToBetaParams(child.getCalibrationCladePrior().mEdge, child.getCalibrationCladePrior().vEdge);
            child.getCalibrationCladePrior().alpha = ab[0];
            child.getCalibrationCladePrior().beta = ab[1];
        }
    }

    // ------------------------------------------------------------------
    // 3. Fit each component

    private RealVector solveQR(RealMatrix A, RealVector b) {
        try {
            return new QRDecomposition(A).getSolver().solve(b);
        } catch (Exception e) {
            return new SingularValueDecomposition(A).getSolver().solve(b);
        }
    }

    private List<CalibrationNode> pathFromRoot(CalibrationNode root, CalibrationNode node) {
        LinkedList<CalibrationNode> path = new LinkedList<>();
        CalibrationNode cur = node;
        while (cur != null && cur != root) {
            path.addFirst(cur);
            cur = cur.parent;
        }
        path.addFirst(root);
        return path;
    }

    // ------------------------------------------------------------------
    // log-moment target from calibration bounds
    private void computeLogTargets(CalibrationNode n) {
        double tLo = n.getCalibrationCladePrior().getLower();
        double tHi = n.getCalibrationCladePrior().getUpper();
        double p = n.getCalibrationCladePrior().getCoverage();

        // Replaced NormalDistribution.inverseCumulativeProbability with Erf.erfInv
        // This calculates the z-score for the given coverage probability
        double z = Math.sqrt(2.0) * Erf.erfInv(p);

        double sigma = (Math.log(tHi) - Math.log(tLo)) / (2 * z);
        n.getCalibrationCladePrior().sigma2 = sigma * sigma;
        n.getCalibrationCladePrior().mu = Math.log(tLo) + z * sigma;
    }

    // Minimum total concentration (alpha+beta) expressed as a multiple of 1/mean(1-mean).
    // A value of 1 ensures each Beta has at least one "pseudo-observation" relative to its mean,
    // preventing near-degenerate point-mass or near-uniform distributions.
    private static final double MIN_CONCENTRATION_FACTOR = 1.0;

    // ------------------------------------------------------------------
    // Utility math
    public static double[] invertLogMomentsToBetaParams(double m, double v) {
        // Asymptotic initialization: for large n=a+b, digamma(a)-digamma(n) ≈ ln(a/n)
        // and trigamma(a)-trigamma(n) ≈ 1/a - 1/n = (1-mu)/(mu*n).
        // Solving these gives mu0 = exp(m) and n0 = (1-mu0)/(mu0*v).
        double mu0 = Math.exp(m);
        mu0 = Math.min(1.0 - 1e-6, Math.max(1e-6, mu0));
        double n0 = Math.max(1e-3, (1.0 - mu0) / (mu0 * Math.max(v, 1e-15)));
        double a = Math.max(1e-6, mu0 * n0);
        double b = Math.max(1e-6, n0 - a);

        // 2×2 Newton on F(a,b) = [digamma(a)-digamma(a+b)-m, trigamma(a)-trigamma(a+b)-v]
        // Jacobian: J = [[triA-triN, -triN], [tetA-tetN, -tetN]]
        // where tetragamma is approximated via central finite difference of trigamma.
        for (int iter = 0; iter < 200; iter++) {
            double n = a + b;
            double f1 = Gamma.digamma(a) - Gamma.digamma(n) - m;
            double f2 = Gamma.trigamma(a) - Gamma.trigamma(n) - v;
            if (Math.abs(f1) < 1e-10 && Math.abs(f2) < 1e-10) break;

            double triA = Gamma.trigamma(a);
            double triN = Gamma.trigamma(n);
            double hA = Math.max(a * 1e-5, 1e-9);
            double hN = Math.max(n * 1e-5, 1e-9);
            double tetA = (Gamma.trigamma(a + hA) - Gamma.trigamma(a - hA)) / (2 * hA);
            double tetN = (Gamma.trigamma(n + hN) - Gamma.trigamma(n - hN)) / (2 * hN);

            double j11 = triA - triN;
            double j12 = -triN;
            double j21 = tetA - tetN;
            double j22 = -tetN;
            double det = j11 * j22 - j12 * j21;
            if (Math.abs(det) < 1e-20) break;

            // J^{-1} * [-f1, -f2]: da = (-j22*f1 + j12*f2)/det, db = (j21*f1 - j11*f2)/det
            double da = (-j22 * f1 + j12 * f2) / det;
            double db = ( j21 * f1 - j11 * f2) / det;

            // Armijo line search: halve step until ||F||^2 decreases
            double res0 = f1 * f1 + f2 * f2;
            double step = 1.0;
            for (int ls = 0; ls < 50; ls++) {
                double an = Math.max(a + step * da, 1e-6);
                double bn = Math.max(b + step * db, 1e-6);
                double nn = an + bn;
                double g1 = Gamma.digamma(an) - Gamma.digamma(nn) - m;
                double g2 = Gamma.trigamma(an) - Gamma.trigamma(nn) - v;
                if (g1 * g1 + g2 * g2 < res0) break;
                step *= 0.5;
            }
            a = Math.max(a + step * da, 1e-6);
            b = Math.max(b + step * db, 1e-6);
        }

        // Enforce minimum concentration relative to the mean.
        // Floor: n >= MIN_CONCENTRATION_FACTOR / (mu*(1-mu)) ensures each Beta has at least
        // MIN_CONCENTRATION_FACTOR effective pseudo-observations scaled to the mean.
        double mu = a / (a + b);
        double minN = MIN_CONCENTRATION_FACTOR / (mu * (1.0 - mu));
        if (a + b < minN) {
            a = minN * mu;
            b = minN * (1.0 - mu);
        }

        return new double[]{a, b};
    }

    private RealVector nnls(RealMatrix A, RealVector b) {
        int n = A.getColumnDimension();
        double[] x = new double[n];
        boolean[] passive = new boolean[n]; // true = in passive set P

        for (int outerIter = 0; outerIter < 3 * n; outerIter++) {
            // Gradient of residual: w = Aᵀ(b - Ax)
            RealVector residual = b.subtract(A.operate(new ArrayRealVector(x)));
            RealVector w = A.transpose().operate(residual);

            // Find best candidate in zero set Z
            int jStar = -1;
            double maxW = 1e-10; // tolerance
            for (int j = 0; j < n; j++) {
                if (!passive[j] && w.getEntry(j) > maxW) {
                    maxW = w.getEntry(j);
                    jStar = j;
                }
            }
            if (jStar == -1) break; // optimality condition met

            passive[jStar] = true;

            // Inner loop: solve on passive set, enforce positivity
            for (int innerIter = 0; innerIter < 3 * n; innerIter++) {
                // Build sub-matrix A_P from passive columns
                int pSize = 0;
                for (boolean p : passive) if (p) pSize++;
                int[] pIdx = new int[pSize];
                int k = 0;
                for (int j = 0; j < n; j++) if (passive[j]) pIdx[k++] = j;

                RealMatrix Ap = new Array2DRowRealMatrix(A.getRowDimension(), pSize);
                for (int c = 0; c < pSize; c++)
                    Ap.setColumnVector(c, A.getColumnVector(pIdx[c]));

                // Unconstrained LS solve on passive set
                RealVector sp;
                try {
                    sp = new QRDecomposition(Ap).getSolver().solve(b);
                } catch (Exception e) {
                    sp = new SingularValueDecomposition(Ap).getSolver().solve(b);
                }

                // Check if all passive entries are positive
                boolean allPositive = true;
                for (int c = 0; c < pSize; c++) {
                    if (sp.getEntry(c) <= 0) { allPositive = false; break; }
                }
                if (allPositive) {
                    for (int c = 0; c < pSize; c++) x[pIdx[c]] = sp.getEntry(c);
                    break;
                }

                // Step to boundary: α = min_{j∈P, s_j≤0} x_j / (x_j - s_j)
                double alpha = Double.POSITIVE_INFINITY;
                for (int c = 0; c < pSize; c++) {
                    if (sp.getEntry(c) <= 0) {
                        double xj = x[pIdx[c]];
                        alpha = Math.min(alpha, xj / (xj - sp.getEntry(c)));
                    }
                }

                // Update x and move near-zero passive entries back to Z
                for (int c = 0; c < pSize; c++) {
                    x[pIdx[c]] += alpha * (sp.getEntry(c) - x[pIdx[c]]);
                }
                for (int c = 0; c < pSize; c++) {
                    if (Math.abs(x[pIdx[c]]) < 1e-12) {
                        x[pIdx[c]] = 0.0;
                        passive[pIdx[c]] = false;
                    }
                }
            }
        }
        return new ArrayRealVector(x);
    }

    private void collectLeafTaxa(Node node, Set<String> leafIDs) {
        if (node.isLeaf()) {
            leafIDs.add(node.getID());
        } else {
            for (int i = 0; i < node.getChildCount(); i++) {
                collectLeafTaxa(node.getChild(i), leafIDs);
            }
        }
    }

    // ------------------------------------------------------------------
    // BEAST interface
    @Override
    public double calculateLogP() {
        TreeInterface tree = treeInput.get();
        double logP = 0;
        for (CalibrationNode n : calibrationNodes) {
            double lp = cladeLogDensity(n, tree);
            if (lp == Double.NEGATIVE_INFINITY) return this.logP = Double.NEGATIVE_INFINITY;
            logP += lp;
        }
        return this.logP = logP;
    }

    private double cladeLogDensity(CalibrationNode n, TreeInterface tree) {
        // Check for monophyly
        Set<String> leafIDs = new HashSet<>();
        Node mrca = n.getCommonAncestor(tree);

        collectLeafTaxa(mrca, leafIDs);
        if (!leafIDs.equals(n.taxa.getTaxaNames())) {
            return Double.NEGATIVE_INFINITY; // clade is not monophyletic!
        }

        double t = mrca.getHeight();
        if (n.parent == null) {
            // root lognormal
            return logNormalLogPdf(t, n.getCalibrationCladePrior().mu, Math.sqrt(n.getCalibrationCladePrior().sigma2));
        }

        Node pNode = n.parent.getCommonAncestor(tree);
        double tp = pNode.getHeight();
        if (n.getCalibrationCladePrior().isOverlapEdge) {
            // overlapping → Beta
            double r = t / tp;
            if (r <= 0 || r >= 1) return Double.NEGATIVE_INFINITY;
            return betaLogPdf(r, n.getCalibrationCladePrior().alpha, n.getCalibrationCladePrior().beta) - Math.log(tp); // Jacobian for conversion of the joint density of tp and r to tp and t
        } else {
            // non-overlapping → truncated lognormal
            double lp = logNormalLogPdf(t, n.getCalibrationCladePrior().mu, Math.sqrt(n.getCalibrationCladePrior().sigma2));
            double lcdf = logNormalLogCdf(tp, n.getCalibrationCladePrior().mu, Math.sqrt(n.getCalibrationCladePrior().sigma2));
            if (Double.isInfinite(lcdf)) return Double.NEGATIVE_INFINITY;
            return lp - lcdf;
        }
    }

    private double logNormalLogPdf(double x, double mu, double sigma) {
        if (x <= 0) return Double.NEGATIVE_INFINITY;
        double z = (Math.log(x) - mu) / sigma;
        return -Math.log(x * sigma * Math.sqrt(2 * Math.PI)) - 0.5 * z * z;
    }

    private double logNormalLogCdf(double x, double mu, double sigma) {
        if (x <= 0) return Double.NEGATIVE_INFINITY;

        // 1. Convert to Z-score
        double z = (Math.log(x) - mu) / sigma;

        // 2. Calculate Cumulative Probability using Erf
        // Formula: 0.5 * (1 + erf(z / sqrt(2)))
        double p = 0.5 * (1.0 + Erf.erf(z / Math.sqrt(2.0)));

        return Math.log(p);
    }

    private double betaLogPdf(double x, double a, double b) {
        if (x <= 0 || x >= 1) return Double.NEGATIVE_INFINITY;
        return (a - 1) * Math.log(x) + (b - 1) * Math.log(1 - x)
                - (Gamma.logGamma(a) + Gamma.logGamma(b) - Gamma.logGamma(a + b));
    }

    @Override
    public void init(final PrintStream out) {
        for (CalibrationNode n : calibrationNodes) {
            String id = n.taxa.getID();
            out.print("logP(mrca(" + id + "))\t");
            out.print("mrca.age(" + id + ")\t");
        }
    }

    @Override
    public void log(final long sample, final PrintStream out) {
        TreeInterface tree = treeInput.get();
        for (CalibrationNode n : calibrationNodes) {
            // Recompute from the current tree rather than reading a cached side-effect map,
            // so the logged value always matches the accepted state (see cladeLogDensity).
            out.print(cladeLogDensity(n, tree) + "\t");
            out.print(n.getCommonAncestor(tree).getHeight() + "\t");
        }
    }

    @Override
    public void close(final PrintStream out) {
        // nothing to do
    }

    public List<CalibrationCladePrior> getCalibrationCladePriors() {
        return cladesInput.get();
    }

    @Override
    public List<String> getArguments() {
        return List.of();
    }

    @Override
    public List<String> getConditions() {
        return List.of();
    }

    @Override
    public void sample(State state, Random random) {
    }

    @Override
    protected boolean requiresRecalculation() {

        return true;
    }

    // ------------------------------------------------------------------
    // 2. Decompose forest into subtrees of clades that have intervals overlapping with their parents
    public static class CalibrationComponent {
        private final List<CalibrationNode> members = new ArrayList<>();
        private CalibrationNode root; // exactly one root per component

        public CalibrationComponent() {
        }

        public static List<CalibrationComponent> partition(CalibrationForest calibrationForest) {
            List<CalibrationNode> all = calibrationForest.getAllNodes();
            List<CalibrationComponent> components = new ArrayList<>();
            Set<CalibrationNode> visited = new HashSet<>();

            for (CalibrationNode start : all) {
                if (visited.contains(start)) continue;

                CalibrationComponent comp = new CalibrationComponent();
                Deque<CalibrationNode> stack = new ArrayDeque<>();
                stack.push(start);
                visited.add(start);

                // --- DFS over overlapping edges ---
                while (!stack.isEmpty()) {
                    CalibrationNode cur = stack.pop();
                    comp.add(cur);

                    for (CalibrationNode other : all) {
                        if (cur == other || visited.contains(other)) continue;
                        if (overlaps(cur, other)) {
                            stack.push(other);
                            visited.add(other);
                        }
                    }
                }

                // --- Identify the highest node in the component ---
                comp.root = comp.findComponentRoot();
                components.add(comp);
            }

            return components;
        }

        private static boolean overlaps(CalibrationNode a, CalibrationNode b) {
            if (a.parent == b)
                return b.getCalibrationCladePrior().getLower() < a.getCalibrationCladePrior().getUpper(); // child inside parent's range
            if (b.parent == a)
                return a.getCalibrationCladePrior().getLower() < b.getCalibrationCladePrior().getUpper();
            return false;
        }

        public void add(CalibrationNode n) {
            members.add(n);
        }

        public CalibrationNode getRoot() {
            return root;
        }

        public List<CalibrationNode> getMembers() {
            return members;
        }

        private CalibrationNode findComponentRoot() {
            CalibrationNode rootCandidate = null;
            for (CalibrationNode n : members) {
                // A root is one whose parent is either null or not in the component,
                // or one whose parent does not overlap with it.
                if (n.parent == null || !members.contains(n.parent) || !overlaps(n, n.parent)) {
                    if (rootCandidate != null)
                        throw new IllegalStateException("Multiple roots found in component!");
                    rootCandidate = n;
                }
            }
            if (rootCandidate == null)
                throw new IllegalStateException("No root found for component!");
            return rootCandidate;
        }
    }
}
