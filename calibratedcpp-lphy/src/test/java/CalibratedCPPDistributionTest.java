import calibratedcpp.lphy.prior.Calibration;
import calibratedcpp.lphy.prior.CalibrationArray;
import calibratedcpp.lphy.tree.CPPTree;
import calibratedcpp.lphy.tree.CalibratedAgeDependentCPPTree;
import calibratedcpp.lphy.tree.CalibratedCPPTree;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.Value;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Supplier;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Statistical equivalence test between {@link CalibratedCPPTree} (constant death rate μ) and
 * {@link CalibratedAgeDependentCPPTree} (constant lifetime = 1/μ).
 *
 * <p>When the age-dependent model uses a constant lifetime = 1/μ, its effective death rate
 * equals μ exactly because both generators call the same {@code CPPUtils} functions with
 * the same numeric argument. The two should therefore produce <em>identical</em> tree
 * distributions, and all two-sample KS tests should comfortably fail to reject H₀.
 *
 * <p>Statistics compared (3 000 trees per generator):
 * <ul>
 *   <li>Total tree length (sum of all branch lengths)</li>
 *   <li>Root age (tree height) — varies only in the stem-conditioned scenario</li>
 *   <li>Mean non-root internal node age</li>
 *   <li>Max non-root internal node age (= second-oldest internal node)</li>
 *   <li>Pybus &amp; Harvey (2000) gamma statistic — measures rate variation through time</li>
 * </ul>
 */
public class CalibratedCPPDistributionTest {

    private static final int N = 3000;

    /**
     * Suite-wide Bonferroni correction. There are ~{@value #NUM_KS_ASSERTIONS} KS assertions across
     * all scenarios in this class; at a per-assertion α of 0.01 the family-wise false-positive rate
     * is ~0.4 (≈ 49 × 0.01), which reds the suite roughly one run in three. Controlling the
     * family-wise error at {@value #FAMILY_ALPHA} sets the per-assertion level to
     * FAMILY_ALPHA / NUM_KS_ASSERTIONS. Every null here is TRUE (identical models), so the reduced
     * per-test power costs nothing; a real regression still produces D far above the critical value.
     * (If you add or remove assertKS calls, update NUM_KS_ASSERTIONS — being off by a few is harmless.)
     */
    private static final int NUM_KS_ASSERTIONS = 49;
    private static final double FAMILY_ALPHA = 1e-3;
    private static final double KS_ALPHA = FAMILY_ALPHA / NUM_KS_ASSERTIONS;

    // =========================================================================
    // Tree summary statistics
    // =========================================================================

    static double treeLength(TimeTree tree) {
        double total = 0;
        double branchLength = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            TimeTreeNode node = tree.getNodeByIndex(i);
            if (!node.isRoot()) total += node.getParent().getAge() - node.getAge();
            if (!node.isRoot()) branchLength += node.getBranchDuration();
        }
        if (total != branchLength){
            throw new RuntimeException("Tree length does not match, tree length: " + tree + ", branch length: " + branchLength);
        }
        return total;
    }

    /** Internal node ages sorted ascending (youngest → root). */
    static double[] sortedInternalAges(TimeTree tree) {
        double[] ages = new double[tree.getNodeCount()];
        int idx = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            TimeTreeNode nd = tree.getNodeByIndex(i);
            if (!nd.isLeaf()) ages[idx++] = nd.getAge();
        }
        ages = Arrays.copyOf(ages, idx);
        Arrays.sort(ages);
        return ages;
    }

    /**
     * Pybus &amp; Harvey (2000) gamma statistic.
     *
     * <p>Measures whether lineage diversification rates speed up (+γ) or slow down (−γ)
     * toward the present compared with a constant-rate expectation.
     * Under a constant-rate BD process, γ → N(0, 1) as n → ∞.
     *
     * <p>Formula: sort internal node ages ascending (a₁ ≤ … ≤ a_{n−1}).
     * At the k-th interval from the present (a_{k−1} to aₖ, a₀ = 0) there are n−k+1 lineages.
     * Lineage-weighted interval wₖ = (n−k+1)·(aₖ − a_{k−1}).
     * Total W = Σwₖ.  Cumulative Tₖ = Σⱼ₌₁ᵏ wⱼ₋₁.
     * γ = [(1/(n−2))·Σₖ₌₁ⁿ⁻² Tₖ − W/2] / [W/√(12(n−2))].
     */
    static double gammaStatistic(TimeTree tree) {
        double[] a = sortedInternalAges(tree);
        int ni = a.length;          // n-1 internal nodes
        int n  = ni + 1;            // number of taxa

        double[] g = new double[ni];
        g[0] = a[0];
        for (int i = 1; i < ni; i++) g[i] = a[i] - a[i - 1];

        double[] w = new double[ni];
        double W = 0;
        for (int i = 0; i < ni; i++) { w[i] = (n - i) * g[i]; W += w[i]; }

        double sumT = 0, Tk = 0;
        for (int k = 1; k <= ni - 1; k++) { Tk += w[k - 1]; sumT += Tk; }

        double num = (1.0 / (n - 2)) * sumT - W / 2.0;
        double den = W / Math.sqrt(12.0 * (n - 2));
        return num / den;
    }

    /** Recursive leaf count for a subtree. */
    static int leafCount(TimeTreeNode node) {
        if (node.isLeaf()) return 1;
        int c = 0;
        for (TimeTreeNode ch : node.getChildren()) c += leafCount(ch);
        return c;
    }

    /**
     * Colless imbalance index (Colless 1982).
     * For each binary internal node: |left_leaves − right_leaves|.  Sum over all internal nodes.
     * Zero for a perfectly balanced tree; maximum for a caterpillar.
     */
    static double collessIndex(TimeTree tree) {
        double sum = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            TimeTreeNode nd = tree.getNodeByIndex(i);
            if (!nd.isLeaf() && nd.getChildren().size() == 2) {
                sum += Math.abs(leafCount(nd.getChildren().get(0))
                              - leafCount(nd.getChildren().get(1)));
            }
        }
        return sum;
    }

    /**
     * Number of cherries: internal nodes with exactly two leaf children.
     * Under the Yule model E[cherries] ≈ n/3 for large n.
     */
    static double numCherries(TimeTree tree) {
        int count = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            TimeTreeNode nd = tree.getNodeByIndex(i);
            if (!nd.isLeaf()) {
                long leafKids = nd.getChildren().stream()
                        .filter(TimeTreeNode::isLeaf).count();
                if (leafKids == 2) count++;
            }
        }
        return count;
    }

    // --- Aldous (1996/2001) beta statistic helpers ---

    /** log Γ(x) via Stirling's series, accurate to < 1e-10 for x > 1. */
    static double logGamma(double x) {
        double r = 0;
        while (x < 8) { r -= Math.log(x); x += 1; }
        double xi = 1.0 / x, xsq = xi * xi;
        return r + (x - 0.5) * Math.log(x) - x + 0.5 * Math.log(2 * Math.PI)
                + xi * (1.0 / 12 - xsq * (1.0 / 360 - xsq / 1260));
    }

    static double digamma(double x) {
        double r = 0;
        while (x < 8) { r -= 1.0 / x; x += 1; }
        double xi = 1.0 / x, xsq = xi * xi;
        return r + Math.log(x) - 0.5 * xi - xsq * (1.0 / 12 - xsq * (1.0 / 120 - xsq / 252));
    }

    static double logSumExp(double[] v) {
        double max = Double.NEGATIVE_INFINITY;
        for (double x : v) if (x > max) max = x;
        double s = 0;
        for (double x : v) s += Math.exp(x - max);
        return max + Math.log(s);
    }

    /**
     * Derivative of the full Aldous beta log-likelihood (with partition-function correction).
     *
     * <p>The beta-splitting model assigns to a split (a, b) of n = a+b nodes:
     *   p_β(a, b | n) ∝ B(a+β+1, b+β+1)
     * The MLE score is:
     *   ℓ'(β) = Σ_v { [ψ(a_v+β+1) + ψ(b_v+β+1)]
     *                 − E_{p_β}[ψ(J+β+1) + ψ(n_v−J+β+1)] }
     * where the expectation is over J ~ p_β(·, n_v−· | n_v).
     * This correction term is essential; without it the score is always negative.
     */
    static double betaScore(double beta, List<int[]> splits) {
        double score = 0;
        for (int[] ab : splits) {
            int a = ab[0], b = ab[1], nv = a + b;
            score += digamma(a + beta + 1) + digamma(b + beta + 1);
            // Compute partition weights log B(j+β+1, nv-j+β+1) for j=1..nv-1
            double[] logW = new double[nv - 1];
            for (int j = 1; j < nv; j++)
                logW[j - 1] = logGamma(j + beta + 1) + logGamma(nv - j + beta + 1)
                             - logGamma(nv + 2 * beta + 2);
            double logZ = logSumExp(logW);
            double expected = 0;
            for (int j = 1; j < nv; j++)
                expected += Math.exp(logW[j - 1] - logZ)
                          * (digamma(j + beta + 1) + digamma(nv - j + beta + 1));
            score -= expected;
        }
        return score;
    }

    /**
     * Aldous (1996/2001) beta statistic: MLE of β in the beta-splitting model via bisection.
     *
     * <p>β ≈ 0  corresponds to the Yule (pure-birth) model;
     * β &lt; 0  indicates more imbalanced trees than Yule;
     * β &gt; 0  indicates more balanced trees.
     * The search range [−1.9, 10] covers almost all trees; at the boundaries the
     * likelihood is flat, so we return the boundary value rather than diverging.
     */
    static double betaStatistic(TimeTree tree) {
        List<int[]> splits = new ArrayList<>();
        for (int i = 0; i < tree.getNodeCount(); i++) {
            TimeTreeNode nd = tree.getNodeByIndex(i);
            if (!nd.isLeaf() && nd.getChildren().size() == 2) {
                int a = leafCount(nd.getChildren().get(0));
                int b = leafCount(nd.getChildren().get(1));
                splits.add(new int[]{a, b});
            }
        }
        if (splits.isEmpty()) return 0.0;

        final double LO = -1.9, HI = 10.0;
        double fLo = betaScore(LO, splits);
        double fHi = betaScore(HI, splits);

        // No sign change → MLE is near the boundary with the larger absolute score
        if (fLo * fHi > 0)
            return Math.abs(fLo) < Math.abs(fHi) ? LO : HI;

        // Bisection: keep [lo, hi] such that score(lo) > 0 and score(hi) < 0
        double lo = fLo > 0 ? LO : HI;
        double hi = fLo > 0 ? HI : LO;
        for (int k = 0; k < 60; k++) {
            double mid = (lo + hi) / 2;
            if (betaScore(mid, splits) > 0) lo = mid; else hi = mid;
        }
        return (lo + hi) / 2;
    }

    // =========================================================================
    // Two-sample Kolmogorov-Smirnov test
    // =========================================================================

    static double ksStatistic(double[] a, double[] b) {
        double[] sa = a.clone(); Arrays.sort(sa);
        double[] sb = b.clone(); Arrays.sort(sb);
        int n = sa.length, m = sb.length;
        int i = 0, j = 0;
        double maxD = 0;
        while (i < n || j < m) {
            double va = i < n ? sa[i] : Double.MAX_VALUE;
            double vb = j < m ? sb[j] : Double.MAX_VALUE;
            if (va <= vb) i++;
            if (vb <= va) j++;
            maxD = Math.max(maxD, Math.abs((double) i / n - (double) j / m));
        }
        return maxD;
    }

    /**
     * Kolmogorov critical constant c such that
     * P(D_{n,m} > c / √N_eff) ≈ α, where N_eff = n·m / (n+m).
     *
     * <p>Closed form of the asymptotic two-sample KS distribution: c(α) = √(-½ ln(α/2)).
     * This reproduces the standard table (α=0.10→1.2238, 0.05→1.3581, 0.01→1.6276,
     * 0.001→1.9499) and, unlike a fixed table, extends to the small per-test α produced by
     * the suite-wide Bonferroni correction.
     */
    static double criticalKS(double alpha) {
        return Math.sqrt(-0.5 * Math.log(alpha / 2.0));
    }

    /**
     * Asserts that the two-sample KS statistic D is below the critical value
     * at significance level {@link #KS_ALPHA}.
     *
     * <p>Using the critical-value form rather than a p-value series avoids
     * numerical instability when D ≈ 0 (identical distributions), where the
     * alternating Kolmogorov series does not converge.
     */
    static void assertKS(double[] cpp, double[] ad, String label) {
        int n = cpp.length, m = ad.length;
        double Neff = (double) n * m / (n + m);
        double d = ksStatistic(cpp, ad);
        double dCrit = criticalKS(KS_ALPHA) / Math.sqrt(Neff);
        assertTrue(d < dCrit,
                String.format("KS test FAILED for '%s': D=%.4f, critical=%.4f (α=%.2f)",
                        label, d, dCrit, KS_ALPHA));
    }

    // =========================================================================
    // Batch sampling
    // =========================================================================

    record Stats(double[] length, double[] rootAge, double[] meanNonRootAge,
                 double[] maxNonRootAge, double[] gamma,
                 double[] colless, double[] cherries, double[] beta) {}

    static Stats collect(int n, Supplier<TimeTree> sampler) {
        double[] len  = new double[n], root    = new double[n],
                 mean = new double[n], maxNR   = new double[n],
                 gam  = new double[n], col     = new double[n],
                 cher = new double[n], bet     = new double[n];
        for (int i = 0; i < n; i++) {
            TimeTree tree = sampler.get();
            len[i]  = treeLength(tree);
            gam[i]  = gammaStatistic(tree);
            col[i]  = collessIndex(tree);
            cher[i] = numCherries(tree);
            bet[i]  = betaStatistic(tree);
            double[] ages = sortedInternalAges(tree);  // ascending; last = root
            root[i] = ages[ages.length - 1];
            int nr = ages.length - 1;
            mean[i] = nr > 0 ? Arrays.stream(ages, 0, nr).average().orElse(0) : 0;
            maxNR[i] = nr > 0 ? ages[nr - 1] : 0;
        }
        return new Stats(len, root, mean, maxNR, gam, col, cher, bet);
    }

    // =========================================================================
    // Tests
    // =========================================================================

    /**
     * Root-conditioned scenario: a single calibration covers all taxa so the root
     * age is fixed.  The free statistics are total tree length, non-root internal
     * node ages, and the gamma statistic.
     *
     * <p>Equivalence condition: {@code lifetime = 1/μ = 2.0} (constant),
     * so {@code effectiveDeathRate = 1/lifetime = μ = 0.5} exactly.
     */
    @Test
    void rootConditionedDistributionsMatchCPP() {
        final double lambda = 2.0, mu = 0.5, rho = 0.3, lifetime = 1.0 / mu;
        final int n = 5;
        final double rootAge = 8.0;
        final String[] all = {"a", "b", "c", "d", "e"};

        CalibratedCPPTree cpp = new CalibratedCPPTree(
                new Value<>("", lambda), new Value<>("", mu), null, null,
                new Value<>("", rho), new Value<>("", n),
                new Value<>("", new CalibrationArray(new Calibration[]{
                        new Calibration(all, rootAge)})),
                null, null, null);

        CalibratedAgeDependentCPPTree ad = new CalibratedAgeDependentCPPTree(
                new Value<>("", lambda), new Value<>("", rho), new Value<>("", n),
                new Value<>("", lifetime),
                new Value<>("", new CalibrationArray(new Calibration[]{
                        new Calibration(all, rootAge)})),
                null, null, null);

        Stats cs = collect(N, () -> cpp.sample().value());
        Stats as = collect(N, () -> ad.sample().value());

        assertKS(cs.length(),         as.length(),         "root-cond: total tree length");
        assertKS(cs.meanNonRootAge(), as.meanNonRootAge(), "root-cond: mean non-root internal age");
        assertKS(cs.maxNonRootAge(),  as.maxNonRootAge(),  "root-cond: max non-root internal age");
        assertKS(cs.gamma(),          as.gamma(),          "root-cond: gamma statistic");
        assertKS(cs.colless(),        as.colless(),        "root-cond: Colless imbalance");
        assertKS(cs.cherries(),       as.cherries(),       "root-cond: number of cherries");
        assertKS(cs.beta(),           as.beta(),           "root-cond: Aldous beta statistic");
    }

    /**
     * Stem-conditioned scenario: a calibration covers a subclade while the stem age
     * is fixed externally.  The root age itself now varies, giving an additional
     * dimension to compare.
     *
     * <p>Same equivalence condition: {@code lifetime = 1/μ}.
     */
    @Test
    void stemConditionedDistributionsMatchCPP() {
        final double lambda = 2.0, mu = 0.5, rho = 0.3, lifetime = 1.0 / mu;
        final int n = 5;
        final double cladeAge = 3.0, stemAge = 10.0;
        final String[] clade = {"a", "b", "c"};

        CalibratedCPPTree cpp = new CalibratedCPPTree(
                new Value<>("", lambda), new Value<>("", mu), null, null,
                new Value<>("", rho), new Value<>("", n),
                new Value<>("", new CalibrationArray(new Calibration[]{
                        new Calibration(clade, cladeAge)})),
                null, new Value<>("", stemAge), null);

        CalibratedAgeDependentCPPTree ad = new CalibratedAgeDependentCPPTree(
                new Value<>("", lambda), new Value<>("", rho), new Value<>("", n),
                new Value<>("", lifetime),
                new Value<>("", new CalibrationArray(new Calibration[]{
                        new Calibration(clade, cladeAge)})),
                null, new Value<>("", stemAge), null);

        Stats cs = collect(N, () -> cpp.sample().value());
        Stats as = collect(N, () -> ad.sample().value());

        assertKS(cs.rootAge(),        as.rootAge(),        "stem-cond: root age (tree height)");
        assertKS(cs.length(),         as.length(),         "stem-cond: total tree length");
        assertKS(cs.meanNonRootAge(), as.meanNonRootAge(), "stem-cond: mean non-root internal age");
        assertKS(cs.maxNonRootAge(),  as.maxNonRootAge(),  "stem-cond: max non-root internal age");
        assertKS(cs.gamma(),          as.gamma(),          "stem-cond: gamma statistic");
        assertKS(cs.colless(),        as.colless(),        "stem-cond: Colless imbalance");
        assertKS(cs.cherries(),       as.cherries(),       "stem-cond: number of cherries");
        assertKS(cs.beta(),           as.beta(),           "stem-cond: Aldous beta statistic");
    }

    /**
     * Nested calibrations: outer clade of 4 taxa at age 4.0, inner clade of 2 at age 2.0,
     * plus 2 uncalibrated taxa (n=6 total).  With 5 internal nodes but only 2 fixed by
     * calibrations, the remaining 3 vary — making the statistics genuinely stochastic.
     * Tests the recursive sub-generator path in {@link CalibratedAgeDependentCPPTree}.
     */
    @Test
    void nestedCalibrationsDistributionsMatchCPP() {
        final double lambda = 2.0, mu = 0.5, rho = 0.4, lifetime = 1.0 / mu;
        final int n = 6;
        final String[] outer = {"1", "2", "3", "4"};
        final String[] inner = {"1", "2"};

        CalibratedCPPTree cpp = new CalibratedCPPTree(
                new Value<>("", lambda), new Value<>("", mu), null, null,
                new Value<>("", rho), new Value<>("", n),
                new Value<>("", new CalibrationArray(new Calibration[]{
                        new Calibration(outer, 4.0),
                        new Calibration(inner, 2.0)})),
                null, new Value<>("", 10.0), null);

        CalibratedAgeDependentCPPTree ad = new CalibratedAgeDependentCPPTree(
                new Value<>("", lambda), new Value<>("", rho), new Value<>("", n),
                new Value<>("", lifetime),
                new Value<>("", new CalibrationArray(new Calibration[]{
                        new Calibration(outer, 4.0),
                        new Calibration(inner, 2.0)})),
                null, new Value<>("", 10.0), null);

        Stats cs = collect(N, () -> cpp.sample().value());
        Stats as = collect(N, () -> ad.sample().value());

        assertKS(cs.rootAge(),        as.rootAge(),        "nested: root age");
        assertKS(cs.length(),         as.length(),         "nested: total tree length");
        assertKS(cs.meanNonRootAge(), as.meanNonRootAge(), "nested: mean non-root internal age");
        assertKS(cs.maxNonRootAge(),  as.maxNonRootAge(),  "nested: max non-root internal age");
        assertKS(cs.gamma(),          as.gamma(),          "nested: gamma statistic");
        assertKS(cs.colless(),        as.colless(),        "nested: Colless imbalance");
        assertKS(cs.cherries(),       as.cherries(),       "nested: number of cherries");
        assertKS(cs.beta(),           as.beta(),           "nested: Aldous beta statistic");
    }

    /**
     * Replicates the fixed_stem_validation.lphy scenario (100 taxa, 4 calibrations, stemAge=3).
     * Mirrors the R validation script lphy_added_fixed_stem_validation_100.R which compares
     * CPP and lphyCPP trees on: tree height, cherries, Colless, gamma, tree length, beta.
     * λ=2, μ=1, ρ=0.1, lifetime=1/μ=1.0  →  effectiveDeathRate = μ exactly.
     */
    @Test
    void fixedStemValidation100TipDistributionsMatchCPP() {
        final double lambda = 2.0, mu = 1.0, rho = 0.1, lifetime = 1.0 / mu;
        final int n = 100;
        final double stemAge = 3.0;

        // clade1: leaf_1, leaf_2  @ age 1.0
        String[] clade1 = {"leaf_1", "leaf_2"};
        // clade2: leaf_45..leaf_51  @ age 1.2
        String[] clade2 = new String[7];
        for (int i = 0; i < 7; i++) clade2[i] = "leaf_" + (45 + i);
        // clade3: leaf_45..leaf_55  @ age 1.5  (superset of clade2)
        String[] clade3 = new String[11];
        for (int i = 0; i < 11; i++) clade3[i] = "leaf_" + (45 + i);
        // clade4: leaf_87..leaf_90  @ age 2.0
        String[] clade4 = {"leaf_87", "leaf_88", "leaf_89", "leaf_90"};

        Calibration[] cals = {
                new Calibration(clade1, 1.0),
                new Calibration(clade2, 1.2),
                new Calibration(clade3, 1.5),
                new Calibration(clade4, 2.0)
        };

        CalibratedCPPTree cpp = new CalibratedCPPTree(
                new Value<>("", lambda), new Value<>("", mu), null, null,
                new Value<>("", rho), new Value<>("", n),
                new Value<>("", new CalibrationArray(cals)),
                null, new Value<>("", stemAge), null);

        CalibratedAgeDependentCPPTree ad = new CalibratedAgeDependentCPPTree(
                new Value<>("", lambda), new Value<>("", rho), new Value<>("", n),
                new Value<>("", lifetime),
                new Value<>("", new CalibrationArray(cals)),
                null, new Value<>("", stemAge), null);

        final int N100 = 1000;
        Stats cs = collect(N100, () -> cpp.sample().value());
        Stats as = collect(N100, () -> ad.sample().value());

        assertKS(cs.rootAge(),        as.rootAge(),        "100tip: root age (tree height)");
        assertKS(cs.length(),         as.length(),         "100tip: total tree length");
        assertKS(cs.meanNonRootAge(), as.meanNonRootAge(), "100tip: mean non-root internal age");
        assertKS(cs.gamma(),          as.gamma(),          "100tip: gamma statistic");
        assertKS(cs.colless(),        as.colless(),        "100tip: Colless imbalance");
        assertKS(cs.cherries(),       as.cherries(),       "100tip: number of cherries");
        assertKS(cs.beta(),           as.beta(),           "100tip: Aldous beta statistic");
    }

    /**
     * No-calibrations path: when calibrations is null and rootAge is provided, both
     * CalibratedCPPTree and CalibratedAgeDependentCPPTree delegate to CPPTree.
     * All three should therefore produce identical tree distributions.
     */
    @Test
    void noCalibrationRootAgeMatchesCPPTree() {
        final double lambda = 2.0, mu = 0.5, rho = 0.3, lifetime = 1.0 / mu;
        final int n = 8;
        final double rootAge = 6.0;

        CPPTree reference = new CPPTree(
                new Value<>("", lambda), new Value<>("", mu), null, null,
                new Value<>("", rho), null, new Value<>("", n),
                new Value<>("", rootAge), null);

        CalibratedCPPTree cpp = new CalibratedCPPTree(
                new Value<>("", lambda), new Value<>("", mu), null, null,
                new Value<>("", rho), new Value<>("", n),
                null, null, null, new Value<>("", rootAge));

        CalibratedAgeDependentCPPTree ad = new CalibratedAgeDependentCPPTree(
                new Value<>("", lambda), new Value<>("", rho), new Value<>("", n),
                new Value<>("", lifetime),
                null, null, null, new Value<>("", rootAge));

        Stats refStats = collect(N, () -> reference.sample().value());
        Stats cppStats = collect(N, () -> cpp.sample().value());
        Stats adStats  = collect(N, () -> ad.sample().value());

        // CalibratedCPPTree (rootAge only) must match CPPTree
        assertKS(refStats.length(),   cppStats.length(),   "no-cal: CPP vs ref tree length");
        assertKS(refStats.gamma(),    cppStats.gamma(),    "no-cal: CPP vs ref gamma");
        assertKS(refStats.colless(),  cppStats.colless(),  "no-cal: CPP vs ref Colless");
        assertKS(refStats.cherries(), cppStats.cherries(), "no-cal: CPP vs ref cherries");
        assertKS(refStats.beta(),     cppStats.beta(),     "no-cal: CPP vs ref beta");

        // CalibratedAgeDependentCPPTree (rootAge only) must also match CPPTree
        assertKS(refStats.length(),   adStats.length(),    "no-cal: AD vs ref tree length");
        assertKS(refStats.gamma(),    adStats.gamma(),     "no-cal: AD vs ref gamma");
        assertKS(refStats.colless(),  adStats.colless(),   "no-cal: AD vs ref Colless");
        assertKS(refStats.cherries(), adStats.cherries(),  "no-cal: AD vs ref cherries");
        assertKS(refStats.beta(),     adStats.beta(),      "no-cal: AD vs ref beta");
    }

    /**
     * Disjoint calibrations (two independent clades).
     * Tests the multi-maximal-calibration code path.
     */
    @Test
    void disjointCalibrationsDistributionsMatchCPP() {
        final double lambda = 2.0, mu = 0.5, rho = 0.4, lifetime = 1.0 / mu;
        final int n = 7;
        final String[] c1 = {"a", "b", "c"};
        final String[] c2 = {"d", "e"};

        CalibratedCPPTree cpp = new CalibratedCPPTree(
                new Value<>("", lambda), new Value<>("", mu), null, null,
                new Value<>("", rho), new Value<>("", n),
                new Value<>("", new CalibrationArray(new Calibration[]{
                        new Calibration(c1, 3.0),
                        new Calibration(c2, 2.0)})),
                null, new Value<>("", 8.0), null);

        CalibratedAgeDependentCPPTree ad = new CalibratedAgeDependentCPPTree(
                new Value<>("", lambda), new Value<>("", rho), new Value<>("", n),
                new Value<>("", lifetime),
                new Value<>("", new CalibrationArray(new Calibration[]{
                        new Calibration(c1, 3.0),
                        new Calibration(c2, 2.0)})),
                null, new Value<>("", 8.0), null);

        Stats cs = collect(N, () -> cpp.sample().value());
        Stats as = collect(N, () -> ad.sample().value());

        assertKS(cs.rootAge(),        as.rootAge(),        "disjoint: root age");
        assertKS(cs.length(),         as.length(),         "disjoint: total tree length");
        assertKS(cs.meanNonRootAge(), as.meanNonRootAge(), "disjoint: mean non-root internal age");
        assertKS(cs.maxNonRootAge(),  as.maxNonRootAge(),  "disjoint: max non-root internal age");
        assertKS(cs.gamma(),          as.gamma(),          "disjoint: gamma statistic");
        assertKS(cs.colless(),        as.colless(),        "disjoint: Colless imbalance");
        assertKS(cs.cherries(),       as.cherries(),       "disjoint: number of cherries");
        assertKS(cs.beta(),           as.beta(),           "disjoint: Aldous beta statistic");
    }
}
