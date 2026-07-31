import calibratedcpp.lphy.prior.Calibration;
import calibratedcpp.lphy.prior.CalibrationArray;
import calibratedcpp.lphy.tree.CalibratedAgeDependentCPPTree;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.Value;
import org.junit.jupiter.api.RepeatedTest;
import org.junit.jupiter.api.Test;

import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Tests for {@link CalibratedAgeDependentCPPTree} in the calibratedcpp-lphy module.
 *
 * <p>Covers:
 * <ul>
 *   <li>Constructor validation (null checks)</li>
 *   <li>Root-conditioned sampling (calibration covers all taxa)</li>
 *   <li>Stem-conditioned sampling (explicit stem age)</li>
 *   <li>Nested calibrations (inner clade age honoured)</li>
 *   <li>Disjoint calibrations (two independent clades)</li>
 *   <li>Outgroup taxa via otherNames</li>
 *   <li>Stem-age warning when root calibration is given alongside a stem age</li>
 *   <li>getParams / setParam round-trip</li>
 *   <li>getOrigin() and getRootCondition() after sampling</li>
 * </ul>
 */
public class CalibratedAgeDependentCPPTreeTest {

    // -------------------------------------------------------------------------
    // Helpers
    // -------------------------------------------------------------------------

    private static CalibratedAgeDependentCPPTree generate(
            double lambda, double rho, int n, double meanLifetime,
            double stemAge, Calibration... cals) {
        return new CalibratedAgeDependentCPPTree(
                new Value<>("", lambda),
                new Value<>("", rho),
                new Value<>("", n),
                new Value<>("", meanLifetime),
                new Value<>("", new CalibrationArray(cals)),
                null,
                new Value<>("", stemAge),
                null);
    }

    private static CalibratedAgeDependentCPPTree generateNoStem(
            double lambda, double rho, int n, double meanLifetime,
            Calibration... cals) {
        return new CalibratedAgeDependentCPPTree(
                new Value<>("", lambda),
                new Value<>("", rho),
                new Value<>("", n),
                new Value<>("", meanLifetime),
                new Value<>("", new CalibrationArray(cals)),
                null, null, null);
    }

    /** Every internal node must be strictly older than each of its children. */
    private static void assertAgesOrder(TimeTree tree) {
        for (int i = 0; i < tree.getNodeCount(); i++) {
            TimeTreeNode node = tree.getNodeByIndex(i);
            for (TimeTreeNode child : node.getChildren()) {
                assertTrue(node.getAge() > child.getAge(),
                        "parent age " + node.getAge() + " must be > child age " + child.getAge());
            }
        }
    }

    /** Returns the MRCA node for the given taxa, or null if not found. */
    private static TimeTreeNode findMRCA(TimeTree tree, String... taxa) {
        List<String> target = List.of(taxa);
        for (int i = 0; i < tree.getNodeCount(); i++) {
            TimeTreeNode node = tree.getNodeByIndex(i);
            if (!node.isLeaf()) {
                List<String> leaves = List.of(node.getLeafNames());
                if (leaves.size() == target.size() && leaves.containsAll(target)) {
                    return node;
                }
            }
        }
        return null;
    }

    // -------------------------------------------------------------------------
    // Constructor validation
    // -------------------------------------------------------------------------

    @Test
    void nullCalibrationsAndRootAgeSamplesRandomStem() {
        // Neither calibrations nor rootAge is now valid: the uncalibrated path samples the stem/origin
        // from the process (n given -> random stem). This replaces the old "must provide calibrations
        // or rootAge" invariant, removed to enable the random-stem / random-N modes.
        CalibratedAgeDependentCPPTree gen = new CalibratedAgeDependentCPPTree(
                new Value<>("", 1.0), new Value<>("", 0.5), new Value<>("", 3),
                new Value<>("", 2.0),
                null, null, null, null);
        TimeTree tree = gen.sample().value();

        int leaves = 0;
        for (int i = 0; i < tree.getNodeCount(); i++) {
            if (tree.getNodeByIndex(i).isLeaf()) leaves++;
        }
        assertEquals(3, leaves, "n=3 tips sampled");
        assertAgesOrder(tree);
    }

    @Test
    void nullLifetimeThrows() {
        assertThrows(NullPointerException.class, () ->
                new CalibratedAgeDependentCPPTree(
                        new Value<>("", 1.0), new Value<>("", 0.5), new Value<>("", 3),
                        null,
                        new Value<>("", new CalibrationArray(new Calibration[]{
                                new Calibration(new String[]{"a", "b", "c"}, 2.0)})),
                        null, null, null));
    }

    // -------------------------------------------------------------------------
    // Root-conditioned: single calibration covering all taxa
    // -------------------------------------------------------------------------

    @RepeatedTest(10)
    void rootConditionedRootAgeIsExact() {
        String[] all = {"a", "b", "c"};
        double rootAge = 5.0;
        CalibratedAgeDependentCPPTree g = generateNoStem(1.5, 0.5, 3, 1.0,
                new Calibration(all, rootAge));

        TimeTree tree = g.sample().value();

        assertTrue(g.getRootCondition(), "should be root-conditioned");
        assertEquals(rootAge, tree.getRoot().getAge(), 1e-8, "root age");
        assertEquals(3, tree.getLeafNodes().size(), "leaf count");
        assertAgesOrder(tree);
    }

    @RepeatedTest(10)
    void rootConditionedLeafNamesComplete() {
        String[] all = {"x", "y", "z", "w"};
        CalibratedAgeDependentCPPTree g = generateNoStem(1.0, 0.4, 4, 2.0,
                new Calibration(all, 6.0));
        TimeTree tree = g.sample().value();

        List<String> names = Arrays.asList(tree.getTaxaNames());
        for (String taxon : all) {
            assertTrue(names.contains(taxon), "taxon " + taxon + " missing from tree");
        }
    }

    // -------------------------------------------------------------------------
    // Stem-conditioned: explicit stem age, single calibration on a subclade
    // -------------------------------------------------------------------------

    @RepeatedTest(10)
    void stemConditionedCladeAgeIsExact() {
        String[] clade = {"a", "b", "c"};
        double cladeAge = 3.0, stemAge = 8.0;
        CalibratedAgeDependentCPPTree g = generate(1.5, 0.5, 5, 1.0, stemAge,
                new Calibration(clade, cladeAge));

        TimeTree tree = g.sample().value();

        assertFalse(g.getRootCondition(), "should NOT be root-conditioned");
        assertEquals(5, tree.getLeafNodes().size(), "leaf count");
        assertAgesOrder(tree);

        TimeTreeNode mrca = findMRCA(tree, clade);
        assertNotNull(mrca, "MRCA of calibrated clade not found");
        assertEquals(cladeAge, mrca.getAge(), 1e-8, "calibrated clade age");
        assertTrue(tree.getRoot().getAge() <= stemAge,
                "root age " + tree.getRoot().getAge() + " must be <= stem age " + stemAge);
    }

    @RepeatedTest(10)
    void stemConditionedOriginEqualsConditionAge() {
        double stemAge = 10.0;
        CalibratedAgeDependentCPPTree g = generate(1.5, 0.5, 5, 2.0, stemAge,
                new Calibration(new String[]{"a", "b"}, 3.0));

        g.sample();

        assertEquals(stemAge, g.getOrigin().value(), 1e-8,
                "getOrigin() must equal the provided stem age");
    }

    // -------------------------------------------------------------------------
    // Nested calibrations
    // -------------------------------------------------------------------------

    @RepeatedTest(10)
    void nestedCalibrationsAgeHonoured() {
        String[] outer = {"1", "2", "3"};
        String[] inner = {"1", "2"};
        CalibratedAgeDependentCPPTree g = generate(1.0, 0.5, 3, 1.5, 8.0,
                new Calibration(outer, 4.0),
                new Calibration(inner, 2.0));

        TimeTree tree = g.sample().value();

        assertEquals(3, tree.getLeafNodes().size(), "leaf count");
        assertAgesOrder(tree);

        TimeTreeNode outerMRCA = findMRCA(tree, outer);
        TimeTreeNode innerMRCA = findMRCA(tree, inner);
        assertNotNull(outerMRCA, "outer MRCA not found");
        assertNotNull(innerMRCA, "inner MRCA not found");
        assertEquals(4.0, outerMRCA.getAge(), 1e-8, "outer clade age");
        assertEquals(2.0, innerMRCA.getAge(), 1e-8, "inner clade age");
    }

    // -------------------------------------------------------------------------
    // Disjoint calibrations
    // -------------------------------------------------------------------------

    @RepeatedTest(20)
    void disjointCalibrationsAgesHonoured() {
        String[] clade1 = {"a", "b", "c"};
        String[] clade2 = {"d", "e"};
        CalibratedAgeDependentCPPTree g = generate(1.0, 0.5, 7, 1.0, 8.0,
                new Calibration(clade1, 3.0),
                new Calibration(clade2, 2.0));

        TimeTree tree = g.sample().value();

        assertEquals(7, tree.getLeafNodes().size(), "leaf count");
        assertAgesOrder(tree);

        TimeTreeNode mrca1 = findMRCA(tree, clade1);
        TimeTreeNode mrca2 = findMRCA(tree, clade2);
        assertNotNull(mrca1, "MRCA of clade1 not found");
        assertNotNull(mrca2, "MRCA of clade2 not found");
        assertEquals(3.0, mrca1.getAge(), 1e-8, "clade1 age");
        assertEquals(2.0, mrca2.getAge(), 1e-8, "clade2 age");
    }

    // -------------------------------------------------------------------------
    // otherNames: explicitly named outgroup taxa
    // -------------------------------------------------------------------------

    @RepeatedTest(10)
    void otherNamesAppearsInTree() {
        String[] clade = {"a", "b"};
        String[] others = {"x", "y", "z"};
        CalibratedAgeDependentCPPTree g = new CalibratedAgeDependentCPPTree(
                new Value<>("", 1.0), new Value<>("", 0.5), new Value<>("", 5),
                new Value<>("", 1.0),
                new Value<>("", new CalibrationArray(new Calibration[]{new Calibration(clade, 2.0)})),
                new Value<>("", others),
                new Value<>("", 8.0),
                null);

        TimeTree tree = g.sample().value();

        assertEquals(5, tree.getLeafNodes().size(), "leaf count");
        assertAgesOrder(tree);

        List<String> names = Arrays.asList(tree.getTaxaNames());
        for (String o : others) {
            assertTrue(names.contains(o), "otherName '" + o + "' missing from tree");
        }
    }


    // -------------------------------------------------------------------------
    // Root condition and origin accessors
    // -------------------------------------------------------------------------

    @Test
    void getRootConditionFalseForPartialCalibration() {
        // calibration covers only 2 of 4 taxa → not root-conditioned
        CalibratedAgeDependentCPPTree g = generate(1.0, 0.5, 4, 1.0, 8.0,
                new Calibration(new String[]{"a", "b"}, 3.0));
        assertFalse(g.getRootCondition());
    }

    @Test
    void getRootConditionTrueForFullCalibration() {
        String[] all = {"a", "b", "c"};
        CalibratedAgeDependentCPPTree g = generateNoStem(1.0, 0.5, 3, 1.0,
                new Calibration(all, 5.0));
        g.sample(); // triggers rootConditioned = true
        assertTrue(g.getRootCondition());
    }

    // -------------------------------------------------------------------------
    // rootAge without calibrations
    // -------------------------------------------------------------------------

    @RepeatedTest(10)
    void rootAgeOnlyNoCalibrationsSamplesCorrectly() {
        double rootAge = 5.0;
        int n = 6;
        CalibratedAgeDependentCPPTree g = new CalibratedAgeDependentCPPTree(
                new Value<>("", 1.5), new Value<>("", 0.5), new Value<>("", n),
                new Value<>("", 1.0),
                null, null, null,
                new Value<>("", rootAge));

        TimeTree tree = g.sample().value();

        assertEquals(n, tree.getLeafNodes().size(), "leaf count");
        assertEquals(rootAge, tree.getRoot().getAge(), 1e-8, "root age must match provided rootAge");
        assertAgesOrder(tree);
    }

    @RepeatedTest(10)
    void rootAgeOnlyWithOtherNamesSamplesCorrectly() {
        double rootAge = 4.0;
        String[] names = {"alpha", "beta", "gamma", "delta"};
        CalibratedAgeDependentCPPTree g = new CalibratedAgeDependentCPPTree(
                new Value<>("", 2.0), new Value<>("", 0.8), new Value<>("", names.length),
                new Value<>("", 0.5),
                null,
                new Value<>("", names),
                null,
                new Value<>("", rootAge));

        TimeTree tree = g.sample().value();

        assertEquals(names.length, tree.getLeafNodes().size(), "leaf count");
        assertEquals(rootAge, tree.getRoot().getAge(), 1e-8, "root age must match provided rootAge");
        List<String> treeNames = Arrays.asList(tree.getTaxaNames());
        for (String name : names) {
            assertTrue(treeNames.contains(name), "taxon '" + name + "' missing from tree");
        }
    }

    // -------------------------------------------------------------------------
    // Subcritical regime (effectiveDeathRate >> lambda)
    // -------------------------------------------------------------------------

    /**
     * When effectiveDeathRate = 1/lifetime >> lambda, the CPP survival CDF saturates
     * well below 1 (asymptote ≈ lambda * lifetime). Previously, the closed-form
     * inverseCDF would blow up to +Inf for p near the asymptote, causing
     * times[0] = conditionAge to become the minimum and crashing coalesce() with
     * IndexOutOfBoundsException: Index -1 out of bounds.
     */
    @RepeatedTest(30)
    void subcriticalRegimeDoesNotCrash() {
        // lambda=0.05, lifetime=0.01 → effectiveDeathRate=100 >> lambda=0.05
        // CDF asymptote ≈ lambda * lifetime = 0.0005 (deeply subcritical)
        CalibratedAgeDependentCPPTree g = generate(
                0.05, 1.0, 5, 0.01, 5.0,
                new Calibration(new String[]{"a", "b", "c"}, 2.0));
        TimeTree tree = assertDoesNotThrow(() -> g.sample().value(),
                "sample() must not throw in subcritical regime");
        assertAgesOrder(tree);
    }
}
