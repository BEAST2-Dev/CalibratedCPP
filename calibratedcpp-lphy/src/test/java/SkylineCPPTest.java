import calibratedcpp.lphy.tree.CalibratedBirthDeathSkylineTree;
import calibratedcpp.lphy.tree.CalibratedCPPTree;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.Value;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Skyline CPP. A single interval must reduce exactly to the constant-rate CalibratedCPP (proved
 * algebraically: Q = ρλ(e^{rt}-1)/[(b-d)+ρλ(e^{rt}-1)] both ways), so their node-age distributions
 * match; a two-interval skyline must still produce valid trees.
 */
public class SkylineCPPTest {

    private static void internalAges(TimeTreeNode node, List<Double> ages) {
        List<TimeTreeNode> children = node.getChildren();
        if (children == null || children.isEmpty()) return;
        if (node.getParent() != null) ages.add(node.getAge());   // non-root internal node
        for (TimeTreeNode child : children) internalAges(child, ages);
    }

    private static int countLeaves(TimeTreeNode node) {
        List<TimeTreeNode> children = node.getChildren();
        if (children == null || children.isEmpty()) return 1;
        int c = 0;
        for (TimeTreeNode child : children) c += countLeaves(child);
        return c;
    }

    private static double ksD(double[] a, double[] b) {
        double[] x = a.clone(), y = b.clone();
        Arrays.sort(x); Arrays.sort(y);
        int i = 0, j = 0;
        double d = 0;
        while (i < x.length && j < y.length) {
            d = Math.max(d, Math.abs((double) i / x.length - (double) j / y.length));
            if (x[i] <= y[j]) i++; else j++;
        }
        return d;
    }

    private static Value<Double[]> arr(double... v) {
        Double[] a = new Double[v.length];
        for (int i = 0; i < v.length; i++) a[i] = v[i];
        return new Value<>("", a);
    }

    @Test
    public void singleIntervalMatchesConstantRate() {
        double b = 2.0, d = 0.5, rho = 1.0;
        int n = 10;
        double rootAge = 6.0;

        CalibratedBirthDeathSkylineTree sky = new CalibratedBirthDeathSkylineTree(
                arr(b), arr(d), null, null, null, /* changeTimes */ null, new Value<>("", rho),
                new Value<>("", n), null, null, null, new Value<>("", rootAge));
        CalibratedCPPTree cst = new CalibratedCPPTree(
                new Value<>("", b), new Value<>("", d), null, null, new Value<>("", rho),
                new Value<>("", n), null, null, null, new Value<>("", rootAge));

        List<Double> skyAges = new ArrayList<>(), cstAges = new ArrayList<>();
        int reps = 800;
        for (int i = 0; i < reps; i++) {
            TimeTree ts = sky.sample().value();
            assertEquals(n, countLeaves(ts.getRoot()), "skyline fixed-root tip count");
            assertEquals(rootAge, ts.getRoot().getAge(), 1e-9, "skyline root age");
            internalAges(ts.getRoot(), skyAges);
            internalAges(cst.sample().value().getRoot(), cstAges);
        }

        double[] x = skyAges.stream().mapToDouble(Double::doubleValue).toArray();
        double[] y = cstAges.stream().mapToDouble(Double::doubleValue).toArray();
        double d0 = ksD(x, y);
        double neff = (double) x.length * y.length / (x.length + y.length);
        double crit = Math.sqrt(-0.5 * Math.log(1e-4 / 2)) / Math.sqrt(neff);   // family-wise ~1e-4
        assertTrue(d0 < crit, "1-interval skyline vs constant node ages: D=" + d0 + " crit=" + crit);
    }

    @Test
    public void subcriticalFixedStemHasFiniteBranchLengths() {
        // Subcritical (R0 = lambda/mu = 0.5 < 1): Q saturates below 1, so a sampled max-of-n stem is
        // infinite with positive probability. A supplied stemAge must condition the origin to a finite
        // value; before the fix this path ignored stemAge, sampled the stem, and produced NaN branches.
        int n = 10;
        double stem = 5.0;
        CalibratedBirthDeathSkylineTree sky = new CalibratedBirthDeathSkylineTree(
                /* birthRate */ null, arr(1.0), null, null, /* reproductiveNumber */ arr(0.5),
                /* changeTimes */ null, new Value<>("", 0.5),
                new Value<>("", n), null, null, /* stemAge */ new Value<>("", stem), null);

        for (int i = 0; i < 50; i++) {
            TimeTree t = sky.sample().value();
            assertEquals(n, countLeaves(t.getRoot()), "tip count");
            assertBranchLengthsFinite(t.getRoot());
            assertTrue(t.getRoot().getAge() < stem + 1e-9, "root must sit below the stem: " + t.getRoot().getAge());
        }
    }

    private static void assertBranchLengthsFinite(TimeTreeNode node) {
        if (node.getParent() != null) {
            double bl = node.getParent().getAge() - node.getAge();
            assertTrue(Double.isFinite(bl), "branch length must be finite, got " + bl);
        }
        if (node.getChildren() != null) for (TimeTreeNode c : node.getChildren()) assertBranchLengthsFinite(c);
    }

    @Test
    public void twoIntervalProducesValidTrees() {
        double rho = 1.0;
        int n = 8;
        double rootAge = 6.0;
        CalibratedBirthDeathSkylineTree sky = new CalibratedBirthDeathSkylineTree(
                arr(2.0, 1.0), arr(0.5, 0.5), null, null, null, arr(3.0), new Value<>("", rho),
                new Value<>("", n), null, null, null, new Value<>("", rootAge));

        for (int i = 0; i < 500; i++) {
            TimeTree t = sky.sample().value();
            assertEquals(n, countLeaves(t.getRoot()), "two-interval tip count");
            assertEquals(rootAge, t.getRoot().getAge(), 1e-9, "two-interval root age");
            List<Double> ages = new ArrayList<>();
            internalAges(t.getRoot(), ages);
            for (double a : ages) assertTrue(a > 0 && a <= rootAge + 1e-9, "internal age in (0, root]: " + a);
        }
    }
}
