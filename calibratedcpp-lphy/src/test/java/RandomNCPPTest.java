import calibratedcpp.lphy.tree.CalibratedCPPTree;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.Value;
import org.junit.jupiter.api.Test;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Random-N mode: a CalibratedCPPTree with n == null and no calibrations grows a random number
 * of tips from the process, conditioned on the origin (rootAge). Verifies the tip count actually
 * varies, its mean is near the theoretical Q(origin)/(1-Q(origin)) + 1, and every tree is valid.
 */
public class RandomNCPPTest {

    private static int countLeaves(TimeTreeNode node) {
        List<TimeTreeNode> children = node.getChildren();
        if (children == null || children.isEmpty()) return 1;
        int c = 0;
        for (TimeTreeNode child : children) c += countLeaves(child);
        return c;
    }

    private static double qOrigin(double b, double d, double origin) {
        // Q(origin) for constant-rate BD with rho=1: (1 - e^{-r*origin}) / (1 - (d/b)*e^{-r*origin})
        double r = b - d;
        return (1 - Math.exp(-r * origin)) / (1 - (d / b) * Math.exp(-r * origin));
    }

    private static CalibratedCPPTree cpp(double b, double d, double rho, Double stemAge, Double rootAge) {
        return new CalibratedCPPTree(
                new Value<>("", b), new Value<>("", d), null, null, new Value<>("", rho),
                /* n */ null, /* calibrations */ null, /* otherNames */ null,
                stemAge == null ? null : new Value<>("", stemAge),
                rootAge == null ? null : new Value<>("", rootAge));
    }

    /**
     * Fixed STEM, random N: one geometric run, so E[tips] = Q/(1-Q) + 1, and the root (deepest depth)
     * sits strictly below the stem.
     */
    @Test
    public void randomNFixedStem() {
        double b = 1.0, d = 0.5, rho = 1.0, stem = 2.0;
        double q = qOrigin(b, d, stem);
        double expectedMeanTips = q / (1 - q) + 1;                 // ~4.4
        CalibratedCPPTree cpp = cpp(b, d, rho, /* stemAge */ stem, /* rootAge */ null);

        int reps = 2000;
        long totalTips = 0;
        Set<Integer> distinct = new HashSet<>();
        for (int i = 0; i < reps; i++) {
            TimeTree tree = cpp.sample().value();
            distinct.add(countLeaves(tree.getRoot()));
            totalTips += countLeaves(tree.getRoot());
            assertTrue(tree.getRoot().getAge() < stem + 1e-9, "root must sit below the stem");
        }
        double meanTips = totalTips / (double) reps;
        assertTrue(distinct.size() > 1, "tip count should vary across replicates");
        assertTrue(Math.abs(meanTips - expectedMeanTips) < 0.7,
                "stem: mean tips " + meanTips + " should be near " + expectedMeanTips);
    }

    /**
     * Fixed ROOT, random N: the root is the crown of two independent sub-clades (two geometric runs with
     * the root placed at the first exceedance), so E[tips] = 2·Q/(1-Q) + 2, and the root age equals the
     * fixed origin exactly.
     */
    @Test
    public void randomNFixedRoot() {
        double b = 1.0, d = 0.5, rho = 1.0, root = 2.0;
        double q = qOrigin(b, d, root);
        double expectedMeanTips = 2.0 * q / (1 - q) + 2.0;        // ~8.9
        CalibratedCPPTree cpp = cpp(b, d, rho, /* stemAge */ null, /* rootAge */ root);

        int reps = 2000;
        long totalTips = 0;
        Set<Integer> distinct = new HashSet<>();
        for (int i = 0; i < reps; i++) {
            TimeTree tree = cpp.sample().value();
            distinct.add(countLeaves(tree.getRoot()));
            totalTips += countLeaves(tree.getRoot());
            assertEquals(root, tree.getRoot().getAge(), 1e-9, "root age must equal the fixed origin");
        }
        double meanTips = totalTips / (double) reps;
        assertTrue(distinct.size() > 1, "tip count should vary across replicates");
        assertTrue(Math.abs(meanTips - expectedMeanTips) < 1.0,
                "root: mean tips " + meanTips + " should be near " + expectedMeanTips);
    }
}
