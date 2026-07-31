import calibratedcpp.lphy.tree.CPPTree;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.Value;
import org.junit.jupiter.api.Test;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Random-stem fixed-n mode (CPPTree with n set, no rootAge, randomStemAge=true): the stem is the
 * max of n i.i.d. node ages, so its CDF is Q(t)^n; the n-1 internal depths are drawn from q/Q(stem)
 * and the root sits below the stem. Uses rho=1 so the closed-form stem sampler matches the analytic
 * Q(t)^n (the closed-form path ignores rho — a known pre-existing issue, harmless at rho=1).
 */
public class RandomStemCPPTest {

    private static int countLeaves(TimeTreeNode node) {
        List<TimeTreeNode> children = node.getChildren();
        if (children == null || children.isEmpty()) return 1;
        int c = 0;
        for (TimeTreeNode child : children) c += countLeaves(child);
        return c;
    }

    /** Constant-rate BD node-age CDF Q(t) with sampling prob rho (matches CPPUtils.CDF). */
    private static double Q(double b, double d, double rho, double t) {
        double r = b - d, A = rho * b, B = b * (1 - rho) - d;
        if (Math.abs(r) < 1e-10) return A * t / (1 + A * t);
        double e = Math.exp(-r * t);
        return A * (1 - e) / (A + B * e);
    }

    @Test
    public void randomStemFixedNMatchesMaxOfN() {
        double b = 1.5, d = 0.5, rho = 1.0;
        int n = 6;
        CPPTree cpp = new CPPTree(
                new Value<>("", b), new Value<>("", d), null, null, new Value<>("", rho),
                /* taxa */ null, /* n */ new Value<>("", n),
                /* rootAge */ null, /* randomStemAge */ new Value<>("", true));

        int reps = 4000;
        double[] stems = new double[reps];
        for (int i = 0; i < reps; i++) {
            TimeTree tree = cpp.sample().value();
            assertEquals(n, countLeaves(tree.getRoot()), "fixed n tips");
            double stem = cpp.getConditionAge();               // set to the sampled stem
            double rootAge = tree.getRoot().getAge();
            assertTrue(rootAge <= stem + 1e-9, "root (" + rootAge + ") must sit below stem (" + stem + ")");
            stems[i] = stem;
        }

        // Empirical stem CDF must match Q(t)^n at several points.
        for (double t : new double[]{0.5, 1.0, 2.0, 4.0}) {
            double emp = 0;
            for (double s : stems) if (s <= t) emp++;
            emp /= reps;
            double theo = Math.pow(Q(b, d, rho, t), n);
            assertTrue(Math.abs(emp - theo) < 0.04,
                    "stem CDF at t=" + t + ": empirical " + emp + " vs Q(t)^n " + theo);
        }
    }
}
