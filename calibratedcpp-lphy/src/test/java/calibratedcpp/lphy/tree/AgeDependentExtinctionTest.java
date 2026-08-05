package calibratedcpp.lphy.tree;

import calibratedcpp.lphy.prior.Calibration;
import calibratedcpp.lphy.prior.CalibrationArray;
import lphy.base.evolution.tree.TimeTree;
import lphy.base.evolution.tree.TimeTreeNode;
import lphy.core.model.Value;
import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.WeibullDistribution;
import org.junit.jupiter.api.Test;

import java.util.List;

import static calibratedcpp.lphy.tree.CPPUtils.CDF;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Age-dependent extinction CPP. The Volterra-IDE node-age law is validated against the analytic
 * oracle: an exponential lifetime is memoryless, so the extinction hazard is constant (= 1/mean) and
 * Q(t) must equal the constant-rate CPP CDF. Same-package so the protected {@code resolveRates}/{@code cdf}
 * hooks are reachable.
 */
public class AgeDependentExtinctionTest {

    private static int countLeaves(TimeTreeNode node) {
        List<TimeTreeNode> children = node.getChildren();
        if (children == null || children.isEmpty()) return 1;
        int c = 0;
        for (TimeTreeNode child : children) c += countLeaves(child);
        return c;
    }

    @Test
    public void exponentialLifetimeMatchesConstantRate() {
        double birth = 1.5, mean = 2.0, rho = 0.5;
        double mu = 1.0 / mean;
        double origin = 6.0;

        CalibratedAgeDependentExtinctionTree ad = new CalibratedAgeDependentExtinctionTree(
                new Value<>("", birth),
                new Value<>("", LifetimeModel.of(new ExponentialDistribution(mean))),
                new Value<>("", rho), new Value<>("", 5),
                null, null, /* stemAge */ new Value<>("", origin), null, null);
        ad.resolveRates();   // solves the VIDE once

        double maxErr = 0;
        for (double t = 0.2; t <= origin + 1e-9; t += 0.2) {
            double q = ad.cdf(t);
            double qConst = CDF(birth, mu, rho, t);
            maxErr = Math.max(maxErr, Math.abs(q - qConst));
        }
        // O(h^4) Richardson + cubic spline: the exp-lifetime case matches the analytic CDF to ~5e-11.
        assertTrue(maxErr < 1e-7, "age-dependent(exp) vs constant-rate CDF: max |diff| = " + maxErr);
    }

    @Test
    public void weibullProducesValidTrees() {
        double rootAge = 6.0;
        int n = 8;
        CalibratedAgeDependentExtinctionTree ad = new CalibratedAgeDependentExtinctionTree(
                new Value<>("", 1.5),
                new Value<>("", LifetimeModel.of(new WeibullDistribution(1.5, 2.0))),
                new Value<>("", 1.0), new Value<>("", n),
                null, null, null, /* rootAge */ new Value<>("", rootAge), null);

        for (int i = 0; i < 10; i++) {
            TimeTree t = ad.sample().value();
            assertEquals(n, countLeaves(t.getRoot()), "tip count");
            assertEquals(rootAge, t.getRoot().getAge(), 1e-9, "root age");
        }
    }

    @Test
    public void weibullMixtureHasExactMeanAndValidTrees() {
        double mean = 1.3;
        WeibullMixtureLifetime f = new WeibullMixtureLifetime(
                new Value<>("", mean), new Value<>("", 0.6), new Value<>("", 2.4),
                new Value<>("", new Double[]{0.3, 0.7}));
        LifetimeModel life = f.apply().value();

        assertEquals(mean, life.mean(), 1e-12, "mixture mean is exact by construction of the shared scale");
        assertEquals(1.0, life.survival(0.0), 1e-9, "survival(0) = 1");
        double integral = 0; double dt = 5e-4;
        for (double t = dt / 2; t < 80; t += dt) integral += life.density(t) * dt;
        assertEquals(1.0, integral, 1e-3, "mixture density integrates to 1");

        // shape1 = 0.6 < 1 exercises the t=0 singularity path inside the VIDE solver
        CalibratedAgeDependentExtinctionTree ad = new CalibratedAgeDependentExtinctionTree(
                new Value<>("", 1.5), f.apply(), new Value<>("", 1.0), new Value<>("", 6),
                null, null, null, new Value<>("", 6.0), null);
        for (int i = 0; i < 5; i++)
            assertEquals(6, countLeaves(ad.sample().value().getRoot()), "tip count");
    }

    @Test
    public void weibullShapeBelowOneStaysFinite() {
        // Weibull shape < 1 has an infinite lifetime density at t=0; the VIDE solver must stay finite
        // (previously produced NaN, misreported as "process does not saturate"). A fixed origin isolates
        // the singularity handling from the random-origin saturation heuristic (which is the inference case).
        double birth = 2.0, rootAge = 8.0;
        int n = 10;
        for (double k : new double[]{0.7, 0.5, 0.3, 0.15}) {
            CalibratedAgeDependentExtinctionTree ad = new CalibratedAgeDependentExtinctionTree(
                    new Value<>("", birth),
                    new Value<>("", LifetimeModel.of(new WeibullDistribution(k, 1.0))),
                    new Value<>("", 0.4), new Value<>("", n),
                    null, null, null, /* rootAge */ new Value<>("", rootAge), null);
            ad.resolveRates();   // solve the VIDE once before querying the CDF
            // CDF must be finite and a valid probability across the grid (this is what NaN'd before the fix).
            for (double t = 0.1; t <= rootAge; t += 0.1) {
                double q = ad.cdf(t);
                assertTrue(Double.isFinite(q) && q >= 0.0 && q <= 1.0, "CDF invalid at k=" + k + " t=" + t + ": " + q);
            }
            for (int i = 0; i < 5; i++) {
                TimeTree t = ad.sample().value();
                assertEquals(n, countLeaves(t.getRoot()), "tip count at k=" + k);
                assertEquals(rootAge, t.getRoot().getAge(), 1e-9, "root age at k=" + k);
            }
        }
    }

    @Test
    public void sampledOriginGrowsGridAndStaysFinite() {
        // calibrated, no stemAge/rootAge: the origin is sampled from Q(t)^n, so the grid must grow
        // until Q saturates. Before the fix this would clamp and the inverse-CDF search would run off
        // to ~1e12; here the sampled origin must be finite and reasonable.
        Calibration clade = new Calibration(new String[]{"a", "b", "c"});
        clade.setAge(2.0);
        CalibrationArray cals = new CalibrationArray(new Calibration[]{clade});

        CalibratedAgeDependentExtinctionTree ad = new CalibratedAgeDependentExtinctionTree(
                new Value<>("", 1.5),
                new Value<>("", LifetimeModel.of(new WeibullDistribution(1.5, 2.0))),
                new Value<>("", 1.0), new Value<>("", 6),
                new Value<>("", cals), null, /* stemAge */ null, /* rootAge */ null, null);

        for (int i = 0; i < 5; i++) {
            TimeTree t = ad.sample().value();
            assertEquals(6, countLeaves(t.getRoot()), "tip count");
            double origin = ad.getOrigin().value();
            assertTrue(Double.isFinite(origin) && origin > 2.0 && origin < 100.0,
                    "sampled origin must be finite and reasonable (not clamped garbage): " + origin);
        }
    }
}
