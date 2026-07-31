package calibratedcpp;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.inference.parameter.RealScalarParam;
import calibratedcpp.distribution.Weibull;
import org.apache.commons.math3.distribution.WeibullDistribution;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * A Weibull lifetime with shape k < 1 has an infinite density at t=0 (integrable singularity). The VIDE
 * solver must stay finite there and still produce the correct node-age law. The oracle is an independent
 * renewal-equation solve, G'(t) = birth*( S(t) + integral_0^t S(s) G'(t-s) ds ), which uses only the
 * smooth survival function and never samples the singular density.
 */
public class AgeDependentWeibullSingularityTest {

    private static Tree smallTree() {
        Tree tree = new TreeParser();
        tree.initByName("newick", "((A:1,B:1):1,C:2);", "IsLabelledNewick", true, "adjustTipHeights", true);
        return tree;
    }

    private static CalibratedAgeDependentExtinctionModel weibullModel(double k, double scale,
                                                                      double birth, double rho, double origin) {
        Weibull w = new Weibull();
        w.initByName("shape", new RealScalarParam<>(k, PositiveReal.INSTANCE),
                     "scale", new RealScalarParam<>(scale, PositiveReal.INSTANCE));
        CalibratedAgeDependentExtinctionModel model = new CalibratedAgeDependentExtinctionModel();
        model.initByName("tree", smallTree(),
                "origin",               new RealScalarParam<>(origin, PositiveReal.INSTANCE),
                "birthRate",            new RealScalarParam<>(birth, PositiveReal.INSTANCE),
                "rho",                  new RealScalarParam<>(rho, UnitInterval.INSTANCE),
                "lifetimeDistribution", w,
                "gridSize",             20000);
        model.updateModel();   // triggers preCalc() -> solveVIDE(), populating the splines
        return model;
    }

    /** Reference: m=G' via the renewal equation (smooth kernel S), then G = cumulative-trapezoid(m). */
    private static double[][] reference(double k, double scale, double birth, double rho, double origin, int N) {
        WeibullDistribution life = new WeibullDistribution(k, scale);
        double dt = origin / N;
        double[] S = new double[N + 1];
        for (int j = 0; j <= N; j++) S[j] = 1.0 - life.cumulativeProbability(j * dt);
        double[] m = new double[N + 1];
        m[0] = birth * S[0];
        double denom = 1.0 - 0.5 * birth * dt * S[0];
        for (int i = 1; i <= N; i++) {
            double sum = 0.0;
            for (int j = 1; j <= i - 1; j++) sum += S[j] * m[i - j];
            sum += 0.5 * S[i] * m[0];
            m[i] = (birth * S[i] + birth * dt * sum) / denom;
        }
        double[] G = new double[N + 1];
        for (int i = 1; i <= N; i++) G[i] = G[i - 1] + 0.5 * dt * (m[i - 1] + m[i]);
        return new double[][]{G, m};
    }

    @Test
    public void weibullShapeBelowOneStaysFiniteAndMatchesReference() {
        double scale = 1.0, birth = 1.0, rho = 0.4, origin = 8.0;
        for (double k : new double[]{0.9, 0.7, 0.5, 0.3}) {
            CalibratedAgeDependentExtinctionModel model = weibullModel(k, scale, birth, rho, origin);

            int N = 40000;
            double[][] ref = reference(k, scale, birth, rho, origin, N);
            double[] G = ref[0], m = ref[1];
            double dt = origin / N;

            double maxCdfErr = 0, maxDenErr = 0;
            for (double t = 0.1; t <= origin - 1e-9; t += 0.1) {
                double logCdf = model.calculateLogNodeAgeCDF(t);
                double logDen = model.calculateLogNodeAgeDensity(t);
                assertTrue(Double.isFinite(logCdf), "logCDF non-finite at k=" + k + " t=" + t);
                assertTrue(Double.isFinite(logDen), "logDensity non-finite at k=" + k + " t=" + t);

                int i = (int) Math.floor(t / dt); double f = t / dt - i;
                double Gt = G[i] * (1 - f) + G[i + 1] * f;
                double mt = m[i] * (1 - f) + m[i + 1] * f;
                double refLogCdf = Math.log(rho * Gt) - Math.log(1.0 + rho * Gt);
                double refLogDen = Math.log(rho * mt) - 2.0 * Math.log(1.0 + rho * Gt);
                maxCdfErr = Math.max(maxCdfErr, Math.abs(Math.exp(logCdf) - Math.exp(refLogCdf)));
                maxDenErr = Math.max(maxDenErr, Math.abs(Math.exp(logDen) - Math.exp(refLogDen)));
            }
            System.out.printf("BEAST Weibull k=%.1f  maxCdfAbs=%.3e  maxDenAbs=%.3e%n", k, maxCdfErr, maxDenErr);
            // Accuracy degrades gracefully toward the strong-singularity limit (and the renewal reference
            // itself loses accuracy near t=0 as k -> 0). The guaranteed properties are finiteness (asserted
            // above) and convergence; CDF is well under 5e-3 and the density under ~1% across k in [0.3, 1).
            assertTrue(maxCdfErr < 5e-3, "CDF error too large at k=" + k + ": " + maxCdfErr);
            assertTrue(maxDenErr < 1e-2, "density error too large at k=" + k + ": " + maxDenErr);
        }
    }
}
