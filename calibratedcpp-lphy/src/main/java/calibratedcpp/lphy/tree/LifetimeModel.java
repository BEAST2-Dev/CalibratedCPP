package calibratedcpp.lphy.tree;

import org.apache.commons.math3.distribution.RealDistribution;

/**
 * Engine-neutral lifetime law for the age-dependent-extinction CPP: an individual's lifetime
 * has this density, and the node-age law is derived from it via a Volterra IDE (see
 * {@link CalibratedAgeDependentExtinctionTree}). Deliberately independent of both LPhy and BEAST — it
 * exposes only what the solver needs ({@code density}, {@code survival}), plus a mean for the horizon
 * heuristic. LPhy functions (e.g. {@code weibullLifetime}) produce these as {@code Value}s; the
 * LPhyBEAST converter maps them to a BEAST {@code ScalarDistribution}.
 */
public interface LifetimeModel {

    /** Lifetime probability density g(t). */
    double density(double t);

    /** Lifetime survival S(t) = 1 - CDF(t). */
    double survival(double t);

    /** Mean lifetime (used only to size the solver grid when no origin is supplied). */
    double mean();

    static LifetimeModel of(RealDistribution d) {
        return of(d, d.getClass().getSimpleName());
    }

    /** Wrap any Apache commons-math3 continuous distribution; {@code label} is what loggers print. */
    static LifetimeModel of(RealDistribution d, String label) {
        return new LifetimeModel() {
            @Override public double density(double t)  { return d.density(t); }
            @Override public double survival(double t) { return 1.0 - d.cumulativeProbability(t); }
            @Override public double mean()             { return d.getNumericalMean(); }
            @Override public String toString()         { return label; }
        };
    }
}
