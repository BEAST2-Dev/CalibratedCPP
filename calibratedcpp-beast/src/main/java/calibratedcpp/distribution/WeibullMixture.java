package calibratedcpp.distribution;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.distribution.ScalarDistribution;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.SimplexParam;
import beast.base.spec.type.RealScalar;
import beast.base.spec.type.Simplex;
import org.apache.commons.math3.special.Gamma;

import java.util.List;

@Description("Two-component Weibull mixture with a specified mean; the common scale is "
        + "derived from the mean, the shapes and the weights.")
public class WeibullMixture extends ScalarDistribution<RealScalar<PositiveReal>, Double> {

    public Input<RealScalar<PositiveReal>> meanInput =
            new Input<>("mean", "Mean of the mixture.", Input.Validate.REQUIRED);
    public Input<RealScalar<PositiveReal>> shape1Input =
            new Input<>("shape1", "Shape (k) of the first component.", Input.Validate.REQUIRED);
    public Input<RealScalar<PositiveReal>> shape2Input =
            new Input<>("shape2", "Shape (k) of the second component.", Input.Validate.REQUIRED);
    public Input<Simplex> weightsInput =
            new Input<>("weights", "Mixture weights; defaults to equal (0.5, 0.5).");

    // One shared, derived scale instance drives both components.
    private final RealScalarParam<PositiveReal> scale =
            new RealScalarParam<>(1.0, PositiveReal.INSTANCE);
    private final Weibull weibull1 = new Weibull();
    private final Weibull weibull2 = new Weibull();
    private final ScalarMixtureDistribution<RealScalar<PositiveReal>, Double> mixture =
            new ScalarMixtureDistribution<>();

    private Simplex weights;

    // Last parameter values the derived scale/components were built for, so refresh() is idempotent:
    // it early-returns when nothing changed, making it cheap to call from every accessor (the VIDE
    // solver hits density()/cumulativeProbability() thousands of times per solve).
    private double lastMean = Double.NaN, lastK1 = Double.NaN, lastK2 = Double.NaN, lastW0 = Double.NaN;

    @Override
    public void initAndValidate() {
        weights = weightsInput.get() != null
                ? weightsInput.get()
                : new SimplexParam(new double[]{0.5, 0.5});

        weibull1.initByName("shape", shape1Input.get(), "scale", scale);
        weibull2.initByName("shape", shape2Input.get(), "scale", scale);
        mixture.initByName(
                "distribution", weibull1,
                "distribution", weibull2,
                "weights", weights);

        super.initAndValidate();
        refresh();
    }

    @Override
    public void refresh() {
        final double mean = meanInput.get().get();
        final double k1 = shape1Input.get().get();
        final double k2 = shape2Input.get().get();
        final double w0 = weights.get(0);

        // Idempotent: the derived scale/components already match the current parameters, so nothing
        // to rebuild. This makes refresh() free to call from every accessor (below).
        if (mean == lastMean && k1 == lastK1 && k2 == lastK2 && w0 == lastW0) return;

        // E[X] = theta * sum_i w_i * Gamma(1 + 1/k_i)  =>  theta = mean / that sum
        final double denom = w0 * Gamma.gamma(1.0 + 1.0 / k1)
                + weights.get(1) * Gamma.gamma(1.0 + 1.0 / k2);
        scale.set(mean / denom);

        mixture.refresh();   // cascades to weibull1.refresh() and weibull2.refresh()
        lastMean = mean; lastK1 = k1; lastK2 = k2; lastW0 = w0;
    }

    // --- delegate the scalar surface to the internal mixture ---
    // Every accessor refreshes first so the mixture always reflects the current parameters; the guard
    // in refresh() makes this cheap. Without it, a consumer that reads the distribution without going
    // through calculateLogP() (e.g. the age-dependent likelihood's VIDE solver) scores a stale Q.

    @Override public double logDensity(double x)            { refresh(); return mixture.logDensity(x); }
    @Override public double density(double x)               { refresh(); return mixture.density(x); }
    @Override public double cumulativeProbability(double x) { refresh(); return mixture.cumulativeProbability(x); }

    @Override protected double calcLogP(Double value)       { refresh(); return mixture.logDensity(value); }

    @Override
    public double calculateLogP() {
        refresh();
        logP = calcLogP(param.get());
        return logP;
    }

    @Override
    public double getMean(){
        return meanInput.get().get();
    }

    @Override public List<Double> sample() { refresh(); return mixture.sample(); }

    @Override public Double getLowerBoundOfParameter() { return mixture.getLowerBoundOfParameter(); }
    @Override public Double getUpperBoundOfParameter() { return mixture.getUpperBoundOfParameter(); }
}