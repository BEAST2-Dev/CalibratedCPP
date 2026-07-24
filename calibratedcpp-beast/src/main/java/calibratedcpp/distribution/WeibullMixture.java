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

        // E[X] = theta * sum_i w_i * Gamma(1 + 1/k_i)  =>  theta = mean / that sum
        final double denom = weights.get(0) * Gamma.gamma(1.0 + 1.0 / k1)
                + weights.get(1) * Gamma.gamma(1.0 + 1.0 / k2);
        scale.set(mean / denom);

        mixture.refresh();   // cascades to weibull1.refresh() and weibull2.refresh()
    }

    // --- delegate the scalar surface to the internal mixture ---

    @Override public double logDensity(double x)            { return mixture.logDensity(x); }
    @Override public double density(double x)               { return mixture.density(x); }
    @Override public double cumulativeProbability(double x) { return mixture.cumulativeProbability(x); }

    @Override protected double calcLogP(Double value)       { return mixture.logDensity(value); }

    @Override
    public double calculateLogP() {
        refresh();
        logP = calcLogP(param.get());
        return logP;
    }

    @Override public List<Double> sample() { return mixture.sample(); }

    @Override public Double getLowerBoundOfParameter() { return mixture.getLowerBoundOfParameter(); }
    @Override public Double getUpperBoundOfParameter() { return mixture.getUpperBoundOfParameter(); }
}