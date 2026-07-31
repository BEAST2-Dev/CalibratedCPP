package calibratedcpp.lphy.tree;

import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import org.apache.commons.math3.distribution.WeibullDistribution;
import org.apache.commons.math3.special.Gamma;

/**
 * Function producing a two-component Weibull mixture {@link LifetimeDistribution} for the age-dependent
 * extinction CPP. Both components share a common scale θ that is derived from the requested mixture mean,
 * the two shapes and the weights, so the mixture mean is exactly {@code mean}:
 * <pre>  E[X] = θ · Σᵢ wᵢ·Γ(1 + 1/kᵢ)   ⇒   θ = mean / Σᵢ wᵢ·Γ(1 + 1/kᵢ).  </pre>
 * Mirrors the BEAST {@code calibratedcpp.distribution.WeibullMixture}; the LPhyBEAST converter maps this
 * function to that distribution.
 */
public class WeibullMixtureLifetime extends DeterministicFunction<LifetimeDistribution> {

    public static final String meanParamName = "mean";
    public static final String shape1ParamName = "shape1";
    public static final String shape2ParamName = "shape2";
    public static final String weightsParamName = "weights";

    public WeibullMixtureLifetime(
            @ParameterInfo(name = meanParamName, description = "Mean of the mixture.") Value<Number> mean,
            @ParameterInfo(name = shape1ParamName, description = "Shape (k) of the first component.") Value<Number> shape1,
            @ParameterInfo(name = shape2ParamName, description = "Shape (k) of the second component.") Value<Number> shape2,
            @ParameterInfo(name = weightsParamName, description = "Mixture weights; default (0.5, 0.5).", optional = true) Value<Double[]> weights) {
        setParam(meanParamName, mean);
        setParam(shape1ParamName, shape1);
        setParam(shape2ParamName, shape2);
        if (weights != null) setParam(weightsParamName, weights);
    }

    @GeneratorInfo(name = "weibullMixtureLifetime",
            description = "A two-component Weibull mixture lifetime distribution with a specified mean for the "
                    + "age-dependent extinction CPP; the common scale is derived from the mean, shapes and weights.")
    @Override
    public Value<LifetimeDistribution> apply() {
        double mean = ((Value<Number>) getParams().get(meanParamName)).value().doubleValue();
        double k1 = ((Value<Number>) getParams().get(shape1ParamName)).value().doubleValue();
        double k2 = ((Value<Number>) getParams().get(shape2ParamName)).value().doubleValue();

        Value<Double[]> weightsValue = (Value<Double[]>) getParams().get(weightsParamName);
        double w1, w2;
        if (weightsValue == null) {
            w1 = 0.5; w2 = 0.5;
        } else {
            Double[] w = weightsValue.value();
            if (w.length != 2)
                throw new IllegalArgumentException("weibullMixtureLifetime weights must have length 2, got " + w.length);
            double sum = w[0] + w[1];
            w1 = w[0] / sum; w2 = w[1] / sum;   // normalise defensively
        }

        // Shared scale so that the mixture mean equals `mean`.
        double denom = w1 * Gamma.gamma(1.0 + 1.0 / k1) + w2 * Gamma.gamma(1.0 + 1.0 / k2);
        double theta = mean / denom;

        WeibullDistribution d1 = new WeibullDistribution(k1, theta);
        WeibullDistribution d2 = new WeibullDistribution(k2, theta);
        final double ww1 = w1, ww2 = w2, m = mean;

        LifetimeDistribution mixture = new LifetimeDistribution() {
            @Override public double density(double t)  { return ww1 * d1.density(t) + ww2 * d2.density(t); }
            @Override public double survival(double t) {
                return ww1 * (1.0 - d1.cumulativeProbability(t)) + ww2 * (1.0 - d2.cumulativeProbability(t));
            }
            @Override public double mean()             { return m; }   // exact by construction of theta
        };
        return new Value<>(null, mixture, this);
    }
}
