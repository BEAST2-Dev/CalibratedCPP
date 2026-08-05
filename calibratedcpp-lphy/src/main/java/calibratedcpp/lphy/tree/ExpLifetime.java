package calibratedcpp.lphy.tree;

import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import org.apache.commons.math3.distribution.ExponentialDistribution;

/**
 * Function producing an Exponential {@link LifetimeModel} value. Because exponential lifetimes
 * are memoryless, this reduces the age-dependent extinction CPP to the constant-rate CPP with death
 * rate 1/mean — useful as a validation baseline.
 */
public class ExpLifetime extends DeterministicFunction<LifetimeModel> {

    public static final String meanParamName = "mean";

    public ExpLifetime(@ParameterInfo(name = meanParamName, description = "mean lifetime (= 1 / death rate).") Value<Number> mean) {
        setParam(meanParamName, mean);
    }

    @GeneratorInfo(name = "expLifetime", description = "An exponential lifetime distribution (memoryless; reduces to constant-rate death).")
    @Override
    public Value<LifetimeModel> apply() {
        double mean = ((Value<Number>) getParams().get(meanParamName)).value().doubleValue();
        return new Value<>(null, LifetimeModel.of(new ExponentialDistribution(mean),
                "Exp(mean=" + mean + ")"), this);
    }
}
