package calibratedcpp.lphy.tree;

import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import org.apache.commons.math3.distribution.WeibullDistribution;

/** Function producing a Weibull {@link LifetimeModel} value for the age-dependent extinction CPP. */
public class WeibullLifetime extends DeterministicFunction<LifetimeModel> {

    public static final String shapeParamName = "shape";
    public static final String scaleParamName = "scale";

    public WeibullLifetime(@ParameterInfo(name = shapeParamName, description = "Weibull shape (k).") Value<Number> shape,
                           @ParameterInfo(name = scaleParamName, description = "Weibull scale (lambda).") Value<Number> scale) {
        setParam(shapeParamName, shape);
        setParam(scaleParamName, scale);
    }

    @GeneratorInfo(name = "weibullLifetime", description = "A Weibull lifetime distribution for the age-dependent extinction CPP.")
    @Override
    public Value<LifetimeModel> apply() {
        double shape = ((Value<Number>) getParams().get(shapeParamName)).value().doubleValue();
        double scale = ((Value<Number>) getParams().get(scaleParamName)).value().doubleValue();
        return new Value<>(null, LifetimeModel.of(new WeibullDistribution(shape, scale)), this);
    }
}
