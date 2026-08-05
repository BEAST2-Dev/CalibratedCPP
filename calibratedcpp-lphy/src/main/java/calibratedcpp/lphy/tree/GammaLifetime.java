package calibratedcpp.lphy.tree;

import lphy.core.model.DeterministicFunction;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import org.apache.commons.math3.distribution.GammaDistribution;

/** Function producing a Gamma {@link LifetimeModel} value for the age-dependent extinction CPP. */
public class GammaLifetime extends DeterministicFunction<LifetimeModel> {

    public static final String shapeParamName = "shape";
    public static final String scaleParamName = "scale";

    public GammaLifetime(@ParameterInfo(name = shapeParamName, description = "Gamma shape.") Value<Number> shape,
                         @ParameterInfo(name = scaleParamName, description = "Gamma scale.") Value<Number> scale) {
        setParam(shapeParamName, shape);
        setParam(scaleParamName, scale);
    }

    @GeneratorInfo(name = "gammaLifetime", description = "A Gamma lifetime distribution for the age-dependent extinction CPP.")
    @Override
    public Value<LifetimeModel> apply() {
        double shape = ((Value<Number>) getParams().get(shapeParamName)).value().doubleValue();
        double scale = ((Value<Number>) getParams().get(scaleParamName)).value().doubleValue();
        return new Value<>(null, LifetimeModel.of(new GammaDistribution(shape, scale),
                "Gamma(shape=" + shape + ", scale=" + scale + ")"), this);
    }
}
