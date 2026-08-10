package calibratedcpp.lphy.prior;

import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

import java.util.Map;
import java.util.TreeMap;

/**
 * A calibration whose MRCA age is offset-exponentially distributed: age = offset + Exp(mean), a hard
 * minimum bound with exponential decay above it. Used for deep, well-preserved fossil calibrations,
 * as an alternative to {@link UniformMRCA}'s hard upper bound. Declared as its own LPhy generator
 * for the same reason as {@link UniformMRCA}.
 */
public class OffsetExponentialMRCA implements GenerativeDistribution<Calibration> {

    public static final String taxaParamName = "taxa";
    public static final String offsetParamName = "offset";
    public static final String meanParamName = "mean";

    private Value<String[]> taxa;
    private Value<Number> offset;
    private Value<Number> mean;

    public OffsetExponentialMRCA(@ParameterInfo(name = taxaParamName, description = "the taxa defining the MRCA node to calibrate") Value<String[]> taxa,
                                  @ParameterInfo(name = offsetParamName, description = "the hard minimum bound (offset) of the age") Value<Number> offset,
                                  @ParameterInfo(name = meanParamName, description = "the mean of the exponential distribution, measured from the offset") Value<Number> mean) {
        if (mean.value().doubleValue() <= 0) {
            throw new IllegalArgumentException("mean (" + mean.value() + ") must be > 0");
        }
        this.taxa = taxa;
        this.offset = offset;
        this.mean = mean;
    }

    @GeneratorInfo(name = "OffsetExponentialMRCA",
            description = "Creates a calibration constraint for the MRCA of the given taxa, with the age drawn from " +
                    "an exponential distribution of the given mean, offset by a hard minimum bound.")
    @Override
    public RandomVariable<Calibration> sample() {
        double off = getOffset().value().doubleValue();
        double m = getMean().value().doubleValue();
        double age = off - m * Math.log(1.0 - Math.random());
        Calibration calibration = new Calibration(getTaxa().value(), age);
        return new RandomVariable<>(null, calibration, this);
    }

    @Override
    public Map<String, Value> getParams() {
        return new TreeMap<>() {{
            put(taxaParamName, taxa);
            put(offsetParamName, offset);
            put(meanParamName, mean);
        }};
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(taxaParamName)) {
            this.taxa = value;
        } else if (paramName.equals(offsetParamName)) {
            this.offset = value;
        } else if (paramName.equals(meanParamName)) {
            this.mean = value;
        } else {
            throw new RuntimeException("Unrecognised parameter name: " + paramName);
        }
    }

    public Value<String[]> getTaxa() {
        return getParams().get(taxaParamName);
    }

    public Value<Number> getOffset() {
        return getParams().get(offsetParamName);
    }

    public Value<Number> getMean() {
        return getParams().get(meanParamName);
    }
}
