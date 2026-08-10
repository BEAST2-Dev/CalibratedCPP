package calibratedcpp.lphy.prior;

import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

import java.util.Map;
import java.util.TreeMap;

/**
 * A calibration constraining the MRCA age of a clade to a lower and upper bound, the age being uniform between them.
 * Bounds are how calibrations are specified here, so uniform is the main shape needed.
 * Declaring it as its own LPhy generator, rather than composing one from a generic distribution,
 * is what lets {@code UniformMRCAToBEAST} emit the matching BEAST prior.
 */
public class UniformMRCA implements GenerativeDistribution<Calibration> {

    public static final String taxaParamName = "taxa";
    public static final String upperParamName = "upper";
    public static final String lowerParamName = "lower";

    private Value<String[]> taxa;
    private Value<Number> upper;
    private Value<Number> lower;

    public UniformMRCA(@ParameterInfo(name = taxaParamName, description = "the taxa defining the MRCA node to calibrate") Value<String[]> taxa,
                        @ParameterInfo(name = upperParamName, description = "the upper bound of the uniform prior on the age") Value<Number> upper,
                        @ParameterInfo(name = lowerParamName, description = "the lower bound of the uniform prior on the age") Value<Number> lower) {
        if (lower.value().doubleValue() > upper.value().doubleValue()) {
            throw new IllegalArgumentException("lower (" + lower.value() + ") must be <= upper (" + upper.value() + ")");
        }
        this.taxa = taxa;
        this.upper = upper;
        this.lower = lower;
    }

    @GeneratorInfo(name = "UniformMRCA",
            description = "Creates a calibration constraint for the MRCA of the given taxa, with the age drawn uniformly between the given bounds.")
    @Override
    public RandomVariable<Calibration> sample() {
        double lo = getLower().value().doubleValue();
        double hi = getUpper().value().doubleValue();
        double age = lo + Math.random() * (hi - lo);
        Calibration calibration = new Calibration(getTaxa().value(), age);
        return new RandomVariable<>(null, calibration, this);
    }

    @Override
    public Map<String, Value> getParams() {
        return new TreeMap<>() {{
            put(taxaParamName, taxa);
            put(upperParamName, upper);
            put(lowerParamName, lower);
        }};
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(taxaParamName)) {
            this.taxa = value;
        } else if (paramName.equals(upperParamName)) {
            this.upper = value;
        } else if (paramName.equals(lowerParamName)) {
            this.lower = value;
        } else {
            throw new RuntimeException("Unrecognised parameter name: " + paramName);
        }
    }

    public Value<String[]> getTaxa() {
        return getParams().get(taxaParamName);
    }

    public Value<Number> getUpper() {
        return getParams().get(upperParamName);
    }

    public Value<Number> getLower() {
        return getParams().get(lowerParamName);
    }
}
