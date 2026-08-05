package calibratedcpp.lphy.tree;

import calibratedcpp.lphy.prior.CalibrationArray;
import lphy.base.distribution.DistributionConstants;
import lphy.base.evolution.birthdeath.BirthDeathConstants;
import lphy.base.evolution.tree.TimeTree;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

import java.util.Arrays;
import java.util.Map;

/**
 * Calibrated CPP whose node ages follow a birth-death SKYLINE law: piecewise-constant birth and
 * death rates with incomplete extant sampling. It supplies {@link #logCDF} and {@link #invertCDF}.
 * The log-CDF is ported directly from the BEAST
 * {@code CalibratedBirthDeathSkylineModel}, and a single interval reduces exactly to the constant-rate.
 *
 * <p>Change times are absolute ages, strictly increasing; rate arrays are ordered present-to-past
 * with one more entry than the change times (values[j] applies to interval [t_j, t_{j+1}]), matching
 * BEAST's {@code timesAreAges=true} convention. Relative / root-to-present orderings are inference-only.
 */
public class CalibratedBirthDeathSkylineTree extends AbstractCalibratedCPPTree {

    public static final String reproductiveNumberName = "reproductiveNumber";
    public static final String changeTimesName = "changeTimes";

    Value<Double[]> birthRate, deathRate, diversification, turnover, reproductiveNumber;
    Value<Double[]> changeTimes;

    private double rhoVal;
    private double[] intervalStartTimes, lambda, r, cumulativeIntegral, cumulativeExpR;

    public CalibratedBirthDeathSkylineTree(
            @ParameterInfo(name = BirthDeathConstants.lambdaParamName, description = "per-interval birth rate (present-to-past).", optional = true) Value<Double[]> birthRate,
            @ParameterInfo(name = BirthDeathConstants.muParamName, description = "per-interval death rate.", optional = true) Value<Double[]> deathRate,
            @ParameterInfo(name = BirthDeathConstants.diversificationParamName, description = "per-interval diversification (lambda - mu).", optional = true) Value<Double[]> diversification,
            @ParameterInfo(name = BirthDeathConstants.turnoverParamName, description = "per-interval turnover (mu/lambda).", optional = true) Value<Double[]> turnover,
            @ParameterInfo(name = reproductiveNumberName, description = "per-interval reproductive number (lambda/mu).", optional = true) Value<Double[]> reproductiveNumber,
            @ParameterInfo(name = changeTimesName, description = "rate change times as absolute ages, strictly increasing; length = (#rate values - 1).", optional = true) Value<Double[]> changeTimes,
            @ParameterInfo(name = BirthDeathConstants.rhoParamName, description = "sampling probability") Value<Number> rho,
            @ParameterInfo(name = DistributionConstants.nParamName, description = "the total number of taxa; omit for a random number of tips (uncalibrated only).", optional = true) Value<Integer> n,
            @ParameterInfo(name = calibrationsName, description = "an array of calibrations generated from a MRCA prior.", optional = true) Value<CalibrationArray> calibrations,
            @ParameterInfo(name = otherTaxaNames, description = "a string array of taxa names for non-calibrated tips.", optional = true) Value<String[]> otherNames,
            @ParameterInfo(name = stemAgeName, description = "the stem age working as condition time.", optional = true) Value<Number> stemAge,
            @ParameterInfo(name = rootAgeName, description = "the root age to condition on when no calibrations are provided.", optional = true) Value<Number> rootAge) {
        super(n, rho, calibrations, otherNames, stemAge, rootAge);

        int count = 0;
        if (birthRate != null) count++;
        if (deathRate != null) count++;
        if (diversification != null) count++;
        if (turnover != null) count++;
        if (reproductiveNumber != null) count++;
        if (count != 2) {
            throw new IllegalArgumentException(
                    "Must specify exactly two of: birthRate, deathRate, diversification, turnover, reproductiveNumber.");
        }
        if (reproductiveNumber != null && turnover != null) {
            throw new IllegalArgumentException("Cannot specify both reproductiveNumber and turnover.");
        }

        this.birthRate = birthRate;
        this.deathRate = deathRate;
        this.diversification = diversification;
        this.turnover = turnover;
        this.reproductiveNumber = reproductiveNumber;
        this.changeTimes = changeTimes;
    }

    @GeneratorInfo(name = "CalibratedBirthDeathSkylineTree", examples = {},
            description = "The Calibrated Coalescent Point Process with piecewise-constant (skyline) birth and "
                    + "death rates. Node ages follow the BDSKY node-age law; a single interval matches the "
                    + "constant-rate CalibratedCPP.")
    @Override
    public RandomVariable<TimeTree> sample() {
        return super.sample();
    }

    @Override
    protected void resolveRates() {
        rhoVal = getSamplingProb().value().doubleValue();

        Double[] ct = (changeTimes != null && changeTimes.value() != null) ? changeTimes.value() : new Double[0];
        int nIntervals = ct.length + 1;

        checkLength(birthRate, nIntervals, BirthDeathConstants.lambdaParamName);
        checkLength(deathRate, nIntervals, BirthDeathConstants.muParamName);
        checkLength(diversification, nIntervals, BirthDeathConstants.diversificationParamName);
        checkLength(turnover, nIntervals, BirthDeathConstants.turnoverParamName);
        checkLength(reproductiveNumber, nIntervals, reproductiveNumberName);

        intervalStartTimes = new double[nIntervals];
        for (int i = 0; i < ct.length; i++) {
            intervalStartTimes[i + 1] = ct[i];
            if (intervalStartTimes[i + 1] <= intervalStartTimes[i]) {
                throw new IllegalArgumentException("changeTimes must be strictly increasing positive ages.");
            }
        }

        lambda = new double[nIntervals];
        r = new double[nIntervals];
        cumulativeIntegral = new double[nIntervals];
        cumulativeExpR = new double[nIntervals];
        cumulativeIntegral[0] = Double.NEGATIVE_INFINITY;
        cumulativeExpR[0] = 0.0;

        double logRunningSum = Double.NEGATIVE_INFINITY;
        double rRunningSum = 0.0;

        for (int j = 0; j < nIntervals; j++) {
            double[] lm = resolveInterval(j);
            lambda[j] = lm[0];
            r[j] = lm[0] - lm[1];
            if (j < nIntervals - 1) {
                double dt = intervalStartTimes[j + 1] - intervalStartTimes[j];
                logRunningSum = logSumExp(logRunningSum, calculateSegment(lambda[j], r[j], dt, rRunningSum));
                rRunningSum += r[j] * dt;
                cumulativeIntegral[j + 1] = logRunningSum;
                cumulativeExpR[j + 1] = rRunningSum;
            }
        }
    }

    /** Resolve (lambda, mu) for interval j from the two specified rate arrays (mirrors the BEAST model). */
    private double[] resolveInterval(int j) {
        double vL = at(birthRate, j), vM = at(deathRate, j), vD = at(diversification, j),
                vT = at(turnover, j), vR = at(reproductiveNumber, j);
        double l, m;
        if (birthRate != null && deathRate != null) { l = vL; m = vM; }
        else if (birthRate != null && diversification != null) { l = vL; m = vL - vD; }
        else if (birthRate != null && reproductiveNumber != null) { l = vL; m = vL / vR; }
        else if (birthRate != null && turnover != null) { l = vL; m = vL * vT; }
        else if (deathRate != null && diversification != null) { m = vM; l = vM + vD; }
        else if (deathRate != null && turnover != null) { m = vM; l = vM / vT; }
        else if (deathRate != null && reproductiveNumber != null) { m = vM; l = m * vR; }
        else if (diversification != null && reproductiveNumber != null) { m = vD / (vR - 1.0); l = m * vR; }
        else if (diversification != null && turnover != null) { l = vD / (1.0 - vT); m = l * vT; }
        else throw new IllegalArgumentException("Invalid skyline rate parameter combination.");
        if (l < 0.0 || m < 0.0) {
            throw new IllegalArgumentException("Negative birth or death rate in interval " + j + " (l=" + l + ", m=" + m + ").");
        }
        return new double[]{l, m};
    }

    private static double at(Value<Double[]> arr, int j) {
        return (arr != null && arr.value() != null) ? arr.value()[j] : 0.0;
    }

    private static void checkLength(Value<Double[]> arr, int nIntervals, String name) {
        if (arr != null && arr.value() != null && arr.value().length != nIntervals) {
            throw new IllegalArgumentException(name + " has " + arr.value().length + " values but expected "
                    + nIntervals + " (changeTimes length + 1).");
        }
    }

    @Override
    protected double logCDF(double t) {
        int m = getInterval(t);
        double logInt = logSumExp(cumulativeIntegral[m], calculateSegment(lambda[m], r[m], t - intervalStartTimes[m], cumulativeExpR[m]));
        return logInt - logSumExp(0.0, logInt);
    }

    /**
     * Closed-form inverse, replacing the base's numerical bisection. Since Q = I/(1+I), a target CDF
     * value p corresponds to accumulated integral I = p/(1-p); locate the interval whose integral
     * bracket [I(t_m), I(t_{m+1})) contains it, then invert that interval's constant-rate integral
     * exactly. Reduces to {@code CPPUtils.inverseCDF} for a single interval.
     */
    @Override
    protected double invertCDF(double p) {
        if (p <= 0.0) return 0.0;
        if (p >= 1.0) return Double.POSITIVE_INFINITY;
        double target = p / (1.0 - p);                       // Q = I/(1+I)  =>  I = Q/(1-Q)

        // interval m whose accumulated-integral bracket contains the target
        int m = intervalStartTimes.length - 1;
        for (int j = 1; j < intervalStartTimes.length; j++) {
            if (Math.exp(cumulativeIntegral[j]) > target) { m = j - 1; break; }
        }

        double I0 = Math.exp(cumulativeIntegral[m]);          // integral at the start of interval m (0 when m==0)
        double base = rhoVal * lambda[m] * Math.exp(cumulativeExpR[m]);
        double dt;
        if (Math.abs(r[m]) < 1e-9) {
            dt = (target - I0) / base;
        } else {
            double arg = (target - I0) * r[m] / base;
            if (arg <= -1.0) return Double.POSITIVE_INFINITY; // beyond this interval's reachable Q (subcritical saturation)
            dt = Math.log1p(arg) / r[m];
        }
        return intervalStartTimes[m] + dt;
    }

    // --- math helpers ported from CalibratedBirthDeathSkylineModel ---

    private double calculateSegment(double l, double r_val, double dt, double expOffset) {
        double logTerm;
        double x = r_val * dt;
        if (Math.abs(r_val) < 1e-9) {
            logTerm = Math.log(rhoVal) + Math.log(l) + Math.log(dt);
        } else if (r_val > 0) {
            logTerm = Math.log(rhoVal) + Math.log(l) - Math.log(r_val) + Math.log(Math.expm1(x));
        } else {
            logTerm = Math.log(rhoVal) + Math.log(l) - Math.log(-r_val) + Math.log(-Math.expm1(x));
        }
        return logTerm + expOffset;
    }

    private int getInterval(double t) {
        int i = Arrays.binarySearch(intervalStartTimes, t);
        return i < 0 ? Math.max(0, -i - 2) : i;
    }

    private static double logSumExp(double a, double b) {
        if (a == Double.NEGATIVE_INFINITY) return b;
        if (b == Double.NEGATIVE_INFINITY) return a;
        return Math.max(a, b) + Math.log1p(Math.exp(-Math.abs(a - b)));
    }

    @Override
    protected AbstractCalibratedCPPTree newSubClade(int nTaxa, CalibrationArray subCalibrations) {
        return new CalibratedBirthDeathSkylineTree(birthRate, deathRate, diversification, turnover,
                reproductiveNumber, changeTimes, getSamplingProb(), new Value<>("n", nTaxa),
                new Value<>("", subCalibrations), null, null, null);
    }

    @Override
    public Map<String, Value> getParams() {
        Map<String, Value> map = super.getParams();
        if (birthRate != null) map.put(BirthDeathConstants.lambdaParamName, birthRate);
        if (deathRate != null) map.put(BirthDeathConstants.muParamName, deathRate);
        if (diversification != null) map.put(BirthDeathConstants.diversificationParamName, diversification);
        if (turnover != null) map.put(BirthDeathConstants.turnoverParamName, turnover);
        if (reproductiveNumber != null) map.put(reproductiveNumberName, reproductiveNumber);
        if (changeTimes != null) map.put(changeTimesName, changeTimes);
        return map;
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(BirthDeathConstants.lambdaParamName)) birthRate = value;
        else if (paramName.equals(BirthDeathConstants.muParamName)) deathRate = value;
        else if (paramName.equals(BirthDeathConstants.diversificationParamName)) diversification = value;
        else if (paramName.equals(BirthDeathConstants.turnoverParamName)) turnover = value;
        else if (paramName.equals(reproductiveNumberName)) reproductiveNumber = value;
        else if (paramName.equals(changeTimesName)) changeTimes = value;
        else super.setParam(paramName, value);
    }

    public Value<Double[]> getBirthRate()           { return getParams().get(BirthDeathConstants.lambdaParamName); }
    public Value<Double[]> getDeathRate()           { return getParams().get(BirthDeathConstants.muParamName); }
    public Value<Double[]> getDiversificationRate() { return getParams().get(BirthDeathConstants.diversificationParamName); }
    public Value<Double[]> getTurnover()            { return getParams().get(BirthDeathConstants.turnoverParamName); }
    public Value<Double[]> getReproductiveNumber()  { return getParams().get(reproductiveNumberName); }
    public Value<Double[]> getChangeTimes()         { return getParams().get(changeTimesName); }
}
