package calibratedcpp.lphy.tree;

import calibratedcpp.lphy.prior.CalibrationArray;
import lphy.base.distribution.DistributionConstants;
import lphy.base.evolution.birthdeath.BirthDeathConstants;
import lphy.base.evolution.tree.TimeTree;
import lphy.core.logger.LoggerUtils;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

import java.util.Map;

import static calibratedcpp.lphy.tree.CPPUtils.CDF;
import static calibratedcpp.lphy.tree.CPPUtils.sampleTimes;
import static calibratedcpp.lphy.tree.CPPUtils.simRandomStem;

/**
 * Age-dependent calibrated CPP: individuals have a lifetime distribution rather
 * than a constant death rate. The sampled lifetime value is used as a point
 * estimate 1/lifetime for an effective constant death rate, so this class feeds
 * the exact same {@link CPPUtils} birth-death law into {@link AbstractCalibratedCPPTree}
 * as {@link CalibratedCPPTree} does; only the resolved rate differs. The full
 * age-dependent likelihood is evaluated in BEAST via CalibratedAgeDependentExtinctionModel.
 */
public class CalibratedAgeDependentCPPTree extends AbstractCalibratedCPPTree {

    Value<Number> birthRate;
    Value<Number> lifetime;

    private double resolvedBirthRate;
    private double resolvedEffectiveDeathRate;

    public static final String lifetimeName = "lifetime";

    public CalibratedAgeDependentCPPTree(
            @ParameterInfo(name = BirthDeathConstants.lambdaParamName,
                    description = "per-lineage birth rate.") Value<Number> birthRate,
            @ParameterInfo(name = BirthDeathConstants.rhoParamName,
                    description = "sampling probability.") Value<Number> rho,
            @ParameterInfo(name = DistributionConstants.nParamName,
                    description = "the total number of taxa.") Value<Integer> n,
            @ParameterInfo(name = lifetimeName,
                    description = "individual lifetime distribution; must come from a continuous distribution " +
                            "(e.g. Gamma, LogNormal, Exp). The sampled value is used as a point estimate of " +
                            "the mean lifetime for prior-predictive tree simulation. The full distribution is " +
                            "passed to BEAST for exact likelihood evaluation.") Value<Number> lifetime,
            @ParameterInfo(name = calibrationsName,
                    description = "an array of calibrations generated from a MRCA prior " +
                            "(i.e. ConditionedMRCAPrior or MRCAPrior).",
                    optional = true) Value<CalibrationArray> calibrations,
            @ParameterInfo(name = otherTaxaNames,
                    description = "a string array of taxa names for non-calibrated tips.",
                    optional = true) Value<String[]> otherNames,
            @ParameterInfo(name = stemAgeName,
                    description = "the stem age working as condition time.",
                    optional = true) Value<Number> stemAge,
            @ParameterInfo(name = rootAgeName,
                    description = "the root age to condition on when no calibrations are provided.",
                    optional = true) Value<Number> rootAge) {
        super(n, rho, calibrations, otherNames, stemAge, rootAge);

        if (lifetime == null) {
            throw new NullPointerException("lifetime should not be null!");
        }

        if (stemAge != null && calibrations != null && calibrations.value().getCalibrationArray()[0].getTaxa().length == n.value()) {
            LoggerUtils.log.warning("Stem age will be ignored if root calibration is provided.");
        }

        this.birthRate = birthRate;
        this.lifetime = lifetime;
    }

    @GeneratorInfo(name = "CalibratedAgeDependentCPP", examples = {},
            description = "The Calibrated Coalescent Point Process with age-dependent individual lifetimes. "
                    + "The lifetime parameter must come from a continuous distribution (Gamma, LogNormal, or Exp). "
                    + "Tree sampling uses effectiveDeathRate = 1/lifetime for initialisation; "
                    + "the exact age-dependent likelihood is evaluated in BEAST via CalibratedAgeDependentExtinctionModel.")
    @Override
    public RandomVariable<TimeTree> sample() {
        return super.sample();
    }

    @Override
    protected void resolveRates() {
        resolvedBirthRate = getBirthRate().value().doubleValue();
        // 1/mean_lifetime as an effective death rate for BD-based tree initialisation
        resolvedEffectiveDeathRate = 1.0 / getLifetime().value().doubleValue();
    }

    @Override
    protected double logCDF(double t) {
        return Math.log(CDF(resolvedBirthRate, resolvedEffectiveDeathRate, getSamplingProb().value().doubleValue(), t));
    }

    // Point-estimate constant-rate law keeps the closed-form samplers.
    @Override
    protected double[] sampleAges(double lowerTime, double upperTime, int nSims) {
        return sampleTimes(resolvedBirthRate, resolvedEffectiveDeathRate, getSamplingProb().value().doubleValue(), lowerTime, upperTime, nSims);
    }

    @Override
    protected double sampleStemAge(double greaterThan, int nTaxa) {
        return simRandomStem(resolvedBirthRate, resolvedEffectiveDeathRate, greaterThan, nTaxa);
    }

    @Override
    protected AbstractCalibratedCPPTree newSubClade(int nTaxa, CalibrationArray subCalibrations) {
        // recursive: subclade also uses the same lifetime distribution
        return new CalibratedAgeDependentCPPTree(new Value<>("", resolvedBirthRate), getSamplingProb(),
                new Value<>("n", nTaxa), getLifetime(), new Value<>("", subCalibrations), null, null, null);
    }

    @Override
    public Map<String, Value> getParams() {
        Map<String, Value> map = super.getParams();
        map.put(BirthDeathConstants.lambdaParamName, birthRate);
        map.put(lifetimeName, lifetime);
        return map;
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(BirthDeathConstants.lambdaParamName)) birthRate = value;
        else if (paramName.equals(lifetimeName)) lifetime = value;
        else super.setParam(paramName, value);
    }

    public Value<Number> getBirthRate() {
        return getParams().get(BirthDeathConstants.lambdaParamName);
    }

    public Value<Number> getLifetime() {
        return getParams().get(lifetimeName);
    }
}
