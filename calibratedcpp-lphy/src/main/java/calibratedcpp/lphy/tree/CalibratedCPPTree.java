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
import static calibratedcpp.lphy.tree.CPPUtils.getBD;
import static calibratedcpp.lphy.tree.CPPUtils.sampleTimes;
import static calibratedcpp.lphy.tree.CPPUtils.simRandomStem;

/**
 * Constant-rate calibrated CPP: node ages come from the birth-death law with
 * constant birthRate/deathRate (see {@link CPPUtils}). All tree-building logic
 * lives in {@link AbstractCalibratedCPPTree}; this class only resolves the
 * rate parameters and exposes them to that law.
 */
public class CalibratedCPPTree extends AbstractCalibratedCPPTree {

    Value<Number> birthRate;
    Value<Number> deathRate;
    Value<Number> diversification;
    Value<Number> turnover;

    private double resolvedBirthRate;
    private double resolvedDeathRate;

    public CalibratedCPPTree(@ParameterInfo(name = BirthDeathConstants.lambdaParamName, description = "per-lineage birth rate.", optional = true) Value<Number> birthRate,
                             @ParameterInfo(name = BirthDeathConstants.muParamName, description = "per-lineage death rate.", optional = true) Value<Number> deathRate,
                             @ParameterInfo(name = BirthDeathConstants.diversificationParamName, description = "diversification rate (lambda - mu), optional alternative to lambda/mu.", optional = true) Value<Number> diversification,
                             @ParameterInfo(name = BirthDeathConstants.turnoverParamName, description = "turnover (mu/lambda), optional alternative to lambda/mu.", optional = true) Value<Number> turnover,
                             @ParameterInfo(name = BirthDeathConstants.rhoParamName, description = "sampling probability") Value<Number> rho,
                             @ParameterInfo(name = DistributionConstants.nParamName, description = "the total number of taxa; omit for a random number of tips (uncalibrated only).", optional = true) Value<Integer> n,
                             @ParameterInfo(name = calibrationsName, description = "an array of calibrations generated from a MRCA prior (i.e. ConditionedMRCAPrior or MRCAPrior)", optional = true) Value<CalibrationArray> calibrations,
                             @ParameterInfo(name = otherTaxaNames, description = "a string array of taxa names for non-calibrated tips", optional = true) Value<String[]> otherNames,
                             @ParameterInfo(name = stemAgeName, description = "the stem age working as condition time", optional = true) Value<Number> stemAge,
                             @ParameterInfo(name = rootAgeName, description = "the root age to condition on when no calibrations are provided", optional = true) Value<Number> rootAge) {
        super(n, rho, calibrations, otherNames, stemAge, rootAge);

        int count = 0;
        if (birthRate != null) count++;
        if (deathRate != null) count++;
        if (diversification != null) count++;
        if (turnover != null) count++;

        if (count != 2) {
            throw new IllegalArgumentException(
                    "Must specify exactly two of: birthRate, deathRate, diversification, turnover."
            );
        }

        if (stemAge != null && calibrations != null && calibrations.value().getCalibrationArray()[0].getTaxa().length == n.value()) {
            LoggerUtils.log.warning("Stem age will be ignored if root calibration is provided.");
        }

        this.birthRate = birthRate;
        this.deathRate = deathRate;
        this.diversification = diversification;
        this.turnover = turnover;
    }

    @GeneratorInfo(name = "CalibratedCPP", examples = {"CalibratedCPPTree.lphy"},
            description = "The Calibrated Coalescent Point Process (calibrated CPP) method accepts one or more clade taxa and generates a tip-labelled time tree. If a root age is provided, the method is conditioned on root age. If the stem age is provided, the origin is the stem age.")
    @Override
    public RandomVariable<TimeTree> sample() {
        return super.sample();
    }

    @Override
    protected void resolveRates() {
        if (getBirthRate() != null && getDeathRate() != null) { //if both lambda and mu are given
            resolvedBirthRate = getBirthRate().value().doubleValue();
            resolvedDeathRate = getDeathRate().value().doubleValue();
        } else if (getDiversificationRate() != null && getTurnover() != null) { //if both diversification and turnover are given
            double diversificationRate = getDiversificationRate().value().doubleValue();
            double turnoverVal = getTurnover().value().doubleValue();
            double[] bd = getBD(diversificationRate, turnoverVal);
            resolvedBirthRate = bd[0];
            resolvedDeathRate = bd[1];
        } else if (getBirthRate() != null && getDiversificationRate() != null) {
            resolvedBirthRate = getBirthRate().value().doubleValue();
            resolvedDeathRate = resolvedBirthRate - getDiversificationRate().value().doubleValue();
        } else if (getBirthRate() != null && getTurnover() != null) {
            resolvedBirthRate = getBirthRate().value().doubleValue();
            resolvedDeathRate = resolvedBirthRate * getTurnover().value().doubleValue();
        } else if (getDeathRate() != null && getDiversificationRate() != null) {
            resolvedDeathRate = getDeathRate().value().doubleValue();
            resolvedBirthRate = getDiversificationRate().value().doubleValue() + resolvedDeathRate;
        } else if (getDeathRate() != null && getTurnover() != null) {
            resolvedDeathRate = getDeathRate().value().doubleValue();
            resolvedBirthRate = resolvedDeathRate / getTurnover().value().doubleValue();
        } else {
            throw new IllegalArgumentException("Invalid parameter combination, should provide either two of birth rate, death rate, diversification, and turnover.");
        }
    }

    @Override
    protected double logCDF(double t) {
        return Math.log(CDF(resolvedBirthRate, resolvedDeathRate, getSamplingProb().value().doubleValue(), t));
    }

    // Constant-rate keeps its tuned closed-form samplers (which handle the subcritical
    // saturation regime) instead of the base's generic numerical inverter.
    @Override
    protected double[] sampleAges(double lowerTime, double upperTime, int nSims) {
        return sampleTimes(resolvedBirthRate, resolvedDeathRate, getSamplingProb().value().doubleValue(), lowerTime, upperTime, nSims);
    }

    @Override
    protected double sampleStemAge(double greaterThan, int nTaxa) {
        return simRandomStem(resolvedBirthRate, resolvedDeathRate, greaterThan, nTaxa);
    }

    @Override
    protected AbstractCalibratedCPPTree newSubClade(int nTaxa, CalibrationArray subCalibrations) {
        return new CalibratedCPPTree(new Value<>("", resolvedBirthRate), new Value<>("", resolvedDeathRate), null, null,
                getSamplingProb(), new Value<>("n", nTaxa), new Value<>("", subCalibrations), null, null, null);
    }

    @Override
    public Map<String, Value> getParams() {
        Map<String, Value> map = super.getParams();
        if (birthRate != null) map.put(BirthDeathConstants.lambdaParamName, birthRate);
        if (deathRate != null) map.put(BirthDeathConstants.muParamName, deathRate);
        if (diversification != null) map.put(BirthDeathConstants.diversificationParamName, diversification);
        if (turnover != null) map.put(BirthDeathConstants.turnoverParamName, turnover);
        return map;
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(BirthDeathConstants.lambdaParamName)) birthRate = value;
        else if (paramName.equals(BirthDeathConstants.muParamName)) deathRate = value;
        else if (paramName.equals(BirthDeathConstants.diversificationParamName)) diversification = value;
        else if (paramName.equals(BirthDeathConstants.turnoverParamName)) turnover = value;
        else super.setParam(paramName, value);
    }

    public Value<Number> getBirthRate() {
        return getParams().get(BirthDeathConstants.lambdaParamName);
    }

    public Value<Number> getDeathRate() {
        return getParams().get(BirthDeathConstants.muParamName);
    }

    public Value<Number> getDiversificationRate() {
        return getParams().get(BirthDeathConstants.diversificationParamName);
    }

    public Value<Number> getTurnover() {
        return getParams().get(BirthDeathConstants.turnoverParamName);
    }
}
