package calibratedcpp.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.Real;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.evolution.tree.MRCAPrior;
import beast.base.spec.inference.distribution.Uniform;
import beast.base.spec.inference.parameter.RealScalarParam;
import calibratedcpp.lphy.prior.Calibration;
import calibratedcpp.lphy.prior.CalibrationArray;
import calibratedcpp.lphy.prior.OffsetExponentialMRCA;
import calibratedcpp.lphy.prior.UniformMRCA;
import calibratedcpp.lphy.prior.toCalibrationArray;
import calibrationprior.CalibrationCladePrior;
import calibrationprior.CalibrationPrior;
import lphy.core.model.Generator;
import lphy.core.model.GraphicalModelNode;
import lphy.core.model.Value;
import lphy.core.vectorization.array.ArrayFunction;
import lphybeast.BEASTContext;

import java.util.ArrayList;
import java.util.List;

/**
 * Support for calibrations built from independent per-clade {@link UniformMRCA} /
 * {@link OffsetExponentialMRCA} generators combined via {@code toArray(calibrations=[...])},
 * as opposed to the joint density produced by {@link calibratedcpp.lphy.prior.ConditionedMRCAPrior}.
 * The array may mix both generator types and each is dispatched to its matching converter.
 */
public class MRCAPriorCalibrationUtils {

    /** True if the calibrations came from {@code toArray(...)} rather than {@code ConditionedMRCAPrior}. */
    public static boolean isIndependentMRCAPriorSource(Value<CalibrationArray> calibrationsValue) {
        return calibrationsValue.getGenerator() instanceof toCalibrationArray;
    }

    /**
     * Returns the per-clade calibration generators behind a {@code toArray(calibrations=[...])}
     * literal, index-aligned with the resulting {@link CalibrationArray} and any TaxonSet list
     * built from it.
     */
    public static List<Generator<?>> collectIndependentCalibrationGenerators(Value<?> calibrationsValue) {
        List<Generator<?>> out = new ArrayList<>();
        Generator<?> gen = calibrationsValue.getGenerator();
        if (!(gen instanceof toCalibrationArray tca)) {
            throw new IllegalArgumentException(
                    "Expected calibrations to come from toArray(calibrations=[...]), got: "
                            + (gen == null ? "constant value" : gen.getClass().getSimpleName()));
        }
        Value<?> arrayValue = tca.getParams().get(toCalibrationArray.calibrationsParamName);
        Generator<?> arrayGen = arrayValue.getGenerator();
        if (!(arrayGen instanceof ArrayFunction<?> arrayFunction)) {
            throw new IllegalArgumentException(
                    "Expected an array literal of calibrations (e.g. [cal1, cal2, ...]), got: "
                            + (arrayGen == null ? "constant value" : arrayGen.getClass().getSimpleName()));
        }
        for (Value<?> element : arrayFunction.getValues()) {
            Generator<?> elementGen = element.getGenerator();
            if (!(elementGen instanceof UniformMRCA) && !(elementGen instanceof OffsetExponentialMRCA)) {
                throw new IllegalArgumentException(
                        "Expected each calibration in the array to come from UniformMRCA or OffsetExponentialMRCA, got: "
                                + (elementGen == null ? "constant value" : elementGen.getClass().getSimpleName()));
            }
            out.add(elementGen);
        }
        return out;
    }

    /**
     * Builds one {@code MRCAPrior} per calibration generator, reuse the correspond {@link TaxonSet}s
     */
    public static void buildIndependentMRCAPriors(
            List<Generator<?>> calibrationGenerators, List<TaxonSet> taxonSets,
            BEASTInterface treeValue, BEASTContext context) {

        UniformMRCAToBEAST uniformConverter = new UniformMRCAToBEAST();
        OffsetExponentialMRCAToBEAST offsetExponentialConverter = new OffsetExponentialMRCAToBEAST();
        for (int i = 0; i < calibrationGenerators.size(); i++) {
            Generator<?> gen = calibrationGenerators.get(i);
            if (gen instanceof UniformMRCA uniformMRCA) {
                uniformConverter.generatorToBEAST(uniformMRCA, treeValue, taxonSets.get(i), context);
            } else if (gen instanceof OffsetExponentialMRCA offsetExponentialMRCA) {
                offsetExponentialConverter.generatorToBEAST(offsetExponentialMRCA, treeValue, taxonSets.get(i), context);
            } else {
                throw new IllegalArgumentException("Unsupported calibration generator: " + gen.getClass().getSimpleName());
            }
        }
    }

    /**
     * Builds a plain BEAST {@code MRCAPrior(monophyletic=true, distr=Uniform(lower,upper))} from
     * raw double bounds, for calibration sources (e.g. {@code ConditionedMRCAPrior}) that only
     * expose bounds as plain numbers rather than as LPhy {@code Value}s.
     */
    public static MRCAPrior buildBoundedMRCAPrior(BEASTInterface treeValue, TaxonSet taxonSet, double lower, double upper) {
        Uniform uniform = new Uniform();
        uniform.setInputValue("lower", new RealScalarParam<>(lower, Real.INSTANCE));
        uniform.setInputValue("upper", new RealScalarParam<>(upper, Real.INSTANCE));
        uniform.initAndValidate();

        MRCAPrior mrcaPrior = new MRCAPrior();
        mrcaPrior.setInputValue("tree", treeValue);
        mrcaPrior.setInputValue("taxonset", taxonSet);
        mrcaPrior.setInputValue("monophyletic", true);
        mrcaPrior.setInputValue("distr", uniform);
        mrcaPrior.initAndValidate();
        return mrcaPrior;
    }

    // Converter switches, set from the -MRCAPrior / -conditionOnCalibrations flags parsed by
    // CalibratedCPPLPhyBeastMain. Read here so all converters share one parsing path.

    /**
     * True if {@code ConditionedMRCAPrior} calibrations should become independent per-clade
     * {@code MRCAPrior(Uniform)}s instead of one joint {@code CalibrationPrior}.
     */
    public static boolean isMrcaPriorMode() {
        return Boolean.parseBoolean(System.getProperty("calibratedcppMRCAPrior", "false"));
    }

    /** The {@code conditionOnCalibrations} override, or null if unset so the caller's default applies. */
    public static Boolean getConditionOnCalibrationsOverride() {
        String raw = System.getProperty("calibratedcppConditionOnCalibrations");
        if (raw == null || raw.isBlank()) return null;
        return Boolean.parseBoolean(raw);
    }

    /**
     * Builds one joint {@code CalibrationPrior} over all clades, one {@code CalibrationCladePrior}
     * per (spec, taxon set) pair, preserving the nested/overlap structure {@code ConditionedMRCAPrior}
     * encodes rather than reducing it to independent per-clade bounds.
     */
    public static void buildCalibrationPrior(
            BEASTInterface treeValue, List<TaxonSet> taxonSets, Calibration[] calibrationSpecs,
            double confidenceLevel, BEASTContext context, GraphicalModelNode lphyRef) {

        List<CalibrationCladePrior> cladePriors = new ArrayList<>();
        for (int i = 0; i < calibrationSpecs.length; i++) {
            CalibrationCladePrior ccp = new CalibrationCladePrior();
            ccp.setInputValue("taxa", taxonSets.get(i));
            ccp.setInputValue("upperAge", new RealScalarParam<>(calibrationSpecs[i].getUpper(), NonNegativeReal.INSTANCE));
            ccp.setInputValue("lowerAge", new RealScalarParam<>(calibrationSpecs[i].getLower(), NonNegativeReal.INSTANCE));
            ccp.setInputValue("confidenceLevel", new RealScalarParam<>(confidenceLevel, UnitInterval.INSTANCE));
            ccp.initAndValidate();
            cladePriors.add(ccp);
        }

        CalibrationPrior calibrationPrior = new CalibrationPrior();
        calibrationPrior.setInputValue("tree", treeValue);
        calibrationPrior.setInputValue("calibration", cladePriors);
        calibrationPrior.initAndValidate();

        context.addBEASTObject(calibrationPrior, lphyRef);
        context.addExtraLoggable(calibrationPrior);
    }
}
