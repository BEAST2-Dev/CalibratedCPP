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
 * {@link OffsetExponentialMRCA} LPhy generators (combined via {@code toArray(calibrations=[cal1,
 * cal2, ...])}, the same array-literal syntax used by {@code ConditionedMRCAPrior(calibrations=
 * [cal1, cal2, ...])}), as opposed to the joint calibration density produced by
 * {@link calibratedcpp.lphy.prior.ConditionedMRCAPrior}. The array can freely mix both generator
 * types -- each element owns its own age-distribution shape (Uniform vs. offset-exponential) and
 * is dispatched to the matching converter independently.
 */
public class MRCAPriorCalibrationUtils {

    /**
     * @return true if {@code calibrationsValue}'s generator is not {@code ConditionedMRCAPrior}
     *         (i.e. it was built via {@code toArray(calibrations=[...])} over independent
     *         {@code UniformMRCA}/{@code OffsetExponentialMRCA} calls instead).
     */
    public static boolean isIndependentMRCAPriorSource(Value<CalibrationArray> calibrationsValue) {
        return calibrationsValue.getGenerator() instanceof toCalibrationArray;
    }

    /**
     * Walks the {@code toArray(calibrations=[...])} array literal feeding {@code calibrationsValue}
     * and returns the individual per-clade calibration generators (each a {@code UniformMRCA} or
     * {@code OffsetExponentialMRCA}) in the same order their calibrations appear in the resulting
     * {@link CalibrationArray} — so index i here lines up with index i of
     * {@code calibrationsValue.value().getCalibrationArray()} and with any TaxonSet list built
     * from that same array.
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
     * Builds one plain BEAST {@code MRCAPrior} per calibration generator (dispatched to
     * {@code UniformMRCAToBEAST} or {@code OffsetExponentialMRCAToBEAST} by its actual type),
     * reusing the caller's already-built {@link TaxonSet} for each clade (index-aligned with
     * {@code calibrationGenerators}) so the {@code taxonset} reference matches the one used in
     * the tree model's {@code calibrations} list, rather than building a second, duplicate
     * TaxonSet for the same clade.
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

    // ------------------------------------------------------------------------------------------
    // calibratedcpp-only converter switches. Not real lphybeast CLI flags -- forwarded as -D
    // system properties by calibratedcpp-lphybeast-launcher/pom.xml's exec-maven-plugin config
    // (e.g. "mvn ... -DcalibratedcppMRCAPrior=true -DcalibratedcppConditionOnCalibrations=false"),
    // read here so both CalibratedCPPToBEAST and CalibratedAgeDependentCPPToBEAST share one
    // parsing path instead of duplicating System.getProperty calls.

    /**
     * @return true if calibrations sourced from {@code ConditionedMRCAPrior} should be converted
     *         to independent per-clade {@code MRCAPrior(Uniform)} objects (today's/the
     *         age-dependent converter's long-standing behaviour), instead of the default: a
     *         single joint {@code CalibrationPrior} preserving the nested/overlap structure.
     */
    public static boolean isMrcaPriorMode() {
        return Boolean.parseBoolean(System.getProperty("calibratedcppMRCAPrior", "false"));
    }

    /**
     * @return the {@code conditionOnCalibrations} override for
     *         {@code CalibratedBirthDeathSkylineModel}/{@code CalibratedAgeDependentBirthDeathModel},
     *         or {@code null} if not overridden (caller should fall back to its own default).
     */
    public static Boolean getConditionOnCalibrationsOverride() {
        String raw = System.getProperty("calibratedcppConditionOnCalibrations");
        if (raw == null || raw.isBlank()) return null;
        return Boolean.parseBoolean(raw);
    }

    /**
     * Builds a single joint {@code CalibrationPrior} over all given clades — one
     * {@code CalibrationCladePrior} per (calibration spec, taxon set) pair — preserving the
     * nested/overlap structure {@code ConditionedMRCAPrior} encodes (LogNormal at each disjoint
     * root, Beta on overlapping child/parent ratios, truncated LogNormal on nested
     * non-overlapping children), instead of throwing that structure away in favour of independent
     * per-clade bounds. Registered into {@code context} the same way as the independent-MRCAPrior
     * path: {@code addBEASTObject} (for provenance/XML wiring) + {@code addExtraLoggable}.
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
