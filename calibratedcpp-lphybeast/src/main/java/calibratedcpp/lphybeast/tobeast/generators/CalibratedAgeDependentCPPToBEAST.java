package calibratedcpp.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeInterface;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.distribution.Exponential;
import beast.base.spec.inference.distribution.Gamma;
import beast.base.spec.inference.distribution.LogNormal;
import beast.base.spec.inference.distribution.ScalarDistribution;
import beast.base.spec.inference.distribution.Uniform;
import beast.base.spec.inference.parameter.RealScalarParam;
import calibratedcpp.CalibratedAgeDependentExtinctionModel;
import calibratedcpp.lphy.prior.Calibration;
import calibratedcpp.lphy.prior.CalibrationArray;
import calibratedcpp.lphy.prior.ConditionedMRCAPrior;
import calibratedcpp.lphy.tree.CalibratedAgeDependentCPPTree;
import lphy.core.model.Generator;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

import java.util.ArrayList;
import java.util.List;

import static lphybeast.tobeast.TaxaUtils.getTaxonSet;

public class CalibratedAgeDependentCPPToBEAST
        implements GeneratorToBEAST<CalibratedAgeDependentCPPTree, CalibratedAgeDependentExtinctionModel> {

    @Override
    public CalibratedAgeDependentExtinctionModel generatorToBEAST(
            CalibratedAgeDependentCPPTree generator, BEASTInterface value, BEASTContext context) {

        List<TaxonSet> taxonSets = new ArrayList<>();
        CalibratedAgeDependentExtinctionModel model = new CalibratedAgeDependentExtinctionModel();
        model.setInputValue("tree", value);

        // When no calibrations are provided rootAge drives conditioning, so conditionOnRoot=true.
        boolean hasCalibrations = generator.getCalibrations() != null;
        boolean rootConditioned = !hasCalibrations || generator.getRootCondition();
        model.setInputValue("conditionOnRoot", rootConditioned);
        Boolean conditionOnCalibrationsOverride = MRCAPriorCalibrationUtils.getConditionOnCalibrationsOverride();
        model.setInputValue("conditionOnCalibrations",
                conditionOnCalibrationsOverride != null ? conditionOnCalibrationsOverride : hasCalibrations);

        if (!rootConditioned) {
            model.setInputValue("origin", new RealScalarParam<>(
                    generator.getOrigin().value(), PositiveReal.INSTANCE));
        }

        model.setInputValue("birthRate", context.getAsRealScalar(generator.getBirthRate()));
        model.setInputValue("rho",       context.getAsRealScalar(generator.getSamplingProb()));

        model.setInputValue("lifetimeDistribution",
                buildLifetimeDistribution(generator.getLifetime(), context));

        // Remove the 'lifetime' state node and its auto-generated prior; the distribution
        // parameters are wired directly into lifetimeDistribution instead.
        Value<Number> lifetimeValue = generator.getLifetime();
        BEASTInterface lifetimeBEAST = context.getBEASTObject(lifetimeValue);
        if (lifetimeBEAST != null) context.removeBEASTObject(lifetimeBEAST);
        if (lifetimeValue.getGenerator() != null) {
            BEASTInterface lifetimePrior = context.getBEASTObject(lifetimeValue.getGenerator());
            if (lifetimePrior != null) context.removeBEASTObject(lifetimePrior);
        }

        if (!hasCalibrations) {
            // rootAge-only mode: no calibration clades or CalibrationPrior needed.
            model.setInputValue("calibrations", new ArrayList<>());
            model.initAndValidate();
            return model;
        }

        // calibration clades
        Calibration[] calibrationsFromGenerator = generator.getCalibrations().value().getCalibrationArray();
        for (Calibration calibration : calibrationsFromGenerator) {
            taxonSets.add(getTaxonSet((TreeInterface) value, calibration.getTaxa()));
        }
        model.setInputValue("calibrations", taxonSets);
        model.initAndValidate();

        Value<CalibrationArray> calibrationsValue = generator.getCalibrations();
        if (MRCAPriorCalibrationUtils.isIndependentMRCAPriorSource(calibrationsValue)) {
            // calibrations came from toArray/joinArray over independent UniformMRCA(taxa=,upper=,lower=)
            // calls, not ConditionedMRCAPrior's joint density — build one plain, independently-
            // bounded MRCAPrior per clade instead of a CalibrationPrior.
            MRCAPriorCalibrationUtils.buildIndependentMRCAPriors(
                    MRCAPriorCalibrationUtils.collectIndependentCalibrationGenerators(calibrationsValue),
                    taxonSets, value, context);
            return model;
        }

        // Default: preserve ConditionedMRCAPrior's joint/nested density as a single BEAST
        // CalibrationPrior. -DcalibratedcppMRCAPrior=true reverts to the older behaviour: each
        // clade gets its own independent monophyly + Uniform(lower,upper) MRCAPrior over its
        // original bounds, same as the UniformMRCA path above (and what this converter always
        // did before this option existed).
        ConditionedMRCAPrior conditionedMRCAPrior =
                (ConditionedMRCAPrior) calibrationsValue.getInputs().get(0);
        Calibration[] calibrationSpecs = conditionedMRCAPrior.getCalibrations().value();
        if (MRCAPriorCalibrationUtils.isMrcaPriorMode()) {
            for (int i = 0; i < calibrationSpecs.length; i++) {
                beast.base.spec.evolution.tree.MRCAPrior mrcaPrior = MRCAPriorCalibrationUtils.buildBoundedMRCAPrior(
                        value, taxonSets.get(i), calibrationSpecs[i].getLower(), calibrationSpecs[i].getUpper());
                context.addBEASTObject(mrcaPrior, conditionedMRCAPrior);
                context.addExtraLoggable(mrcaPrior);
            }
        } else {
            double coverage = conditionedMRCAPrior.getCoverage() != null
                    ? conditionedMRCAPrior.getCoverage().value().doubleValue() : 0.9;
            MRCAPriorCalibrationUtils.buildCalibrationPrior(
                    value, taxonSets, calibrationSpecs, coverage, context, conditionedMRCAPrior);
        }

        return model;
    }

    /**
     * Inspects the generator of the LPhy {@code lifetime} value and constructs the
     * corresponding BEAST {@link ScalarDistribution} with its parameters wired from
     * the BEAST context. Named LPhy parameters (e.g. {@code shape.lifetime ~ LogNormal(...)})
     * become live BEAST state nodes; anonymous constants are inlined.
     *
     * <p>Supported LPhy distributions:
     * <ul>
     *   <li>{@code Gamma(shape, scale)}      → BEAST {@code Gamma(alpha=shape, theta=scale)}</li>
     *   <li>{@code LogNormal(meanlog, sdlog)} → BEAST {@code LogNormal(M=meanlog, S=sdlog)}</li>
     *   <li>{@code Exp(mean)}                → BEAST {@code Exponential(mean=mean)}</li>
     *   <li>constant {@code v}               → BEAST {@code Uniform(lower=v−0.01, upper=v+0.01)}</li>
     * </ul>
     */
    private ScalarDistribution buildLifetimeDistribution(Value<Number> lifetime, BEASTContext context) {
        Generator<?> gen = lifetime.getGenerator();
        if (gen == null) {
            double v = lifetime.value().doubleValue();
            return buildUniform(v - 0.01, v + 0.01);
        }

        if (gen instanceof lphy.base.distribution.Gamma g) {
            return buildGamma(g, context);
        } else if (gen instanceof lphy.base.distribution.LogNormal ln) {
            return buildLogNormal(ln, context);
        } else if (gen instanceof lphy.base.distribution.Exp e) {
            return buildExponential(e, context);
        } else {
            throw new IllegalArgumentException(
                    "Unsupported lifetime distribution: " + gen.getClass().getSimpleName()
                    + ". Supported distributions: Gamma, LogNormal, Exp.");
        }
    }

    private Uniform buildUniform(double lower, double upper) {
        Uniform u = new Uniform();
        u.setInputValue("lower", new RealScalarParam<>(lower, Real.INSTANCE));
        u.setInputValue("upper", new RealScalarParam<>(upper, Real.INSTANCE));
        u.initAndValidate();
        return u;
    }

    private Gamma buildGamma(lphy.base.distribution.Gamma g, BEASTContext context) {
        Gamma beastGamma = new Gamma();
        beastGamma.setInputValue("alpha", context.getAsRealScalar(g.getShape()));
        beastGamma.setInputValue("theta", context.getAsRealScalar(g.getScale()));
        beastGamma.initAndValidate();
        return beastGamma;
    }

    private LogNormal buildLogNormal(lphy.base.distribution.LogNormal ln, BEASTContext context) {
        LogNormal beastLN = new LogNormal();
        beastLN.setInputValue("M", context.getAsRealScalar(ln.getMeanLog()));
        beastLN.setInputValue("S", context.getAsRealScalar(ln.getSDLog()));
        beastLN.initAndValidate();
        return beastLN;
    }

    private Exponential buildExponential(lphy.base.distribution.Exp e, BEASTContext context) {
        Exponential beastExp = new Exponential();
        beastExp.setInputValue("mean", context.getAsRealScalar(
                (Value<?>) e.getParams().get("mean")));
        beastExp.initAndValidate();
        return beastExp;
    }

    @Override
    public Class<CalibratedAgeDependentCPPTree> getGeneratorClass() {
        return CalibratedAgeDependentCPPTree.class;
    }

    @Override
    public Class<CalibratedAgeDependentExtinctionModel> getBEASTClass() {
        return CalibratedAgeDependentExtinctionModel.class;
    }
}
