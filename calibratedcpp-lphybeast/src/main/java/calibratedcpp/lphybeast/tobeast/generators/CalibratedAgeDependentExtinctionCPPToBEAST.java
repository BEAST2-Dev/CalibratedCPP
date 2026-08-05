package calibratedcpp.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeInterface;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.distribution.Exponential;
import beast.base.spec.inference.distribution.Gamma;
import beast.base.spec.inference.distribution.ScalarDistribution;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.SimplexParam;
import calibratedcpp.CalibratedAgeDependentExtinctionModel;
import calibratedcpp.distribution.Weibull;
import calibratedcpp.distribution.WeibullMixture;
import calibratedcpp.lphy.prior.Calibration;
import calibratedcpp.lphy.prior.CalibrationArray;
import calibratedcpp.lphy.prior.ConditionedMRCAPrior;
import calibratedcpp.lphy.tree.CalibratedAgeDependentExtinctionTree;
import calibratedcpp.lphy.tree.ExpLifetime;
import calibratedcpp.lphy.tree.GammaLifetime;
import calibratedcpp.lphy.tree.WeibullLifetime;
import calibratedcpp.lphy.tree.WeibullMixtureLifetime;
import lphy.core.model.Generator;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

import java.util.ArrayList;
import java.util.List;

import static lphybeast.tobeast.TaxaUtils.getTaxonSet;

/**
 * Maps the LPhy age-dependent-extinction generator to the BEAST {@link CalibratedAgeDependentExtinctionModel}.
 * The engine-neutral {@code LifetimeDistribution} value is mapped back to a BEAST {@link ScalarDistribution}
 * by inspecting the lifetime function that produced it (weibullLifetime → BEAST Weibull, etc.). The
 * origin / conditioning / calibration wiring matches {@link CalibratedCPPToBEAST}.
 */
public class CalibratedAgeDependentExtinctionCPPToBEAST
        implements GeneratorToBEAST<CalibratedAgeDependentExtinctionTree, CalibratedAgeDependentExtinctionModel> {

    @Override
    public CalibratedAgeDependentExtinctionModel generatorToBEAST(
            CalibratedAgeDependentExtinctionTree generator, BEASTInterface value, BEASTContext context) {

        CalibratedAgeDependentExtinctionModel model = new CalibratedAgeDependentExtinctionModel();
        model.setInputValue("tree", value);

        boolean hasCalibrations = generator.getCalibrations() != null;
        boolean rootConditioned = !hasCalibrations || generator.getRootCondition();
        model.setInputValue("conditionOnRoot", rootConditioned);
        model.setInputValue("conditionOnCalibrations", hasCalibrations);

        if (!rootConditioned) {
            model.setInputValue("origin", new RealScalarParam<>(generator.getOrigin().value(), PositiveReal.INSTANCE));
        }

        if (generator.getBirthRate() != null) {
            model.setInputValue("birthRate", context.getAsRealScalar(generator.getBirthRate()));
        } else {
            // reproductiveNumber given: the BEAST model takes birthRate, so pass the derived value
            // lambda = R0 / mean(lifetime). Note this fixes birthRate at that value in the XML rather
            // than inferring R0 directly (a derived-parameter conversion would be needed for that).
            double lambda = generator.getReproductiveNumber().value().doubleValue()
                    / generator.getLifetime().value().mean();
            model.setInputValue("birthRate", new RealScalarParam<>(lambda, PositiveReal.INSTANCE));
        }
        model.setInputValue("rho", context.getAsRealScalar(generator.getSamplingProb()));
        model.setInputValue("lifetimeDistribution", buildLifetimeDistribution(generator.getLifetime(), context));

        if (!hasCalibrations) {
            model.setInputValue("calibrations", new ArrayList<>());
            model.initAndValidate();
            return model;
        }

        List<TaxonSet> taxonSets = new ArrayList<>();
        for (Calibration calibration : generator.getCalibrations().value().getCalibrationArray()) {
            taxonSets.add(getTaxonSet((TreeInterface) value, calibration.getTaxa()));
        }
        model.setInputValue("calibrations", taxonSets);
        model.initAndValidate();

        Value<CalibrationArray> calibrationsValue = generator.getCalibrations();
        if (MRCAPriorCalibrationUtils.isIndependentMRCAPriorSource(calibrationsValue)) {
            MRCAPriorCalibrationUtils.buildIndependentMRCAPriors(
                    MRCAPriorCalibrationUtils.collectIndependentCalibrationGenerators(calibrationsValue),
                    taxonSets, value, context);
            return model;
        }

        ConditionedMRCAPrior conditionedMRCAPrior = (ConditionedMRCAPrior) calibrationsValue.getInputs().get(0);
        Calibration[] calibrationSpecs = conditionedMRCAPrior.getCalibrations().value();
        for (int i = 0; i < calibrationSpecs.length; i++) {
            beast.base.spec.evolution.tree.MRCAPrior mrcaPrior = MRCAPriorCalibrationUtils.buildBoundedMRCAPrior(
                    value, taxonSets.get(i), calibrationSpecs[i].getLower(), calibrationSpecs[i].getUpper());
            context.addBEASTObject(mrcaPrior, conditionedMRCAPrior);
            context.addExtraLoggable(mrcaPrior);
        }

        return model;
    }

    /**
     * Recover a BEAST {@link ScalarDistribution} from the LPhy {@code LifetimeDistribution} value by
     * inspecting the function that produced it. Its named parameters (e.g. {@code shape.lifetime ~ ...})
     * become live BEAST state nodes; constants are inlined — both handled by {@code getAsRealScalar}.
     */
    private ScalarDistribution buildLifetimeDistribution(Value<?> lifetimeDist, BEASTContext context) {
        Generator<?> gen = lifetimeDist.getGenerator();
        if (gen == null) {
            throw new IllegalArgumentException("lifetimeDist must come from a lifetime function "
                    + "(weibullLifetime, gammaLifetime or expLifetime).");
        }
        if (gen instanceof WeibullLifetime) {
            Weibull w = new Weibull();
            w.setInputValue("shape", context.getAsRealScalar(gen.getParams().get(WeibullLifetime.shapeParamName)));
            w.setInputValue("scale", context.getAsRealScalar(gen.getParams().get(WeibullLifetime.scaleParamName)));
            w.initAndValidate();
            return w;
        } else if (gen instanceof GammaLifetime) {
            Gamma g = new Gamma();
            g.setInputValue("alpha", context.getAsRealScalar(gen.getParams().get(GammaLifetime.shapeParamName)));
            g.setInputValue("theta", context.getAsRealScalar(gen.getParams().get(GammaLifetime.scaleParamName)));
            g.initAndValidate();
            return g;
        } else if (gen instanceof ExpLifetime) {
            Exponential e = new Exponential();
            e.setInputValue("mean", context.getAsRealScalar(gen.getParams().get(ExpLifetime.meanParamName)));
            e.initAndValidate();
            return e;
        } else if (gen instanceof WeibullMixtureLifetime) {
            WeibullMixture wm = new WeibullMixture();
            wm.setInputValue("mean",   context.getAsRealScalar(gen.getParams().get(WeibullMixtureLifetime.meanParamName)));
            wm.setInputValue("shape1", context.getAsRealScalar(gen.getParams().get(WeibullMixtureLifetime.shape1ParamName)));
            wm.setInputValue("shape2", context.getAsRealScalar(gen.getParams().get(WeibullMixtureLifetime.shape2ParamName)));
            Value<?> weights = gen.getParams().get(WeibullMixtureLifetime.weightsParamName);
            if (weights != null) {
                // Constant weights become a fixed Simplex; equal weights (0.5, 0.5) are BEAST's default.
                Double[] w = (Double[]) weights.value();
                double sum = 0.0; for (Double wi : w) sum += wi;
                double[] wd = new double[w.length];
                for (int i = 0; i < w.length; i++) wd[i] = w[i] / sum;
                wm.setInputValue("weights", new SimplexParam(wd));
            }
            wm.initAndValidate();
            return wm;
        }
        throw new IllegalArgumentException("Unsupported lifetime function: " + gen.getClass().getSimpleName()
                + ". Supported: weibullLifetime, gammaLifetime, expLifetime, weibullMixtureLifetime.");
    }

    @Override
    public Class<CalibratedAgeDependentExtinctionTree> getGeneratorClass() {
        return CalibratedAgeDependentExtinctionTree.class;
    }

    @Override
    public Class<CalibratedAgeDependentExtinctionModel> getBEASTClass() {
        return CalibratedAgeDependentExtinctionModel.class;
    }
}
