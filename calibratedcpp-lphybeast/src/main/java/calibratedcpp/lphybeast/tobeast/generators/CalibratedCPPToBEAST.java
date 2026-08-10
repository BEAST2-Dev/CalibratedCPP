package calibratedcpp.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeInterface;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import calibratedcpp.CalibratedBirthDeathSkylineModel;
import calibratedcpp.SkylineParameter;
import calibratedcpp.lphy.prior.Calibration;
import calibratedcpp.lphy.prior.CalibrationArray;
import calibratedcpp.lphy.prior.ConditionedMRCAPrior;
import calibratedcpp.lphy.tree.CalibratedCPPTree;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

import java.util.ArrayList;
import java.util.List;

import static lphybeast.tobeast.TaxaUtils.getTaxonSet;

public class CalibratedCPPToBEAST implements GeneratorToBEAST<CalibratedCPPTree, CalibratedBirthDeathSkylineModel> {
    @Override
    public CalibratedBirthDeathSkylineModel generatorToBEAST(CalibratedCPPTree generator, BEASTInterface value, BEASTContext context) {
        List<TaxonSet> taxonSets = new ArrayList<>();
        CalibratedBirthDeathSkylineModel model = new CalibratedBirthDeathSkylineModel();
        model.setInputValue("tree", value);

        // without calibrations, rootAge drives conditioning
        boolean hasCalibrations = generator.getCalibrations() != null;
        boolean rootConditioned = !hasCalibrations || generator.getRootCondition();
        model.setInputValue("conditionOnRoot", rootConditioned);
        Boolean conditionOnCalibrationsOverride = MRCAPriorCalibrationUtils.getConditionOnCalibrationsOverride();
        model.setInputValue("conditionOnCalibrations",
                conditionOnCalibrationsOverride != null ? conditionOnCalibrationsOverride : hasCalibrations);

        if (!rootConditioned) {
            model.setInputValue("origin", new RealScalarParam<>(generator.getOrigin().value(), PositiveReal.INSTANCE));
        }

        if (generator.getBirthRate() != null) {
            SkylineParameter b = new SkylineParameter();
            b.setInputValue("values", context.getAsRealTensor(generator.getBirthRate()));
            b.initAndValidate();
            model.setInputValue("birthRate", b);
        }

        if (generator.getDeathRate() != null) {
            SkylineParameter d = new SkylineParameter();
            d.setInputValue("values", context.getAsRealTensor(generator.getDeathRate()));
            d.initAndValidate();
            model.setInputValue("deathRate", d);
        }

        if (generator.getTurnover() != null) {
            SkylineParameter t = new SkylineParameter();
            t.setInputValue("values", context.getAsRealTensor(generator.getTurnover()));
            t.initAndValidate();
            model.setInputValue("turnover", t);
        }

        if (generator.getDiversificationRate() != null) {
            SkylineParameter d = new SkylineParameter();
            d.setInputValue("values", context.getAsRealTensor(generator.getDiversificationRate()));
            d.initAndValidate();
            model.setInputValue("diversificationRate", d);
        }

        model.setInputValue("rho", context.getAsRealScalar(generator.getSamplingProb()));

        if (!hasCalibrations) {
            // rootAge-only mode: no calibration clades or CalibrationPrior needed.
            model.setInputValue("calibrations", new ArrayList<>());
            model.initAndValidate();
            return model;
        }

        // get clade calibrations
        Calibration[] calibrationsFromGenerator = generator.getCalibrations().value().getCalibrationArray();

        for (Calibration calibration : calibrationsFromGenerator) {
            TaxonSet taxonSet = getTaxonSet((TreeInterface) value, calibration.getTaxa());
            taxonSets.add(taxonSet);
        }

        model.setInputValue("calibrations", taxonSets);
        model.initAndValidate();

        Value<CalibrationArray> calibrationsValue = generator.getCalibrations();
        if (MRCAPriorCalibrationUtils.isIndependentMRCAPriorSource(calibrationsValue)) {
            // calibrations came from array of independent UniformMRCA(taxa=,upper=,lower=)
            // build one plain, independently-bounded MRCAPrior per clade instead of a CalibrationPrior.
            MRCAPriorCalibrationUtils.buildIndependentMRCAPriors(
                    MRCAPriorCalibrationUtils.collectIndependentCalibrationGenerators(calibrationsValue),
                    taxonSets, value, context);
            return model;
        }

        // default keeps the joint density as one CalibrationPrior
        // -MRCAPrior splits it per clade and build Uniform distributed node age for each clade
        ConditionedMRCAPrior conditionedMRCAPrior = (ConditionedMRCAPrior) calibrationsValue.getInputs().get(0);
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

    @Override
    public Class<CalibratedCPPTree> getGeneratorClass() {
        return CalibratedCPPTree.class;
    }

    @Override
    public Class<CalibratedBirthDeathSkylineModel> getBEASTClass() {
        return CalibratedBirthDeathSkylineModel.class;
    }
}
