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
import calibratedcpp.lphy.tree.CalibratedBirthDeathSkylineTree;
import lphy.core.model.Value;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

import java.util.ArrayList;
import java.util.List;

import static lphybeast.tobeast.TaxaUtils.getTaxonSet;

/**
 * Maps the LPhy skyline generator to the BEAST {@link CalibratedBirthDeathSkylineModel}. Each specified
 * rate array becomes a {@link SkylineParameter} (values + optional changeTimes, timesAreAges=true to
 * match the generator's present-to-past ordering). The origin / conditioning / calibration wiring is
 * identical to {@link CalibratedCPPToBEAST}.
 */
public class CalibratedSkylineCPPToBEAST
        implements GeneratorToBEAST<CalibratedBirthDeathSkylineTree, CalibratedBirthDeathSkylineModel> {

    @Override
    public CalibratedBirthDeathSkylineModel generatorToBEAST(
            CalibratedBirthDeathSkylineTree generator, BEASTInterface value, BEASTContext context) {

        CalibratedBirthDeathSkylineModel model = new CalibratedBirthDeathSkylineModel();
        model.setInputValue("tree", value);

        boolean hasCalibrations = generator.getCalibrations() != null;
        boolean rootConditioned = !hasCalibrations || generator.getRootCondition();
        model.setInputValue("conditionOnRoot", rootConditioned);
        model.setInputValue("conditionOnCalibrations", hasCalibrations);

        if (!rootConditioned) {
            model.setInputValue("origin", new RealScalarParam<>(generator.getOrigin().value(), PositiveReal.INSTANCE));
        }

        Value<Double[]> changeTimes = generator.getChangeTimes();
        setRate(model, "birthRate", generator.getBirthRate(), changeTimes, context);
        setRate(model, "deathRate", generator.getDeathRate(), changeTimes, context);
        setRate(model, "diversificationRate", generator.getDiversificationRate(), changeTimes, context);
        setRate(model, "turnover", generator.getTurnover(), changeTimes, context);
        setRate(model, "reproductiveNumber", generator.getReproductiveNumber(), changeTimes, context);

        model.setInputValue("rho", context.getAsRealScalar(generator.getSamplingProb()));

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

    /** Build a SkylineParameter from a per-interval rate array + shared change times, present-to-past. */
    private void setRate(CalibratedBirthDeathSkylineModel model, String inputName,
                         Value<Double[]> values, Value<Double[]> changeTimes, BEASTContext context) {
        if (values == null) return;
        SkylineParameter sp = new SkylineParameter();
        sp.setInputValue("values", context.getAsRealTensor(values));
        if (changeTimes != null) {
            sp.setInputValue("changeTimes", context.getAsRealTensor(changeTimes));
        }
        sp.setInputValue("timesAreAges", true);   // generator arrays are ordered present-to-past
        sp.initAndValidate();
        model.setInputValue(inputName, sp);
    }

    @Override
    public Class<CalibratedBirthDeathSkylineTree> getGeneratorClass() {
        return CalibratedBirthDeathSkylineTree.class;
    }

    @Override
    public Class<CalibratedBirthDeathSkylineModel> getBEASTClass() {
        return CalibratedBirthDeathSkylineModel.class;
    }
}
