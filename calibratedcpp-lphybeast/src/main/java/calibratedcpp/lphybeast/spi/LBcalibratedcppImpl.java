package calibratedcpp.lphybeast.spi;

import beast.base.evolution.datatype.DataType;
import calibratedcpp.lphy.prior.Calibration;
import calibratedcpp.lphy.prior.CalibrationArray;
import calibratedcpp.lphy.prior.CalibrationFunction;
import calibratedcpp.lphy.prior.ConditionedMRCAPrior;
import calibratedcpp.lphy.prior.OffsetExponentialMRCA;
import calibratedcpp.lphy.prior.UniformMRCA;
import calibratedcpp.lphy.prior.toCalibrationArray;
import calibratedcpp.lphy.tree.ExpLifetime;
import calibratedcpp.lphy.tree.GammaLifetime;
import calibratedcpp.lphy.tree.LifetimeModel;
import calibratedcpp.lphy.tree.WeibullLifetime;
import calibratedcpp.lphy.tree.WeibullMixtureLifetime;
import calibratedcpp.lphy.util.TruncatedLogNormal;
import calibratedcpp.lphybeast.tobeast.generators.CPPToBEAST;
import calibratedcpp.lphybeast.tobeast.generators.CalibratedAgeDependentExtinctionCPPToBEAST;
import calibratedcpp.lphybeast.tobeast.generators.CalibratedCPPToBEAST;
import calibratedcpp.lphybeast.tobeast.generators.CalibratedSkylineCPPToBEAST;
import jebl.evolution.sequences.SequenceType;
import lphy.core.model.Generator;
import lphybeast.GeneratorToBEAST;
import lphybeast.ValueToBEAST;
import lphybeast.spi.LPhyBEASTMapping;

import java.util.List;
import java.util.Map;

public class LBcalibratedcppImpl implements LPhyBEASTMapping {
    @Override
    public List<Class<? extends GeneratorToBEAST>> getGeneratorToBEASTs() {
        return List.of(
                CalibratedCPPToBEAST.class, CPPToBEAST.class,
                CalibratedSkylineCPPToBEAST.class, CalibratedAgeDependentExtinctionCPPToBEAST.class
        );
    }

    @Override
    public List<Class<? extends ValueToBEAST>> getValuesToBEASTs() {
        return List.of();
    }

    @Override
    public List<Class<? extends Generator>> getExcludedGenerator() {
        return List.of(
                TruncatedLogNormal.class, ConditionedMRCAPrior.class, CalibrationFunction.class,
                UniformMRCA.class, OffsetExponentialMRCA.class, toCalibrationArray.class,
                WeibullLifetime.class, GammaLifetime.class, ExpLifetime.class,
                WeibullMixtureLifetime.class
        );
    }

    @Override
    public Map<SequenceType, DataType> getDataTypeMap() {
        return Map.of();
    }

    @Override
    public List<Class> getExcludedValueType() {
        return List.of(
                Calibration.class, CalibrationArray.class, LifetimeModel.class
        );
    }
}
