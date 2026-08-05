package calibratedcpp.lphy.spi;

import calibratedcpp.lphy.prior.CalibrationFunction;
import calibratedcpp.lphy.prior.ConditionedMRCAPrior;
import calibratedcpp.lphy.prior.OffsetExponentialMRCA;
import calibratedcpp.lphy.prior.UniformMRCA;
import calibratedcpp.lphy.prior.toCalibrationArray;
import calibratedcpp.lphy.tree.CPPTree;
import calibratedcpp.lphy.tree.CalibratedAgeDependentExtinctionTree;
import calibratedcpp.lphy.tree.CalibratedBirthDeathSkylineTree;
import calibratedcpp.lphy.tree.CalibratedCPPTree;
import calibratedcpp.lphy.tree.ExpLifetime;
import calibratedcpp.lphy.tree.GammaLifetime;
import calibratedcpp.lphy.tree.WeibullLifetime;
import calibratedcpp.lphy.tree.WeibullMixtureLifetime;
import calibratedcpp.lphy.util.TruncatedLogNormal;
import lphy.base.spi.LPhyBaseImpl;
import lphy.core.model.BasicFunction;
import lphy.core.model.GenerativeDistribution;

import java.util.Arrays;
import java.util.List;

public class CalibratedcppImpl extends LPhyBaseImpl {
    public CalibratedcppImpl() {}

    @Override
    public List<Class<? extends GenerativeDistribution>> declareDistributions() {
        return Arrays.asList(
            CPPTree.class, CalibratedCPPTree.class,
                CalibratedBirthDeathSkylineTree.class, CalibratedAgeDependentExtinctionTree.class,
                TruncatedLogNormal.class, ConditionedMRCAPrior.class,
                UniformMRCA.class, OffsetExponentialMRCA.class
        );
    }

    @Override
    public List<Class<? extends BasicFunction>> declareFunctions() {
        return Arrays.asList(
                toCalibrationArray.class, CalibrationFunction.class,
                WeibullLifetime.class, GammaLifetime.class, ExpLifetime.class,
                WeibullMixtureLifetime.class
        );
    }

    @Override
    public String getExtensionName() {
        return "calibratedcpp lphy library";
    }
}
