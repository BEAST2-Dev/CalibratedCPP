open module calibratedcpp.beast {
    requires beast.pkgmgmt;
    requires beast.base;
    requires beast.fx;
    requires hipparchus.core;

    // JavaFX modules used by CalibratedCPPInputEditor (beast.fx does not re-export these)
    requires javafx.graphics;
    requires javafx.controls;
    requires javafx.web;
    requires jdk.jsobject;
    requires org.apache.commons.statistics.distribution;
    requires commons.math3;

    exports calibratedcpp;
    exports calibration;
    exports calibrationprior;
    exports calibratedcpp.beauti;
    exports calibratedcpp.distribution;

    provides beastfx.app.inputeditor.InputEditor with
        calibratedcpp.beauti.CalibratedBirthDeathSkylineInputEditor,
        calibratedcpp.beauti.CalibratedAgeDependentExtinctionInputEditor,
        calibratedcpp.beauti.CalibrationDistributionInputEditor,
        calibratedcpp.beauti.SkylineParameterInputEditor;

    // CalibratedCoalescentPointProcess (abstract) and CalibrationNode/Forest (no default ctor)
    // are registered via version.xml only.
    provides beast.base.core.BEASTInterface with
        calibratedcpp.CalibratedBirthDeathModel,
        calibratedcpp.CalibratedBirthDeathSkylineModel,
        calibratedcpp.CalibratedAgeDependentExtinctionModel,
        calibratedcpp.CalibratedCPPTreeInitialiser,
        calibratedcpp.SkylineParameter,
        calibratedcpp.distribution.Weibull,
        calibratedcpp.distribution.Erlang,
        calibratedcpp.distribution.ScalarMixtureDistribution,
        calibratedcpp.distribution.WeibullMixture,
        calibratedcpp.operators.ChangeTimeOperator,
        calibrationprior.CalibrationPrior,
        calibrationprior.CalibrationDistribution,
        calibrationprior.CalibrationCladePrior,
        calibration.CalibrationForestParser;
}
