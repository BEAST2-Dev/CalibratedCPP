open module calibratedcpp.lphybeast {
    requires lphy.beast;
    requires calibratedcpp.beast;
    requires calibratedcpp.lphy;
    requires beast.base;
    requires lphy.base;

    exports calibratedcpp.lphybeast;
    exports calibratedcpp.lphybeast.spi;
    exports calibratedcpp.lphybeast.tobeast.generators;

    provides lphybeast.spi.LPhyBEASTMapping with calibratedcpp.lphybeast.spi.LBcalibratedcppImpl;
}
