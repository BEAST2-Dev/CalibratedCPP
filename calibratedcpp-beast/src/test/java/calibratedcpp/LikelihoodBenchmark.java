package calibratedcpp;

import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.speciation.CalibrationPoint;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.distribution.Uniform;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.parameter.RealScalarParam;
import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;

import org.openjdk.jmh.results.format.ResultFormatType;

@State(Scope.Thread) // Each thread gets its own instance of the state
@BenchmarkMode(Mode.AverageTime) // Measure average time per operation
@OutputTimeUnit(TimeUnit.MICROSECONDS) // Output results in microseconds (us)
@Warmup(iterations = 3, time = 1) // 3 warm-up iterations of 1 second each
@Measurement(iterations = 5, time = 1) // 5 measurement iterations of 1 second each
@Fork(1) // Run in a separate JVM process to ensure clean environment
public class LikelihoodBenchmark {

    @Param({"10"})
    public int nCalibrations;

    CalibratedBirthDeathModel cpp;
    beast.base.evolution.speciation.CalibratedBirthDeathModel heled_and_drummond;
    Tree tree;
    List<TaxonSet> calibrationsClades;
    List<CalibrationPoint> calibrationPoints;

    // This runs once before the benchmarks start
    // We do initialization here so it doesn't count towards the execution time
    @Setup(Level.Trial)
    public void setup() {
        cpp = new CalibratedBirthDeathModel();
        heled_and_drummond = new beast.base.evolution.speciation.CalibratedBirthDeathModel();
        calibrationsClades = new ArrayList<>();
        calibrationPoints = new ArrayList<>();
        tree = new TreeParser();
        String newick = "(((((((leaf_1:0.0053482383689566,leaf_2:0.0053482383689566):0.00513631470300003,(leaf_3:0.00259161943948299,leaf_4:0.00259161943948299):0.00789293363247364):0.408273300539468,((leaf_5:0.343142440158881,leaf_6:0.343142440158881):0.0323374070437082,(leaf_7:0.214886773898683,leaf_8:0.214886773898683):0.160593073303907):0.0432780064088353):0.0717911514651002,(((leaf_9:0.00124139759074146,leaf_10:0.00124139759074146):0.00116279666983845,(leaf_11:0.00073018772695743,leaf_12:0.00073018772695743):0.00167400653362248):0.0015081805093331,((leaf_13:0.000432465313251912,leaf_14:0.000432465313251912):0.00166092495567772,(leaf_15:0.000925903870300905,leaf_16:0.000925903870300905):0.00116748639862872):0.00181898450098339):0.486636630306612):0.45882358903282,((((leaf_17:0.0250415381283358,leaf_18:0.0250415381283358):0.0360078258337303,(leaf_19:0.00501374371563791,leaf_20:0.00501374371563791):0.0560356202464282):0.410764850132822,((leaf_21:0.206990382368896,leaf_22:0.206990382368896):0.0153716271567664,(leaf_23:0.181968552131459,leaf_24:0.181968552131459):0.0403934573942037):0.249452204569225):0.296669902402815,(((leaf_25:0.0622075176933535,leaf_26:0.0622075176933535):0.0377742785779953,(leaf_27:0.0888783956157675,leaf_28:0.0888783956157675):0.0111034006555813):0.400930636489299,((leaf_29:0.142880344938527,leaf_30:0.142880344938527):0.171019789692806,(leaf_31:0.285445724395697,leaf_32:0.285445724395697):0.0284544102356365):0.187012298129315):0.267571683737055):0.180888477611642):1.04530810834413,(((((leaf_33:0.000548778697062376,leaf_34:0.000548778697062376):0.00595711326807519,(leaf_35:0.0061631439553151,leaf_36:0.0061631439553151):0.000342748009822466):0.00672585735347118,((leaf_37:0.00269852623305871,leaf_38:0.00269852623305871):0.000821969486263879,(leaf_39:0.00337591049444399,leaf_40:0.00337591049444399):0.000144585224878599):0.00971125359928615):0.252005129636834,(((leaf_41:0.0036007729390144,leaf_42:0.0036007729390144):0.00519531868051424,(leaf_43:0.0033142990662333,leaf_44:0.0033142990662333):0.00548179255329534):0.001071381527686,((leaf_45:0.00936098774886417,leaf_46:0.00936098774886417):0.000425700724809363,(leaf_47:0.00494817329477705,leaf_48:0.00494817329477705):0.00483851517889648):8.07846735411122e-05):0.255369405808228):0.0908009489275084,((((leaf_49:0.00959976713057155,leaf_50:0.00959976713057155):0.00193663653154847,(leaf_51:0.0111719537895163,leaf_52:0.0111719537895163):0.000364449872603702):0.00182190329717111,((leaf_53:0.0122817793969458,leaf_54:0.0122817793969458):0.00106498051034262,(leaf_55:0.00209579550389016,leaf_56:0.00209579550389016):0.0112509644033983):1.15470520026858e-05):0.25655489307181,(((leaf_57:0.00130176990143963,leaf_58:0.00130176990143963):0.000498934724185128,(leaf_59:0.000104709377330967,leaf_60:0.000104709377330967):0.00169599524829379):0.0339763114629944,((leaf_61:0.00504640237017586,leaf_62:0.00504640237017586):0.0247332543001393,(leaf_63:0.0113779123438897,leaf_64:0.0113779123438897):0.0184017443264255):0.00599735941830395):0.234136183942482):0.08612462785185):1.63864287457053):5.00531929754652,((((((leaf_65:0.0993376277852445,leaf_66:0.0993376277852445):0.0277324841725708,(leaf_67:0.0427918520979902,leaf_68:0.0427918520979902):0.0842782598598251):0.282930426652657,((leaf_69:0.197699420460351,leaf_70:0.197699420460351):0.0568358111006935,(leaf_71:0.0304582405965438,leaf_72:0.0304582405965438):0.2240769909645):0.155465307049428):0.249216777826801,(((leaf_73:0.000504410793363742,leaf_74:0.000504410793363742):0.0127404238089494,(leaf_75:0.0108515809933659,leaf_76:0.0108515809933659):0.0023932536089473):0.306026405292782,((leaf_77:0.0489001165334496,leaf_78:0.0489001165334496):0.0525987905835901,(leaf_79:0.0229843977360284,leaf_80:0.0229843977360284):0.0785145093810113):0.217772332778056):0.339946076542177):0.673268058849571,((((leaf_81:0.0791057071654503,leaf_82:0.0791057071654503):0.21888372219784,(leaf_83:0.186062876711128,leaf_84:0.186062876711128):0.111926552652162):0.00159840422063279,((leaf_85:0.134196290056542,leaf_86:0.134196290056542):0.143282492889759,(leaf_87:0.145301213568814,leaf_88:0.145301213568814):0.132177569377487):0.0221090506376222):0.565101906463337,(((leaf_89:0.0282689879745863,leaf_90:0.0282689879745863):0.0432786600984429,(leaf_91:0.0220527755129882,leaf_92:0.0220527755129882):0.049494872560041):0.0519719182668528,((leaf_93:0.091523652919298,leaf_94:0.091523652919298):0.0147940722108269,(leaf_95:0.0454212577456006,leaf_96:0.0454212577456006):0.0608964673845243):0.0172018412097571):0.741170173707378):0.467795635239583):4.11476138670573,(((((leaf_97:0.0775272810841276,leaf_98:0.0775272810841276):0.0146643419457146,(leaf_99:0.00603036413870096,leaf_100:0.00603036413870096):0.0861612588911412):0.0419611316921203,((leaf_101:0.084473471609642,leaf_102:0.084473471609642):0.0111589837893506,(leaf_103:0.0692929707261133,leaf_104:0.0692929707261133):0.0263394846728793):0.0385202993229698):0.838567645728143,(((leaf_105:0.0311156832535027,leaf_106:0.0311156832535027):0.00868929046436534,(leaf_107:0.00804186975004645,leaf_108:0.00804186975004645):0.0317631039678215):0.452535942306971,((leaf_109:0.0408987016152422,leaf_110:0.0408987016152422):0.00856266252056692,(leaf_111:0.00107066206097734,leaf_112:0.00107066206097734):0.0483907020748318):0.44287955188903):0.480379484425267):1.2984854592467,((((leaf_113:0.0406253345706209,leaf_114:0.0406253345706209):0.0333254809351503,(leaf_115:0.0182626770973543,leaf_116:0.0182626770973543):0.0556881384084169):0.0939277980361145,((leaf_117:0.0196800346285998,leaf_118:0.0196800346285998):0.0704843554687012,(leaf_119:0.0578484251316978,leaf_120:0.0578484251316978):0.0323159649656032):0.0777142234445846):0.755916334078699,(((leaf_121:0.301802373784873,leaf_122:0.301802373784873):0.11202585744297,(leaf_123:0.201427100746006,leaf_124:0.201427100746006):0.212401130481838):0.0356259544543635,((leaf_125:0.106412613757358,leaf_126:0.106412613757358):0.314952453116538,(leaf_127:0.348804165150462,leaf_128:0.348804165150462):0.0725609017234336):0.0280891188083107):0.474340761938378):1.34741091207622):3.17604090229577):1.55275323800743);";
        tree.initByName("newick", newick,
                "IsLabelledNewick", true,
                "adjustTipHeights", false);

        TaxonSet taxa = tree.getTaxonset();
        ParametricDistribution ageDist = new Uniform();
        ageDist.initByName("upper", 1.0, "lower", 0.0);

        CalibrationPoint rootCalibrationPoint = new CalibrationPoint();
        ParametricDistribution rootAgeDist = new Uniform();
        rootAgeDist.initByName("lower", 6.5, "upper", 7.5);
        rootCalibrationPoint.initByName("taxonset", taxa,
                "distr", rootAgeDist);
        calibrationPoints.add(rootCalibrationPoint);

        for (int i = 0; i < nCalibrations; i++) {
            String t1Name = "leaf_" + (2 * i + 1);
            String t2Name = "leaf_" + (2 * (i + 1));

            TaxonSet taxonSet = new TaxonSet();
            List<Taxon> taxonList = new ArrayList<>();

            taxonList.add(taxa.getTaxon(t1Name));
            taxonList.add(taxa.getTaxon(t2Name));

            taxonSet.initByName("taxon", taxonList);

            calibrationsClades.add(taxonSet);

            CalibrationPoint calibrationPoint = new CalibrationPoint();
            calibrationPoint.initByName("taxonset", taxonSet, "distr", ageDist);
            calibrationPoints.add(calibrationPoint);
        }

        RealParameter turnoverRP = new RealParameter("0.0");
        RealParameter birthRateRP = new RealParameter("2.0");
        RealParameter rhoRP = new RealParameter("1.0");

        RealScalarParam<Real> turnover = new RealScalarParam<>(0.0, Real.INSTANCE);
        RealScalarParam<Real> birthRate = new RealScalarParam<>(2.0, Real.INSTANCE);
        RealScalarParam<Real> rho = new RealScalarParam<>(1.0, Real.INSTANCE);

        cpp.initByName("tree", tree,
                "birthRate", birthRate,
                "turnover", turnover,
                "rho", rho,
                "calibrations", calibrationsClades,
                "conditionOnRoot", true);
        heled_and_drummond.initByName("tree", tree,
                "calibrations", calibrationPoints,
                "birthRate", birthRateRP,
                "relativeDeathRate", turnoverRP,
                "sampleProbability", rhoRP);

        double cppVal = cpp.calculateTreeLogLikelihood(tree);
        double hdVal = heled_and_drummond.calculateTreeLogLikelihood(tree);
        double epsilon = 1e-9;

        double logNfactorial = 0.0;
        for (int i = 1; i <= taxa.getTaxonCount(); i++) {
            logNfactorial += Math.log(i);
        }

        if (Math.abs(cppVal - logNfactorial - hdVal) > epsilon) {
            throw new RuntimeException(String.format(
                    "Likelihood mismatch! CPP: %f, BD: %f", cppVal, hdVal
            ));
        }

        System.out.println("Verification passed: Both models produce " + cppVal);
    }

    // Benchmark for CalibratedBirthDeathModel
    @Benchmark
    public double measureHeledAndDrummond() {
        // 1. Change the parameter slightly to invalidate the cache
        // We toggle between 2.0 and 2.0000001
        double newVal = (cpp.birthRateInput.get().get() > 2.0) ? 2.0 : 2.0000001;
        heled_and_drummond.birthRateInput.get().setValue(newVal);

        // 2. Measure the calculation
        return heled_and_drummond.calculateTreeLogLikelihood(tree);
    }

    //     Benchmark for CalibratedCoalescentPointProcess
    @Benchmark
    public double measureCPP() {
        return cpp.calculateTreeLogLikelihood(tree);
    }

    // Standard main method to launch the benchmark
    public static void main(String[] args) throws Exception {
        // -------------------------------------------------------------
        // 3. GENERATE PARAMETERS PROGRAMMATICALLY
        // -------------------------------------------------------------
        // Create an array of strings: ["1", "2", "3", ... "50"]
        String[] calibrationCounts = IntStream.rangeClosed(1, 6)
                .mapToObj(String::valueOf)
                .toArray(String[]::new);

        Options opt = new OptionsBuilder()
                .include(LikelihoodBenchmark.class.getSimpleName())
                // Pass the sequence 1..50 to the "nCalibrations" field
                .param("nCalibrations", calibrationCounts)
                .resultFormat(ResultFormatType.CSV)
                .result("validation/calibratedcpp/heled-and-drummond_validation_and_benchmark/benchmark_results.csv")
                .build();

        new Runner(opt).run();
    }
}