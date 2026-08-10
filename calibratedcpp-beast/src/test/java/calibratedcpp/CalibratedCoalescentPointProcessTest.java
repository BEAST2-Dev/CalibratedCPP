package calibratedcpp;

import beast.base.evolution.speciation.CalibratedBirthDeathModel;
import beast.base.evolution.speciation.CalibrationPoint;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.inference.distribution.Uniform;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.alignment.Taxon;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.parameter.RealScalarParam;
import calibration.CalibrationForest;
import calibration.CalibrationNode;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class CalibratedCoalescentPointProcessTest {

    private calibratedcpp.CalibratedBirthDeathModel cpp;
    private final Tree tree;
    TaxonSet taxaABC;
    TaxonSet taxaDE;
    TaxonSet taxaABCDE;
    TaxonSet taxaHI;
    TaxonSet taxaHIJ;
    TaxonSet taxaFGHIJ;
    private final List<TaxonSet> calibrations;
    private CalibratedCoalescentPointProcess rootConditionedCPP;

    public CalibratedCoalescentPointProcessTest() {
        tree = new TreeParser();
        tree.initByName("newick", "((((A:2,B:2):1,C:3):1,(D:1.5,E:1.5):2.5):2,((F:1,G:1):4,((H:0.5,I:0.5):2,J:2.5):2.5):1):0;",
                "adjustTipHeights", false,
                "IsLabelledNewick", true);


        cpp = new calibratedcpp.CalibratedBirthDeathModel();
        cpp.initByName("tree", tree,
                "origin", new RealScalarParam<>(6.5, Real.INSTANCE),
                "birthRate", new RealScalarParam<>(3.0, Real.INSTANCE),
                "deathRate", new RealScalarParam<>(2.0, Real.INSTANCE),
                "rho", new RealScalarParam<>(0.1, Real.INSTANCE)
        );

        rootConditionedCPP = new calibratedcpp.CalibratedBirthDeathModel();
        rootConditionedCPP.initByName("tree", tree,
                "conditionOnRoot", true,
                "birthRate", new RealScalarParam<>(3.0, Real.INSTANCE),
                "deathRate", new RealScalarParam<>(2.0, Real.INSTANCE),
                "rho", new RealScalarParam<>(0.1, Real.INSTANCE)
        );

        // Create taxa objects
        Taxon A = new Taxon();
        A.setID("A");
        Taxon B = new Taxon();
        B.setID("B");
        Taxon C = new Taxon();
        C.setID("C");
        Taxon D = new Taxon();
        D.setID("D");
        Taxon E = new Taxon();
        E.setID("E");
        Taxon F = new Taxon();
        F.setID("F");
        Taxon G = new Taxon();
        G.setID("G");
        Taxon H = new Taxon();
        H.setID("H");
        Taxon I = new Taxon();
        I.setID("I");
        Taxon J = new Taxon();
        J.setID("J");

        // Create TaxonSets
        taxaABC = new TaxonSet(Arrays.asList(A, B, C));
        taxaDE = new TaxonSet(Arrays.asList(D, E));
        taxaABCDE = new TaxonSet(Arrays.asList(A, B, C, D, E));
        taxaHI = new TaxonSet(Arrays.asList(H, I));
        taxaHIJ = new TaxonSet(Arrays.asList(H, I, J));
        taxaFGHIJ = new TaxonSet(Arrays.asList(F, G, H, I, J));

        calibrations = Arrays.asList(taxaDE, taxaABC, taxaABCDE, taxaHI, taxaHIJ, taxaFGHIJ);
    }

    @Test
    public void calculateUnConditionedTreeLogLikelihood() {

        int n = tree.getLeafNodeCount();

        double log_n_factorial = 0;
        for (int i = 1; i <= n; i++) {
            log_n_factorial += Math.log(i);
        }

        assertEquals(-25.05062 + (n - 1) * Math.log(2) - log_n_factorial,
                cpp.calculateUnConditionedTreeLogLikelihood(tree), 1e-4, "Unconditioned density of the tree is incorrect.");
    }

    @Test
    public void calibrationForest() {
        CalibrationForest calibrationForest = new CalibrationForest(calibrations);
        List<CalibrationNode> calibrationNodes = calibrationForest.getAllNodes();

        for (CalibrationNode node : calibrationNodes) {
            if (!node.isRoot) {
                assertTrue(node.parent.taxa.getTaxonCount() > node.taxa.getTaxonCount());
                assertTrue(node.parent.getCommonAncestor(tree).getHeight() > node.getCommonAncestor(tree).getHeight());
            }
        }

        int n_roots = 0;
        for (CalibrationNode node : calibrationNodes) {
            n_roots += node.isRoot ? 1 : 0;
        }

        assertEquals(2, n_roots, "Calibration forest should have 2 roots but was " + n_roots);
    }

    @Test
    public void calculateLogMarginalDensityOfCalibrations() {
        CalibrationForest calibrationForest = new CalibrationForest(calibrations);
        List<CalibrationNode> calibrationNodes = calibrationForest.getAllNodes();

        CalibrationNode cpFGHIJ = CalibrationNode.getByTaxa(calibrationNodes, taxaFGHIJ);
        CalibrationNode cpABCDE = CalibrationNode.getByTaxa(calibrationNodes, taxaABCDE);

        double labellings = 0;
        for (int i = 1; i <= 5; i++) {
            labellings += Math.log(i);
        }
        labellings = 2 * labellings;
        for (int i = 1; i <= 10; i++) {
            labellings -= Math.log(i);
        }

        assert cpABCDE != null;
        assert cpFGHIJ != null;
        assertEquals(rootConditionedCPP.computeCalibrationDensity(tree, cpABCDE) +
                        rootConditionedCPP.computeCalibrationDensity(tree, cpFGHIJ) +
                        Math.log(2.0) + rootConditionedCPP.calculateLogNodeAgeDensity(6.0) + labellings,
                rootConditionedCPP.calculateMarginalLogDensityOfCalibrations(tree, calibrationForest), 1e-4,
                "Marginal density of the calibrations and root is incorrect.");

        assertEquals(cpp.computeCalibrationDensity(tree, cpABCDE) + cpp.computeCalibrationDensity(tree, cpFGHIJ) +
                        Math.log(2.0) + cpp.calculateLogNodeAgeCDF(6.5) + Math.log1p(-Math.exp(cpp.calculateLogNodeAgeCDF(5.0) - cpp.calculateLogNodeAgeCDF(6.5))) +
                        labellings,
                cpp.calculateMarginalLogDensityOfCalibrations(tree, calibrationForest), 1e-4,
                "MarginalDensity of the calibrations is incorrect.");
    }

    @Test
    public void computeCalibrationDensity() {
        CalibrationForest calibrationForest = new CalibrationForest(calibrations);
        List<CalibrationNode> calibrationNodes = calibrationForest.getAllNodes();

        CalibrationNode cpHI = CalibrationNode.getByTaxa(calibrationNodes, taxaHI);
        CalibrationNode cpDE = CalibrationNode.getByTaxa(calibrationNodes, taxaDE);
        CalibrationNode cpABC = CalibrationNode.getByTaxa(calibrationNodes, taxaABC);
        CalibrationNode cpHIJ = CalibrationNode.getByTaxa(calibrationNodes, taxaHIJ);
        CalibrationNode cpFGHIJ = CalibrationNode.getByTaxa(calibrationNodes, taxaFGHIJ);
        CalibrationNode cpABCDE = CalibrationNode.getByTaxa(calibrationNodes, taxaABCDE);

        assert cpHI != null;
        assertEquals(cpp.calculateLogNodeAgeDensity(0.5),
                cpp.computeCalibrationDensity(tree, cpHI), 1e-6, "Density for calibration HI is incorrect.");
        assert cpDE != null;
        assertEquals(cpp.calculateLogNodeAgeDensity(1.5),
                cpp.computeCalibrationDensity(tree, cpDE), 1e-6, "Density for calibration DE is incorrect.");
        assert cpABC != null;
        assertEquals(cpp.calculateLogNodeAgeDensity(3.0) + cpp.calculateLogNodeAgeCDF(3.0) + Math.log(2.0),
                cpp.computeCalibrationDensity(tree, cpABC), 1e-6, "Density for calibration ABC is incorrect.");
        assert cpHIJ != null;
        assertEquals(cpp.calculateLogNodeAgeDensity(0.5) + cpp.calculateLogNodeAgeDensity(2.5) + Math.log(2.0)
                        - (Math.log(3.0)),
                cpp.computeCalibrationDensity(tree, cpHIJ), 1e-6, "Density for calibration HIJ is incorrect.");
        assert cpABCDE != null;
        assertEquals(cpp.calculateLogNodeAgeDensity(3.0) + cpp.calculateLogNodeAgeCDF(3.0) + Math.log(2.0) +
                        cpp.calculateLogNodeAgeDensity(1.5) +
                        cpp.calculateLogNodeAgeDensity(4.0) +
                        2 * Math.log(2.0) - (Math.log(5.0) + Math.log(4.0)),
                cpp.computeCalibrationDensity(tree, cpABCDE), 1e-6, "Density for calibration ABCDE is incorrect.");
        assert cpFGHIJ != null;
        assertEquals(cpp.calculateLogNodeAgeDensity(0.5) + cpp.calculateLogNodeAgeDensity(2.5) + Math.log(2.0)
                        - (Math.log(3.0)) + // density of HIJ
                        cpp.calculateLogNodeAgeDensity(5.0) +
                        Math.log(4.0 * (Math.exp(cpp.calculateLogNodeAgeCDF(5.0)) - Math.exp(cpp.calculateLogNodeAgeCDF(2.5)))
                                + 2.0 * Math.exp(cpp.calculateLogNodeAgeCDF(5.0))) + Math.log(2.0) - (Math.log(5.0) + Math.log(4.0)),
                cpp.computeCalibrationDensity(tree, cpFGHIJ), 1e-6, "Density for calibration FGHIJ is incorrect.");
    }

    @Test
    public void calculateTreeLogLikelihood() {
        CalibrationForest calibrationForest = new CalibrationForest(calibrations);
        cpp = new calibratedcpp.CalibratedBirthDeathModel();
        cpp.initByName("tree", tree,
                "origin", new RealScalarParam<>(6.5, Real.INSTANCE),
                "calibrations", calibrations,
                "birthRate", new RealScalarParam<>(3.0, Real.INSTANCE),
                "deathRate", new RealScalarParam<>(2.0, Real.INSTANCE),
                "rho", new RealScalarParam<>(0.1, Real.INSTANCE));

        assertEquals(cpp.calculateUnConditionedTreeLogLikelihood(tree) -
                        cpp.calculateMarginalLogDensityOfCalibrations(tree, calibrationForest) -
                        Math.log1p(-Math.exp(cpp.calculateLogNodeAgeCDF(6.5))),
                cpp.calculateTreeLogLikelihood(tree), 1e-4, "Tree log likelihood incorrect.");

        rootConditionedCPP = new calibratedcpp.CalibratedBirthDeathModel();
        rootConditionedCPP.initByName("tree", tree,
                "origin", new RealScalarParam<>(6.5, Real.INSTANCE),
                "calibrations", calibrations,
                "conditionOnRoot", true,
                "conditionOnCalibrations", true,
                "birthRate", new RealScalarParam<>(3.0, Real.INSTANCE),
                "deathRate", new RealScalarParam<>(2.0, Real.INSTANCE),
                "rho", new RealScalarParam<>(0.1, Real.INSTANCE));

        assertEquals(rootConditionedCPP.calculateUnConditionedTreeLogLikelihood(tree) -
                        rootConditionedCPP.calculateMarginalLogDensityOfCalibrations(tree, calibrationForest) -
                        Math.log1p(-Math.exp(rootConditionedCPP.calculateLogNodeAgeCDF(6.0))),
                rootConditionedCPP.calculateTreeLogLikelihood(tree), 1e-4, "Tree log likelihood incorrect.");

        CalibratedCoalescentPointProcess nonmonophyleticCPP = new calibratedcpp.CalibratedBirthDeathModel();
        TaxonSet badClade = new TaxonSet(Arrays.asList(new Taxon("A"), new Taxon("B"), new Taxon("J")));
        nonmonophyleticCPP.initByName("tree", tree,
                "origin", new RealScalarParam<>(6.5, Real.INSTANCE),
                "calibrations", badClade,
                "conditionOnRoot", true,
                "conditionOnCalibrations", true,
                "birthRate", new RealScalarParam<>(3.0, Real.INSTANCE),
                "deathRate", new RealScalarParam<>(2.0, Real.INSTANCE),
                "rho", new RealScalarParam<>(0.1, Real.INSTANCE));

        assertEquals(Double.NEGATIVE_INFINITY, nonmonophyleticCPP.calculateTreeLogLikelihood(tree), 1e-4);

        cpp = new calibratedcpp.CalibratedBirthDeathModel();
        CalibratedBirthDeathModel heled_and_drummond = new CalibratedBirthDeathModel();
        List<TaxonSet> calibrationsClades = new ArrayList<>();
        List<CalibrationPoint> calibrationPoints = new ArrayList<>();
        Tree tree = new TreeParser();
        String newick = "(((((((leaf_1:0.0053482383689566,leaf_2:0.0053482383689566):0.00513631470300003,(leaf_3:0.00259161943948299,leaf_4:0.00259161943948299):0.00789293363247364):0.408273300539468,((leaf_5:0.343142440158881,leaf_6:0.343142440158881):0.0323374070437082,(leaf_7:0.214886773898683,leaf_8:0.214886773898683):0.160593073303907):0.0432780064088353):0.0717911514651002,(((leaf_9:0.00124139759074146,leaf_10:0.00124139759074146):0.00116279666983845,(leaf_11:0.00073018772695743,leaf_12:0.00073018772695743):0.00167400653362248):0.0015081805093331,((leaf_13:0.000432465313251912,leaf_14:0.000432465313251912):0.00166092495567772,(leaf_15:0.000925903870300905,leaf_16:0.000925903870300905):0.00116748639862872):0.00181898450098339):0.486636630306612):0.45882358903282,((((leaf_17:0.0250415381283358,leaf_18:0.0250415381283358):0.0360078258337303,(leaf_19:0.00501374371563791,leaf_20:0.00501374371563791):0.0560356202464282):0.410764850132822,((leaf_21:0.206990382368896,leaf_22:0.206990382368896):0.0153716271567664,(leaf_23:0.181968552131459,leaf_24:0.181968552131459):0.0403934573942037):0.249452204569225):0.296669902402815,(((leaf_25:0.0622075176933535,leaf_26:0.0622075176933535):0.0377742785779953,(leaf_27:0.0888783956157675,leaf_28:0.0888783956157675):0.0111034006555813):0.400930636489299,((leaf_29:0.142880344938527,leaf_30:0.142880344938527):0.171019789692806,(leaf_31:0.285445724395697,leaf_32:0.285445724395697):0.0284544102356365):0.187012298129315):0.267571683737055):0.180888477611642):1.04530810834413,(((((leaf_33:0.000548778697062376,leaf_34:0.000548778697062376):0.00595711326807519,(leaf_35:0.0061631439553151,leaf_36:0.0061631439553151):0.000342748009822466):0.00672585735347118,((leaf_37:0.00269852623305871,leaf_38:0.00269852623305871):0.000821969486263879,(leaf_39:0.00337591049444399,leaf_40:0.00337591049444399):0.000144585224878599):0.00971125359928615):0.252005129636834,(((leaf_41:0.0036007729390144,leaf_42:0.0036007729390144):0.00519531868051424,(leaf_43:0.0033142990662333,leaf_44:0.0033142990662333):0.00548179255329534):0.001071381527686,((leaf_45:0.00936098774886417,leaf_46:0.00936098774886417):0.000425700724809363,(leaf_47:0.00494817329477705,leaf_48:0.00494817329477705):0.00483851517889648):8.07846735411122e-05):0.255369405808228):0.0908009489275084,((((leaf_49:0.00959976713057155,leaf_50:0.00959976713057155):0.00193663653154847,(leaf_51:0.0111719537895163,leaf_52:0.0111719537895163):0.000364449872603702):0.00182190329717111,((leaf_53:0.0122817793969458,leaf_54:0.0122817793969458):0.00106498051034262,(leaf_55:0.00209579550389016,leaf_56:0.00209579550389016):0.0112509644033983):1.15470520026858e-05):0.25655489307181,(((leaf_57:0.00130176990143963,leaf_58:0.00130176990143963):0.000498934724185128,(leaf_59:0.000104709377330967,leaf_60:0.000104709377330967):0.00169599524829379):0.0339763114629944,((leaf_61:0.00504640237017586,leaf_62:0.00504640237017586):0.0247332543001393,(leaf_63:0.0113779123438897,leaf_64:0.0113779123438897):0.0184017443264255):0.00599735941830395):0.234136183942482):0.08612462785185):1.63864287457053):5.00531929754652,((((((leaf_65:0.0993376277852445,leaf_66:0.0993376277852445):0.0277324841725708,(leaf_67:0.0427918520979902,leaf_68:0.0427918520979902):0.0842782598598251):0.282930426652657,((leaf_69:0.197699420460351,leaf_70:0.197699420460351):0.0568358111006935,(leaf_71:0.0304582405965438,leaf_72:0.0304582405965438):0.2240769909645):0.155465307049428):0.249216777826801,(((leaf_73:0.000504410793363742,leaf_74:0.000504410793363742):0.0127404238089494,(leaf_75:0.0108515809933659,leaf_76:0.0108515809933659):0.0023932536089473):0.306026405292782,((leaf_77:0.0489001165334496,leaf_78:0.0489001165334496):0.0525987905835901,(leaf_79:0.0229843977360284,leaf_80:0.0229843977360284):0.0785145093810113):0.217772332778056):0.339946076542177):0.673268058849571,((((leaf_81:0.0791057071654503,leaf_82:0.0791057071654503):0.21888372219784,(leaf_83:0.186062876711128,leaf_84:0.186062876711128):0.111926552652162):0.00159840422063279,((leaf_85:0.134196290056542,leaf_86:0.134196290056542):0.143282492889759,(leaf_87:0.145301213568814,leaf_88:0.145301213568814):0.132177569377487):0.0221090506376222):0.565101906463337,(((leaf_89:0.0282689879745863,leaf_90:0.0282689879745863):0.0432786600984429,(leaf_91:0.0220527755129882,leaf_92:0.0220527755129882):0.049494872560041):0.0519719182668528,((leaf_93:0.091523652919298,leaf_94:0.091523652919298):0.0147940722108269,(leaf_95:0.0454212577456006,leaf_96:0.0454212577456006):0.0608964673845243):0.0172018412097571):0.741170173707378):0.467795635239583):4.11476138670573,(((((leaf_97:0.0775272810841276,leaf_98:0.0775272810841276):0.0146643419457146,(leaf_99:0.00603036413870096,leaf_100:0.00603036413870096):0.0861612588911412):0.0419611316921203,((leaf_101:0.084473471609642,leaf_102:0.084473471609642):0.0111589837893506,(leaf_103:0.0692929707261133,leaf_104:0.0692929707261133):0.0263394846728793):0.0385202993229698):0.838567645728143,(((leaf_105:0.0311156832535027,leaf_106:0.0311156832535027):0.00868929046436534,(leaf_107:0.00804186975004645,leaf_108:0.00804186975004645):0.0317631039678215):0.452535942306971,((leaf_109:0.0408987016152422,leaf_110:0.0408987016152422):0.00856266252056692,(leaf_111:0.00107066206097734,leaf_112:0.00107066206097734):0.0483907020748318):0.44287955188903):0.480379484425267):1.2984854592467,((((leaf_113:0.0406253345706209,leaf_114:0.0406253345706209):0.0333254809351503,(leaf_115:0.0182626770973543,leaf_116:0.0182626770973543):0.0556881384084169):0.0939277980361145,((leaf_117:0.0196800346285998,leaf_118:0.0196800346285998):0.0704843554687012,(leaf_119:0.0578484251316978,leaf_120:0.0578484251316978):0.0323159649656032):0.0777142234445846):0.755916334078699,(((leaf_121:0.301802373784873,leaf_122:0.301802373784873):0.11202585744297,(leaf_123:0.201427100746006,leaf_124:0.201427100746006):0.212401130481838):0.0356259544543635,((leaf_125:0.106412613757358,leaf_126:0.106412613757358):0.314952453116538,(leaf_127:0.348804165150462,leaf_128:0.348804165150462):0.0725609017234336):0.0280891188083107):0.474340761938378):1.34741091207622):3.17604090229577):1.55275323800743);";
        tree.initByName("newick", newick,
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        TaxonSet taxa = tree.getTaxonset();
        ParametricDistribution ageDist = new Uniform();
        ageDist.initByName("upper", 1.0, "lower", 0.0);

        CalibrationPoint rootCalibrationPoint = new CalibrationPoint();
        ParametricDistribution rootAgeDist = new Uniform();
        rootAgeDist.initByName("lower", 6.5, "upper", 7.5);
        rootCalibrationPoint.initByName("taxonset", taxa,
                "distr", rootAgeDist);
        calibrationPoints.add(rootCalibrationPoint);

        for (int i = 0; i < 3; i++) {
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
        double logNfactorial = 0.0;
        for (int i = 1; i <= taxa.getTaxonCount(); i++) {
            logNfactorial += Math.log(i);
        }

        RealScalarParam<Real> cppTurnover = new RealScalarParam<>(0.0, Real.INSTANCE);
        RealScalarParam<Real> cppBirthRate = new RealScalarParam<>(2.0, Real.INSTANCE);
        RealScalarParam<Real> cppRho = new RealScalarParam<>(1.0, Real.INSTANCE);
        RealParameter hdTurnover = new RealParameter("0.0");
        RealParameter hdBirthRate = new RealParameter("2.0");
        RealParameter hdRho = new RealParameter("1.0");

        cpp.initByName("tree", tree,
                "birthRate", cppBirthRate,
                "turnover", cppTurnover,
                "rho", cppRho,
                "calibrations", calibrationsClades,
                "conditionOnRoot", true);
        heled_and_drummond.initByName("tree", tree,
                "calibrations", calibrationPoints,
                "birthRate", hdBirthRate,
                "relativeDeathRate", hdTurnover,
                "sampleProbability", hdRho);
        assertEquals(cpp.calculateTreeLogLikelihood(tree),
                heled_and_drummond.calculateTreeLogLikelihood(tree) + logNfactorial, 1e-10,
                "Calibrated CPP tree log likelihood does not math Heled and Drummond Tree log likelihood.");
    }

    @Test
    public void heledAndDrummondComparison() {
        cpp = new calibratedcpp.CalibratedBirthDeathModel();
        CalibratedBirthDeathModel heled_and_drummond = new CalibratedBirthDeathModel();
        RealParameter hdRho = new RealParameter("1.0");
        RealParameter hdTurnover = new RealParameter("0.0");
        RealParameter hdBirthRate = new RealParameter("1.0");
        RealScalarParam<Real> cppRho = new RealScalarParam<>(1.0, Real.INSTANCE);
        RealScalarParam<Real> cppTurnover = new RealScalarParam<>(0.0, Real.INSTANCE);
        RealScalarParam<Real> cppBirthRate = new RealScalarParam<>(1.0, Real.INSTANCE);
        Tree tree = new TreeParser();
        tree.initByName("newick", "((((((((((((leaf_36:0.0257564750822701,leaf_40:0.0257564750822701):0.0408766525906903,leaf_8:0.0666331276729604):0.33784921072051,(((leaf_69:0.161465975796346,leaf_41:0.161465975796346):0.167664175344224,(leaf_82:0.0546036773493572,leaf_14:0.0546036773493572):0.274526473791213):0.0628134879891438,leaf_76:0.391943639129714):0.0125386992637559):0.40429517583965,leaf_97:0.808777514233121):0.0797578949429036,((((leaf_25:0.0265891826814647,leaf_93:0.0265891826814647):0.408929550049297,(leaf_60:0.0668468415852558,leaf_4:0.0668468415852558):0.368671891145506):0.103624686904037,((leaf_3:0.00856185800216441,leaf_44:0.00856185800216441):0.261366839599589,(leaf_32:0.0900877169305987,leaf_19:0.0900877169305987):0.179840980671155):0.269214722033045):0.246537911198444,((leaf_23:0.434973653613092,((((leaf_24:0.0441479803365146,leaf_27:0.0441479803365146):0.19461883594281,((leaf_39:0.0629976251906552,leaf_64:0.0629976251906552):0.125878956425076,leaf_35:0.188876581615731):0.0498902346635939):0.0920611869103335,leaf_84:0.330828003189658):0.0938282048312802,(leaf_66:0.346435089735903,leaf_95:0.346435089735903):0.0782211182850354):0.0103174455921539):0.204772868083289,(leaf_6:0.555436362376936,(leaf_58:0.165879919804176,(leaf_57:0.160837145962567,((leaf_98:0.00789646768020438,leaf_42:0.00789646768020438):0.00303307664476509,leaf_56:0.0109295443249695):0.149907601637597):0.00504277384160948):0.38955644257276):0.0843101593194453):0.14593480913686):0.102854078342782):0.0435164129237676,(((leaf_78:0.30658437310568,leaf_75:0.30658437310568):0.196584430365379,(((leaf_20:0.277245201712736,(leaf_16:0.0604484643090907,leaf_30:0.0604484643090907):0.216796737403645):0.0313505036545378,leaf_81:0.308595705367274):0.0551891550134793,(((((leaf_74:0.126414472170799,(leaf_31:0.0281822278476253,leaf_18:0.0281822278476253):0.0982322443231737):0.137886945558104,(leaf_77:0.206138359238225,leaf_80:0.206138359238225):0.0581630584906775):0.0378316760404452,leaf_9:0.302133093769348):0.0131602232339091,(leaf_43:0.0777015544992117,leaf_65:0.0777015544992117):0.237591762504045):0.020637879755451,(leaf_70:0.0398131209514894,leaf_1:0.0398131209514894):0.296118075807219):0.0278536636220449):0.139383943090305):0.293294993665746,leaf_12:0.796463797136805):0.135588024962987):0.12642423310104,((((leaf_17:0.336027213306863,leaf_21:0.336027213306863):0.160126882550966,((((leaf_83:0.0796139233613961,leaf_99:0.0796139233613961):0.0605244795635852,leaf_11:0.140138402924981):0.0120162169621776,leaf_91:0.152154619887159):0.0710648878320053,(leaf_59:0.118183002992202,(leaf_34:0.021817066471393,leaf_92:0.021817066471393):0.0963659365208089):0.105036504726962):0.272934588138665):0.00899466370199514,(leaf_2:0.467853907337959,leaf_15:0.467853907337959):0.0372948522218648):0.0138473059376347,leaf_10:0.518996065497459):0.539479989703373):0.0123539804084682,(((leaf_13:0.0340387079715563,leaf_63:0.0340387079715563):0.0134229850965503,leaf_33:0.0474616930681066):0.0477769201156261,leaf_96:0.0952386131837327):0.975591422425567):0.185759267653266,(leaf_100:0.141102866301929,leaf_71:0.141102866301929):1.11548643696064):0.289572462272809,(((leaf_62:0.288881440067835,leaf_29:0.288881440067835):0.246345551450539,((((leaf_85:0.128008687886708,leaf_38:0.128008687886708):0.0159728161929468,leaf_67:0.143981504079655):0.0784370353526523,leaf_94:0.222418539432307):0.241808948635236,((leaf_68:0.0958255274549801,leaf_28:0.0958255274549801):0.35486637290526,(leaf_73:0.20833072778014,leaf_7:0.20833072778014):0.2423611725801):0.0135355877073033):0.0709995034508305):0.494591261350899,(((leaf_26:0.165867740396315,(leaf_5:0.0157430030157796,leaf_79:0.0157430030157796):0.150124737380536):0.129225280882201,leaf_37:0.295093021278516):0.322675115432227,(leaf_86:0.428185040694139,(leaf_22:0.414453988225309,(leaf_61:0.176570020060729,leaf_72:0.176570020060729):0.23788396816458):0.0137310524688301):0.189583096016604):0.41205011615853):0.516343512666102):0.0763871229425035,(((leaf_54:0.133197288572896,leaf_55:0.133197288572896):0.276552982984382,leaf_52:0.409750271557278):1.09024972844272,(leaf_53:1.4230695138519,(((leaf_45:0.189963610986554,(leaf_51:0.151644855185039,(leaf_47:0.131318550225172,(leaf_50:0.0204962210877418,leaf_46:0.0204962210877418):0.11082232913743):0.0203263049598669):0.038318755801515):0.30109173844668,leaf_48:0.491055349433234):0.708944650566766,leaf_49:1.2):0.223069513851897):0.0769304861481035):0.122548888477878):0.624676687025677,((leaf_89:0.2621349808486,leaf_88:0.2621349808486):1.7378650191514,(leaf_90:1.03897286705537,leaf_87:1.03897286705537):0.961027132944629):0.247225575503555);",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);

        TaxonSet taxonSet = tree.getTaxonset();

        List<Taxon> taxonList1 = new ArrayList<>();
        for (int i = 45; i <= 51; i++) {
            String TaxonID = "leaf_" + i;
            Taxon taxon = new Taxon();
            taxon.setID(TaxonID);
            taxonList1.add(taxon);
        }
        TaxonSet taxonSet1 = new TaxonSet();
        taxonSet1.initByName("taxon", taxonList1);

        List<Taxon> taxonList2 = new ArrayList<>();
        for (int i = 45; i <= 55; i++) {
            String TaxonID = "leaf_" + i;
            Taxon taxon = new Taxon();
            taxon.setID(TaxonID);
            taxonList2.add(taxon);
        }
        TaxonSet taxonSet2 = new TaxonSet();
        taxonSet2.initByName("taxon", taxonList2);

        List<Taxon> taxonList3 = new ArrayList<>();
        for (int i = 87; i <= 90; i++) {
            String TaxonID = "leaf_" + i;
            Taxon taxon = new Taxon();
            taxon.setID(TaxonID);
            taxonList3.add(taxon);
        }
        TaxonSet taxonSet3 = new TaxonSet();
        taxonSet3.initByName("taxon", taxonList3);


        ParametricDistribution distr1 = new Uniform();
        distr1.initByName("lower", 1.0,
                "upper", 2.0);

        CalibrationPoint calibrationPoint1 = new CalibrationPoint();
        calibrationPoint1.initByName("taxonset", taxonSet1,
                "distr", distr1);

        ParametricDistribution distr2 = new Uniform();
        distr2.initByName("lower", 1.0,
                "upper", 2.0);

        CalibrationPoint calibrationPoint2 = new CalibrationPoint();
        calibrationPoint2.initByName("taxonset", taxonSet2,
                "distr", distr2);

        ParametricDistribution distr3 = new Uniform();
        distr3.initByName("lower", 1.5,
                "upper", 2.5);

        CalibrationPoint calibrationPoint3 = new CalibrationPoint();
        calibrationPoint3.initByName("taxonset", taxonSet3,
                "distr", distr3);

        ParametricDistribution distr4 = new Uniform();
        distr4.initByName("lower", 2.0,
                "upper", 3.0);

        CalibrationPoint calibrationPoint4 = new CalibrationPoint();
        calibrationPoint4.initByName("taxonset", taxonSet,
                "distr", distr4);

        List<TaxonSet> calibrationClades = new ArrayList<>();
        calibrationClades.add(taxonSet1);
        calibrationClades.add(taxonSet2);
        calibrationClades.add(taxonSet3);
        calibrationClades.add(taxonSet);

        List<CalibrationPoint> calibrationPoints = new ArrayList<>();
        calibrationPoints.add(calibrationPoint1);
        calibrationPoints.add(calibrationPoint2);
        calibrationPoints.add(calibrationPoint3);
        calibrationPoints.add(calibrationPoint4);

        cpp.initByName("birthRate", cppBirthRate,
                "turnover", cppTurnover,
                "rho", cppRho,
                "tree", tree,
                "calibrations", calibrationClades,
                "conditionOnRoot", true);

        heled_and_drummond.initByName("birthRate", hdBirthRate,
                "tree", tree,
                "relativeDeathRate", hdTurnover,
                "sampleProbability", hdRho,
                "calibrations", calibrationPoints);

        double logNFactorial = 0.0;
        for (int i = 1; i <= 100; i++) {
            logNFactorial += Math.log(i);
        }
        try (PrintWriter writer = new PrintWriter(new FileWriter("../validation/calibratedcpp/heled-and-drummond_validation_and_benchmark/validation_results.csv"))) {

            // Write the CSV Header
            writer.println("Iteration,BirthRate,CPP_LogLikelihood,HeledAndDrummond_LogLikelihood");

            for (int i = 0; i < 41; i++) {
                double incr = (double) i / 10;
                double currentBirthRate = 1.0 + incr;
                String birthRateValue = String.valueOf(currentBirthRate);

                cppBirthRate.set(currentBirthRate);
                cpp.updateParameters();

                hdBirthRate.setValue(currentBirthRate);

                double cppVal = cpp.calculateTreeLogLikelihood(tree);
                double heled_and_drummondVal = heled_and_drummond.calculateTreeLogLikelihood(tree) + logNFactorial;

                // Write the values to the CSV file
                // Format: Iteration, BirthRate, CPP Val, HD Val
                writer.printf("%d,%s,%.8f,%.8f%n", i, birthRateValue, cppVal, heled_and_drummondVal);

                assertEquals(cppVal, heled_and_drummondVal, 1e-9,
                        "Tree Log-Likelihoods do not match for birthRate(" + currentBirthRate + ")");
            }
        } catch (IOException e) {
            e.printStackTrace();
            // Optional: Fail the test if file writing fails
            throw new RuntimeException("Failed to write to CSV file", e);
        }
    }
}