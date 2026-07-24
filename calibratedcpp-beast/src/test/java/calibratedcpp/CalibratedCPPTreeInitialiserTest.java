package calibratedcpp;

import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.Real;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.evolution.tree.MRCAPrior;
import beast.base.spec.inference.distribution.Exponential;
import beast.base.spec.inference.distribution.LogNormal;
import beast.base.spec.inference.distribution.Uniform;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.inference.Distribution;
import beast.base.parser.XMLParser;
import beast.base.util.Randomizer;
import calibration.CalibrationForest;
import calibrationprior.CalibrationCladePrior;
import calibrationprior.CalibrationPrior;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

public class CalibratedCPPTreeInitialiserTest {

    private static final String[] TAXA = {"A", "B", "C", "D", "E", "F"};

    private TaxonSet taxonSet(String... ids) {
        List<Taxon> taxa = new ArrayList<>();
        for (String id : ids) taxa.add(new Taxon(id));
        TaxonSet ts = new TaxonSet(taxa);
        ts.setID(String.join("_", ids));
        return ts;
    }

    private Tree initialTree() {
        Tree tree = new TreeParser();
        tree.initByName("newick", "(((A:1,B:1):1,C:2):1,(D:1,(E:1,F:1):1):1);",
                "IsLabelledNewick", true, "adjustTipHeights", true);
        return tree;
    }

    private CalibrationCladePrior cladePrior(TaxonSet ts, double lower, double upper) {
        CalibrationCladePrior p = new CalibrationCladePrior();
        p.initByName("taxa", ts,
                "lowerAge", new RealScalarParam<>(lower, NonNegativeReal.INSTANCE),
                "upperAge", new RealScalarParam<>(upper, NonNegativeReal.INSTANCE));
        p.setID("clade_" + ts.getID());
        return p;
    }

    /** MRCAPrior with a Uniform(lower, upper) age distribution — the common BEAUti calibration. */
    private MRCAPrior uniformMRCAPrior(Tree tree, TaxonSet ts, double lower, double upper) {
        Uniform uniform = new Uniform();
        uniform.initByName("lower", new RealScalarParam<>(lower, Real.INSTANCE),
                "upper", new RealScalarParam<>(upper, Real.INSTANCE));
        MRCAPrior prior = new MRCAPrior();
        prior.initByName("tree", tree, "taxonset", ts, "monophyletic", true, "distr", uniform);
        prior.setID("mrca_" + ts.getID());
        return prior;
    }

    /** MRCAPrior with a LogNormal age distribution — open upper support (0, +inf). */
    private MRCAPrior logNormalMRCAPrior(Tree tree, TaxonSet ts, double M, double S) {
        LogNormal ln = new LogNormal();
        ln.initByName("M", new RealScalarParam<>(M, Real.INSTANCE),
                "S", new RealScalarParam<>(S, PositiveReal.INSTANCE));
        MRCAPrior prior = new MRCAPrior();
        prior.initByName("tree", tree, "taxonset", ts, "monophyletic", true, "distr", ln);
        prior.setID("mrca_" + ts.getID());
        return prior;
    }

    private Set<String> leafSet(Node node) {
        Set<String> ids = new HashSet<>();
        collect(node, ids);
        return ids;
    }

    private void collect(Node node, Set<String> ids) {
        if (node.isLeaf()) ids.add(node.getID());
        else for (Node c : node.getChildren()) collect(c, ids);
    }

    private Node mrca(Tree tree, Set<String> target) {
        for (Node n : tree.getNodesAsArray()) {
            if (leafSet(n).equals(target)) return n;
        }
        return null;
    }

    private void assertValidTree(Tree tree) {
        // right number of taxa, all unique IDs, valid heights
        assertEquals(TAXA.length, tree.getLeafNodeCount());
        assertEquals(2 * TAXA.length - 1, tree.getNodeCount());
        Set<String> ids = new HashSet<>();
        for (Node leaf : tree.getExternalNodes()) ids.add(leaf.getID());
        assertEquals(TAXA.length, ids.size(), "leaf IDs must be unique and complete");
        for (Node n : tree.getInternalNodes()) {
            for (Node c : n.getChildren()) {
                assertTrue(n.getHeight() > c.getHeight(),
                        "parent must be strictly older than child");
            }
        }
    }

    @Test
    public void nestedCalibrationsAreMonophyleticAndWithinBounds() {
        Randomizer.setSeed(42);

        TaxonSet full = taxonSet(TAXA);
        TaxonSet ab = taxonSet("A", "B");
        TaxonSet abc = taxonSet("A", "B", "C");
        TaxonSet ef = taxonSet("E", "F");

        Tree tree = initialTree();
        CalibratedBirthDeathModel model = new CalibratedBirthDeathModel();
        model.initByName("tree", tree,
                "origin", new RealScalarParam<>(50.0, PositiveReal.INSTANCE),
                "birthRate", new RealScalarParam<>(2.0, Real.INSTANCE),
                "deathRate", new RealScalarParam<>(1.0, Real.INSTANCE),
                "rho", new RealScalarParam<>(0.1, Real.INSTANCE),
                "calibrations", new ArrayList<>(List.of(ab, abc, ef)));

        CalibrationForest forest = CalibrationForest.fromPriors(List.of(
                cladePrior(ab, 1.0, 3.0),
                cladePrior(abc, 4.0, 8.0),
                cladePrior(ef, 2.0, 5.0)));

        CalibratedCPPTreeInitialiser init = new CalibratedCPPTreeInitialiser();
        init.initByName("initial", tree, "taxonset", full,
                "model", model, "calibrationForest", forest);

        // Repeat to shake out sampling-dependent edge cases.
        for (int rep = 0; rep < 200; rep++) {
            init.initStateNodes();
            assertValidTree(tree);

            assertClade(tree, ab, 1.0, 3.0);
            assertClade(tree, abc, 4.0, 8.0);
            assertClade(tree, ef, 2.0, 5.0);

            // nesting: {A,B} MRCA younger than {A,B,C} MRCA
            double abAge = mrca(tree, Set.of("A", "B")).getHeight();
            double abcAge = mrca(tree, Set.of("A", "B", "C")).getHeight();
            assertTrue(abAge < abcAge, "nested clade must be younger");

            // root below origin
            assertTrue(tree.getRoot().getHeight() < 50.0, "root must be below origin");
        }
    }

    private void assertClade(Tree tree, TaxonSet ts, double lower, double upper) {
        Set<String> target = new HashSet<>(ts.getTaxaNames());
        Node mrca = mrca(tree, target);
        assertNotNull(mrca, "clade " + ts.getID() + " is not monophyletic");
        double age = mrca.getHeight();
        assertTrue(age >= lower - 1e-6 && age <= upper + 1e-6,
                "clade " + ts.getID() + " age " + age + " outside [" + lower + "," + upper + "]");
    }

    @Test
    public void mrcaPriorBoundsAreRespected() {
        Randomizer.setSeed(123);

        TaxonSet full = taxonSet(TAXA);
        TaxonSet ab = taxonSet("A", "B");
        TaxonSet abc = taxonSet("A", "B", "C");
        TaxonSet ef = taxonSet("E", "F");

        Tree tree = initialTree();
        CalibratedBirthDeathModel model = new CalibratedBirthDeathModel();
        model.initByName("tree", tree,
                "origin", new RealScalarParam<>(200.0, PositiveReal.INSTANCE),
                "birthRate", new RealScalarParam<>(2.0, Real.INSTANCE),
                "deathRate", new RealScalarParam<>(1.0, Real.INSTANCE),
                "rho", new RealScalarParam<>(0.1, Real.INSTANCE),
                "calibrations", new ArrayList<>(List.of(ab, abc, ef)));

        // Tight Uniform bounds (nested), plus a LogNormal (open upper) on {E,F}.
        List<MRCAPrior> constraints = List.of(
                uniformMRCAPrior(tree, ab, 10.0, 12.0),
                uniformMRCAPrior(tree, abc, 40.0, 45.0),
                logNormalMRCAPrior(tree, ef, 3.0, 0.5)); // median e^3 ≈ 20

        CalibratedCPPTreeInitialiser init = new CalibratedCPPTreeInitialiser();
        init.initByName("initial", tree, "taxonset", full, "model", model,
                "constraint", new ArrayList<>(constraints));

        for (int rep = 0; rep < 200; rep++) {
            init.initStateNodes();
            assertValidTree(tree);

            // Uniform supports must be honoured exactly.
            assertClade(tree, ab, 10.0, 12.0);
            assertClade(tree, abc, 40.0, 45.0);

            // LogNormal: open upper, support (0, +inf) — must at least be monophyletic and positive,
            // and (being nested under nothing here) below the root.
            Node efMrca = mrca(tree, Set.of("E", "F"));
            assertNotNull(efMrca, "{E,F} must be monophyletic");
            assertTrue(efMrca.getHeight() > 0.0);

            double abAge = mrca(tree, Set.of("A", "B")).getHeight();
            double abcAge = mrca(tree, Set.of("A", "B", "C")).getHeight();
            assertTrue(abAge < abcAge, "nested {A,B} must be younger than {A,B,C}");
            assertTrue(tree.getRoot().getHeight() < 200.0, "root below origin");
        }
    }

    @Test
    public void boundsUnpackedFromCalibrationDistributionWrapper() {
        Randomizer.setSeed(99);

        TaxonSet full = taxonSet(TAXA);
        TaxonSet ab = taxonSet("A", "B");
        TaxonSet ef = taxonSet("E", "F");

        Tree tree = initialTree();
        CalibratedBirthDeathModel model = new CalibratedBirthDeathModel();
        model.initByName("tree", tree,
                "origin", new RealScalarParam<>(100.0, PositiveReal.INSTANCE),
                "birthRate", new RealScalarParam<>(2.0, Real.INSTANCE),
                "deathRate", new RealScalarParam<>(1.0, Real.INSTANCE),
                "rho", new RealScalarParam<>(0.1, Real.INSTANCE),
                "calibrations", new ArrayList<>(List.of(ab, ef)));

        // Wrap two MRCAPriors in a CalibrationDistribution (CompoundDistribution) as BEAUti does.
        calibrationprior.CalibrationDistribution wrapper = new calibrationprior.CalibrationDistribution();
        wrapper.initByName("distribution", new ArrayList<>(List.of(
                uniformMRCAPrior(tree, ab, 15.0, 18.0),
                uniformMRCAPrior(tree, ef, 6.0, 9.0))));

        CalibratedCPPTreeInitialiser init = new CalibratedCPPTreeInitialiser();
        init.initByName("initial", tree, "taxonset", full, "model", model,
                "calibration", wrapper);

        for (int rep = 0; rep < 100; rep++) {
            init.initStateNodes();
            assertValidTree(tree);
            assertClade(tree, ab, 15.0, 18.0);
            assertClade(tree, ef, 6.0, 9.0);
            assertTrue(tree.getRoot().getHeight() < 100.0);
        }
    }

    /**
     * End-to-end check of the wiring the BEAUti template emits: parse a generated-style BEAST XML
     * (Tree + model + CalibrationDistribution(MRCAPriors) + the initialiser as an init object,
     * connected by idref) through the real XMLParser, then confirm the initialiser wrote a starting
     * tree that gives both the calibrated-CPP prior and the MRCA priors a finite starting density.
     */
    @Test
    public void generatedXmlWiringInitialisesFiniteStartingState() throws Exception {
        Randomizer.setSeed(2026);

        String P = "beast.base.spec.inference.parameter.RealScalarParam";
        String U = "beast.base.spec.inference.distribution.Uniform";
        String M = "beast.base.spec.evolution.tree.MRCAPrior";

        String xml =
            "<beast version='2.7'>\n" +
            "  <taxonset id='fullTaxa' spec='beast.base.evolution.alignment.TaxonSet'>\n" +
            "    <taxon id='A' spec='beast.base.evolution.alignment.Taxon'/>\n" +
            "    <taxon id='B' spec='beast.base.evolution.alignment.Taxon'/>\n" +
            "    <taxon id='C' spec='beast.base.evolution.alignment.Taxon'/>\n" +
            "    <taxon id='D' spec='beast.base.evolution.alignment.Taxon'/>\n" +
            "    <taxon id='E' spec='beast.base.evolution.alignment.Taxon'/>\n" +
            "    <taxon id='F' spec='beast.base.evolution.alignment.Taxon'/>\n" +
            "  </taxonset>\n" +
            "  <tree id='Tree' spec='beast.base.evolution.tree.Tree' taxonset='@fullTaxa'/>\n" +
            "  <taxonset id='ab' spec='beast.base.evolution.alignment.TaxonSet'>\n" +
            "    <taxon idref='A'/><taxon idref='B'/></taxonset>\n" +
            "  <taxonset id='ef' spec='beast.base.evolution.alignment.TaxonSet'>\n" +
            "    <taxon idref='E'/><taxon idref='F'/></taxonset>\n" +
            "  <distribution id='cpp' spec='calibratedcpp.CalibratedBirthDeathModel' conditionOnCalibrations='true' conditionOnRoot='true'>\n" +
            "    <tree idref='Tree'/>\n" +
            "    <birthRate spec='" + P + "' domain='PositiveReal' value='1.0'/>\n" +
            "    <deathRate spec='" + P + "' domain='NonNegativeReal' value='0.2'/>\n" +
            "    <rho spec='" + P + "' domain='UnitInterval' value='0.5'/>\n" +
            "    <calibrations idref='ab'/>\n" +
            "    <calibrations idref='ef'/>\n" +
            "  </distribution>\n" +
            "  <distribution id='calib' spec='calibrationprior.CalibrationDistribution'>\n" +
            "    <distribution spec='" + M + "' monophyletic='true' tree='@Tree' taxonset='@ab'>\n" +
            "      <distr spec='" + U + "'>\n" +
            "        <lower spec='" + P + "' domain='Real' value='1.5'/>\n" +
            "        <upper spec='" + P + "' domain='Real' value='2.0'/>\n" +
            "      </distr>\n" +
            "    </distribution>\n" +
            "    <distribution spec='" + M + "' monophyletic='true' tree='@Tree' taxonset='@ef'>\n" +
            "      <distr spec='" + U + "'>\n" +
            "        <lower spec='" + P + "' domain='Real' value='0.8'/>\n" +
            "        <upper spec='" + P + "' domain='Real' value='1.2'/>\n" +
            "      </distr>\n" +
            "    </distribution>\n" +
            "  </distribution>\n" +
            "  <tree id='cppInit' spec='calibratedcpp.CalibratedCPPTreeInitialiser' estimate='false' initial='@Tree'>\n" +
            "    <model idref='cpp'/>\n" +
            "    <calibration idref='calib'/>\n" +
            "  </tree>\n" +
            "</beast>\n";

        List<beast.base.core.BEASTInterface> objects = new XMLParser().parseBareFragments(xml, true);
        Map<String, beast.base.core.BEASTInterface> byId = new java.util.HashMap<>();
        for (beast.base.core.BEASTInterface o : objects) byId.put(o.getID(), o);

        Tree tree = (Tree) byId.get("Tree");
        Distribution cpp = (Distribution) byId.get("cpp");
        Distribution calib = (Distribution) byId.get("calib");

        // The initialiser (run during parse via initAndValidate) produced a constraint-satisfying tree.
        assertValidTree(tree);
        assertClade(tree, taxonSet("A", "B"), 1.5, 2.0);
        assertClade(tree, taxonSet("E", "F"), 0.8, 1.2);
        assertTrue(tree.getRoot().getHeight() >= 2.0, "root above the oldest calibration clade");

        double cppLogP = cpp.calculateLogP();
        assertTrue(Double.isFinite(cppLogP), "calibrated-CPP prior should be finite at the start, was " + cppLogP);

        // Each clade's MRCA prior has finite density at the starting tree. (We sum the wrapper's
        // children directly rather than calling CompoundDistribution.calculateLogP(), whose dirty-
        // flag caching only behaves correctly inside an MCMC/State recalculation.)
        double calibLogP = 0.0;
        for (Distribution child : ((beast.base.inference.CompoundDistribution) calib).pDistributions.get()) {
            double lp = child.calculateLogP();
            assertTrue(Double.isFinite(lp), "MRCA prior should be finite at the start, was " + lp);
            calibLogP += lp;
        }
        assertTrue(Double.isFinite(calibLogP), "MRCA priors should be finite at the start, was " + calibLogP);
    }

    /**
     * The BEAUti template now emits a bare {@code <init initial="@Tree"/>} with no model or
     * calibration reference — the initialiser discovers the calibrated-CPP prior and the MRCA priors
     * from the tree's outputs, so they stay defined only in the posterior. This parses such an XML
     * (no {@code <model>}/{@code <calibration>} on the init) and confirms discovery yields the same
     * finite starting state. initStateNodes() is called explicitly, as MCMC.run() does after parsing.
     */
    @Test
    public void discoversModelAndCalibrationsFromTree() throws Exception {
        Randomizer.setSeed(7);

        String P = "beast.base.spec.inference.parameter.RealScalarParam";
        String U = "beast.base.spec.inference.distribution.Uniform";
        String M = "beast.base.spec.evolution.tree.MRCAPrior";

        String xml =
            "<beast version='2.7'>\n" +
            "  <taxonset id='fullTaxa' spec='beast.base.evolution.alignment.TaxonSet'>\n" +
            "    <taxon id='A' spec='beast.base.evolution.alignment.Taxon'/>\n" +
            "    <taxon id='B' spec='beast.base.evolution.alignment.Taxon'/>\n" +
            "    <taxon id='C' spec='beast.base.evolution.alignment.Taxon'/>\n" +
            "    <taxon id='D' spec='beast.base.evolution.alignment.Taxon'/>\n" +
            "    <taxon id='E' spec='beast.base.evolution.alignment.Taxon'/>\n" +
            "    <taxon id='F' spec='beast.base.evolution.alignment.Taxon'/>\n" +
            "  </taxonset>\n" +
            "  <tree id='Tree' spec='beast.base.evolution.tree.Tree' taxonset='@fullTaxa'/>\n" +
            "  <taxonset id='ab' spec='beast.base.evolution.alignment.TaxonSet'>\n" +
            "    <taxon idref='A'/><taxon idref='B'/></taxonset>\n" +
            "  <taxonset id='ef' spec='beast.base.evolution.alignment.TaxonSet'>\n" +
            "    <taxon idref='E'/><taxon idref='F'/></taxonset>\n" +
            "  <distribution id='cpp' spec='calibratedcpp.CalibratedBirthDeathModel' conditionOnCalibrations='true' conditionOnRoot='true'>\n" +
            "    <tree idref='Tree'/>\n" +
            "    <birthRate spec='" + P + "' domain='PositiveReal' value='1.0'/>\n" +
            "    <deathRate spec='" + P + "' domain='NonNegativeReal' value='0.2'/>\n" +
            "    <rho spec='" + P + "' domain='UnitInterval' value='0.5'/>\n" +
            "    <calibrations idref='ab'/>\n" +
            "    <calibrations idref='ef'/>\n" +
            "  </distribution>\n" +
            "  <distribution id='mrcaAB' spec='" + M + "' monophyletic='true' tree='@Tree' taxonset='@ab'>\n" +
            "    <distr spec='" + U + "'><lower spec='" + P + "' domain='Real' value='1.5'/>" +
            "      <upper spec='" + P + "' domain='Real' value='2.0'/></distr>\n" +
            "  </distribution>\n" +
            "  <distribution id='mrcaEF' spec='" + M + "' monophyletic='true' tree='@Tree' taxonset='@ef'>\n" +
            "    <distr spec='" + U + "'><lower spec='" + P + "' domain='Real' value='0.8'/>" +
            "      <upper spec='" + P + "' domain='Real' value='1.2'/></distr>\n" +
            "  </distribution>\n" +
            "  <tree id='cppInit' spec='calibratedcpp.CalibratedCPPTreeInitialiser' estimate='false' initial='@Tree'/>\n" +
            "</beast>\n";

        List<beast.base.core.BEASTInterface> objects = new XMLParser().parseBareFragments(xml, true);
        Map<String, beast.base.core.BEASTInterface> byId = new java.util.HashMap<>();
        for (beast.base.core.BEASTInterface o : objects) byId.put(o.getID(), o);

        Tree tree = (Tree) byId.get("Tree");
        // Mimic MCMC.run(): call the initialiser after the whole graph is built and outputs are wired.
        ((beast.base.inference.StateNodeInitialiser) byId.get("cppInit")).initStateNodes();

        assertValidTree(tree);
        assertClade(tree, taxonSet("A", "B"), 1.5, 2.0);
        assertClade(tree, taxonSet("E", "F"), 0.8, 1.2);

        assertTrue(Double.isFinite(((Distribution) byId.get("cpp")).calculateLogP()),
                "discovered calibrated-CPP prior should be finite");
        assertTrue(Double.isFinite(((Distribution) byId.get("mrcaAB")).calculateLogP()), "mrcaAB finite");
        assertTrue(Double.isFinite(((Distribution) byId.get("mrcaEF")).calculateLogP()), "mrcaEF finite");
    }

    /**
     * Regression for BEAUti-generated XMLs where a taxon shared by several clades is emitted as
     * duplicate suffixed taxa ({@code A1}, {@code B1}, ...) absent from the tree. The initialiser
     * must resolve these back to the real taxa, place every real taxon exactly once, and keep the
     * resolved clade monophyletic — never crash or drop taxa (which previously threw an
     * ArrayIndexOutOfBounds during Tree.initArrays).
     */
    @Test
    public void duplicateSuffixedCladeTaxaAreResolved() {
        Randomizer.setSeed(5);

        TaxonSet full = taxonSet(TAXA); // A..F
        TaxonSet abDup = taxonSet("A1", "B1"); // phantom duplicates of A and B

        Tree tree = initialTree();
        CalibratedBirthDeathModel model = new CalibratedBirthDeathModel();
        model.initByName("tree", tree,
                "conditionOnRoot", true,
                "birthRate", new RealScalarParam<>(2.0, Real.INSTANCE),
                "deathRate", new RealScalarParam<>(1.0, Real.INSTANCE),
                "rho", new RealScalarParam<>(0.5, Real.INSTANCE),
                "calibrations", new ArrayList<>(List.of(abDup)));

        CalibratedCPPTreeInitialiser init = new CalibratedCPPTreeInitialiser();
        init.initByName("initial", tree, "taxonset", full, "model", model);

        for (int rep = 0; rep < 100; rep++) {
            init.initStateNodes();
            assertValidTree(tree); // exactly 6 taxa, 11 nodes, strictly increasing heights
            Set<String> leaves = new HashSet<>();
            for (Node l : tree.getExternalNodes()) leaves.add(l.getID());
            assertEquals(Set.of("A", "B", "C", "D", "E", "F"), leaves, "all real taxa present, no phantoms");
            assertNotNull(mrca(tree, Set.of("A", "B")), "{A1,B1} clade resolved to a monophyletic {A,B}");
        }
    }

    /**
     * The age-dependent model builds its CDF spline lazily in updateModel(), which the MCMC first
     * calls during likelihood evaluation — after initialisation. The initialiser must prime the law
     * and, failing that, fall back to uniform draws rather than NPE on the null spline.
     */
    @Test
    public void ageDependentModelDoesNotNpeOnLazySpline() {
        Randomizer.setSeed(11);

        TaxonSet full = taxonSet(TAXA);
        TaxonSet ab = taxonSet("A", "B");

        Tree tree = initialTree();
        Exponential lifetime = new Exponential();
        lifetime.initByName("mean", new RealScalarParam<>(1.0, PositiveReal.INSTANCE));

        CalibratedAgeDependentExtinctionModel model = new CalibratedAgeDependentExtinctionModel();
        model.initByName("tree", tree,
                "conditionOnRoot", true,
                "birthRate", new RealScalarParam<>(1.5, PositiveReal.INSTANCE),
                "rho", new RealScalarParam<>(0.5, UnitInterval.INSTANCE),
                "lifetimeDistribution", lifetime,
                "calibrations", new ArrayList<>(List.of(ab)));

        CalibratedCPPTreeInitialiser init = new CalibratedCPPTreeInitialiser();
        init.initByName("initial", tree, "taxonset", full, "model", model);

        for (int rep = 0; rep < 20; rep++) {
            init.initStateNodes(); // must not throw
            assertValidTree(tree);
            assertNotNull(mrca(tree, Set.of("A", "B")), "{A,B} clade monophyletic");
        }
    }

    /**
     * With no model wired (as during BEAUti's createSubNet, before the model's rates are connected)
     * the initialiser must be a harmless no-op — not throw, and leave the target tree untouched.
     */
    /**
     * With a (fitted) CalibrationPrior, clade MRCA ages are drawn from its hierarchy — the top clade
     * from its LogNormal, and an overlapping nested clade from Beta(α,β) on the child/parent ratio.
     * Confirms the resulting starting tree is nested (child younger than parent) and gives the
     * CalibrationPrior a finite starting density.
     */
    @Test
    public void calibrationPriorHierarchicalSampling() {
        Randomizer.setSeed(55);

        TaxonSet full = taxonSet(TAXA);
        TaxonSet ab = taxonSet("A", "B");
        TaxonSet abc = taxonSet("A", "B", "C");

        Tree tree = initialTree();
        CalibratedBirthDeathModel model = new CalibratedBirthDeathModel();
        model.initByName("tree", tree, "conditionOnRoot", true,
                "birthRate", new RealScalarParam<>(1.0, Real.INSTANCE),
                "deathRate", new RealScalarParam<>(0.5, Real.INSTANCE),
                "rho", new RealScalarParam<>(0.5, UnitInterval.INSTANCE),
                "calibrations", new ArrayList<>(List.of(ab, abc)));

        // Overlapping nested intervals ({A,B}⊂{A,B,C}, [3,6] overlaps [4,8]) → Beta on the ratio.
        CalibrationPrior cp = new CalibrationPrior();
        cp.initByName("tree", tree, "calibration", new ArrayList<>(List.of(
                cladePrior(ab, 3.0, 6.0), cladePrior(abc, 4.0, 8.0))));

        CalibratedCPPTreeInitialiser init = new CalibratedCPPTreeInitialiser();
        init.initByName("initial", tree, "taxonset", full, "model", model, "calibration", cp);

        for (int rep = 0; rep < 200; rep++) {
            init.initStateNodes();
            assertValidTree(tree);
            double abAge = mrca(tree, Set.of("A", "B")).getHeight();
            double abcAge = mrca(tree, Set.of("A", "B", "C")).getHeight();
            assertTrue(abAge > 0.0, "positive clade age");
            assertTrue(abAge < abcAge, "nested {A,B} must be younger than {A,B,C}");
            assertTrue(Double.isFinite(cp.calculateLogP()), "CalibrationPrior finite at the start");
        }
    }

    @Test
    public void nullModelIsNoOp() {
        TaxonSet full = taxonSet(TAXA);
        Tree tree = initialTree();

        CalibratedCPPTreeInitialiser init = new CalibratedCPPTreeInitialiser();
        init.initByName("initial", tree, "taxonset", full); // no 'model'

        init.initStateNodes(); // must not throw
        assertEquals(TAXA.length, tree.getLeafNodeCount(), "target tree left intact");
    }

    @Test
    public void noCalibrationsConditionOnRoot() {
        Randomizer.setSeed(7);

        TaxonSet full = taxonSet(TAXA);
        Tree tree = initialTree();
        CalibratedBirthDeathModel model = new CalibratedBirthDeathModel();
        model.initByName("tree", tree,
                "conditionOnRoot", true,
                "birthRate", new RealScalarParam<>(2.0, Real.INSTANCE),
                "deathRate", new RealScalarParam<>(1.0, Real.INSTANCE),
                "rho", new RealScalarParam<>(0.5, Real.INSTANCE));

        CalibratedCPPTreeInitialiser init = new CalibratedCPPTreeInitialiser();
        init.initByName("initial", tree, "taxonset", full, "model", model,
                "rootHeight", 5.0);

        init.initStateNodes();
        assertValidTree(tree);
        assertEquals(5.0, tree.getRoot().getHeight(), 1e-9, "root height should honour rootHeight input");
    }
}
