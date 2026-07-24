package calibratedcpp;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveInt;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.inference.distribution.Exponential;
import beast.base.spec.inference.parameter.IntScalarParam;
import beast.base.spec.inference.parameter.RealScalarParam;
import calibratedcpp.distribution.Erlang;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 * Validates CalibratedAgeDependentBirthDeathModel against CalibratedBirthDeathModel.
 *
 * <p>Erlang(n=1, scale=1/mu) is identical to Exponential(rate=mu), which in turn
 * matches a constant-rate birth-death process with death rate mu. Both the closed-form
 * Erlang path and the numerical VIDE path are tested against the reference BD model.</p>
 */
public class CalibratedAgeDependentBirthDeathModelTest {

    // Tree with root age 2 (tips adjusted to height 0)
    private static Tree smallTree() {
        Tree tree = new TreeParser();
        tree.initByName("newick", "((A:1,B:1):1,C:2);",
                "IsLabelledNewick", true,
                "adjustTipHeights", true);
        return tree;
    }

    private static CalibratedBirthDeathModel bdModel(
            Tree tree, double lambda, double mu, double rho, double origin) {
        CalibratedBirthDeathModel bd = new CalibratedBirthDeathModel();
        bd.initByName("tree", tree,
                "origin",    new RealScalarParam<>(origin, PositiveReal.INSTANCE),
                "birthRate", new RealScalarParam<>(lambda, PositiveReal.INSTANCE),
                "deathRate", new RealScalarParam<>(mu, NonNegativeReal.INSTANCE),
                "rho",       new RealScalarParam<>(rho, UnitInterval.INSTANCE));
        return bd;
    }

    /**
     * Erlang(alpha=1, beta=1/mu) closed-form path.
     * beta is BEAST3 Gamma's scale parameter, so theta = 1/beta = mu = death rate.
     */
    private static CalibratedAgeDependentExtinctionModel erlangModel(
            Tree tree, double lambda, double mu, double rho, double origin) {
        Erlang erlang = new Erlang();
        erlang.initByName("shape", new IntScalarParam<>(1, PositiveInt.INSTANCE),
                         "scale",  new RealScalarParam<>(1.0 / mu, PositiveReal.INSTANCE));
        CalibratedAgeDependentExtinctionModel model = new CalibratedAgeDependentExtinctionModel();
        model.initByName("tree", tree,
                "origin",              new RealScalarParam<>(origin, PositiveReal.INSTANCE),
                "birthRate",           new RealScalarParam<>(lambda, PositiveReal.INSTANCE),
                "rho",                 new RealScalarParam<>(rho, UnitInterval.INSTANCE),
                "lifetimeDistribution", erlang);
        return model;
    }

    /**
     * Exponential(mean=1/mu) VIDE numerical path.
     * Uses a large origin so the solver covers the full test range [0, 20].
     */
    private static CalibratedAgeDependentExtinctionModel videModel(
            Tree tree, double lambda, double mu, double rho, double origin) {
        Exponential expDist = new Exponential();
        expDist.initByName("mean", new RealScalarParam<>(1.0 / mu, PositiveReal.INSTANCE));
        CalibratedAgeDependentExtinctionModel model = new CalibratedAgeDependentExtinctionModel();
        model.initByName("tree", tree,
                "origin",              new RealScalarParam<>(origin, PositiveReal.INSTANCE),
                "birthRate",           new RealScalarParam<>(lambda, PositiveReal.INSTANCE),
                "rho",                 new RealScalarParam<>(rho, UnitInterval.INSTANCE),
                "lifetimeDistribution", expDist,
                "gridSize",            20000);
        return model;
    }

    // -------------------------------------------------------------------------
    // Erlang (closed-form) vs CalibratedBirthDeathModel
    // -------------------------------------------------------------------------

    @Test
    public void erlangShape1MatchesBDLogDensityTest() {
        double[] times = {0.5, 1.0, 2.0, 5.0, 10.0};

        // Supercritical: lambda > mu
        Tree tree = smallTree();
        CalibratedBirthDeathModel bdSuper = bdModel(tree, 2.0, 1.0, 0.1, 20.0);
        CalibratedAgeDependentExtinctionModel erlangSuper = erlangModel(tree, 2.0, 1.0, 0.1, 20.0);
        erlangSuper.calculateTreeLogLikelihood(tree); // initialise roots and residues

        for (double t : times) {
            assertEquals(bdSuper.calculateLogNodeAgeDensity(t),
                         erlangSuper.calculateLogNodeAgeDensity(t),
                         1e-6, "Erlang(1) density ≠ BD (supercritical) at t=" + t);
        }

        // Subcritical: lambda < mu
        Tree tree2 = smallTree();
        CalibratedBirthDeathModel bdSub = bdModel(tree2, 1.0, 2.0, 0.1, 20.0);
        CalibratedAgeDependentExtinctionModel erlangSub = erlangModel(tree2, 1.0, 2.0, 0.1,20.0);
        erlangSub.calculateTreeLogLikelihood(tree2);

        for (double t : times) {
            assertEquals(bdSub.calculateLogNodeAgeDensity(t),
                         erlangSub.calculateLogNodeAgeDensity(t),
                         1e-6, "Erlang(1) density ≠ BD (subcritical) at t=" + t);
        }
    }

    @Test
    public void erlangShape1MatchesBDLogCDFTest() {
        double[] times = {0.5, 1.0, 2.0, 5.0, 10.0};

        Tree tree = smallTree();
        CalibratedBirthDeathModel bdSuper = bdModel(tree, 2.0, 1.0, 0.1, 20.0);
        CalibratedAgeDependentExtinctionModel erlangSuper = erlangModel(tree, 2.0, 1.0, 0.1,20.0);
        erlangSuper.calculateTreeLogLikelihood(tree);

        for (double t : times) {
            assertEquals(bdSuper.calculateLogNodeAgeCDF(t),
                         erlangSuper.calculateLogNodeAgeCDF(t),
                         1e-6, "Erlang(1) CDF ≠ BD (supercritical) at t=" + t);
        }

        Tree tree2 = smallTree();
        CalibratedBirthDeathModel bdSub = bdModel(tree2, 1.0, 2.0, 0.1, 20.0);
        CalibratedAgeDependentExtinctionModel erlangSub = erlangModel(tree2, 1.0, 2.0, 0.1,20.0);
        erlangSub.calculateTreeLogLikelihood(tree2);

        for (double t : times) {
            assertEquals(bdSub.calculateLogNodeAgeCDF(t),
                         erlangSub.calculateLogNodeAgeCDF(t),
                         1e-6, "Erlang(1) CDF ≠ BD (subcritical) at t=" + t);
        }
    }

    @Test
    public void erlangShape1MatchesBDTreeLogLikelihoodTest() {
        // Supercritical
        Tree tree = smallTree();
        CalibratedBirthDeathModel bdSuper = bdModel(tree, 2.0, 1.0, 0.1,20.0);
        CalibratedAgeDependentExtinctionModel erlangSuper = erlangModel(tree, 2.0, 1.0, 0.1,20.0);

        assertEquals(bdSuper.calculateTreeLogLikelihood(tree),
                     erlangSuper.calculateTreeLogLikelihood(tree),
                     1e-6, "Erlang(1) tree log likelihood ≠ BD (supercritical)");

        // Subcritical
        Tree tree2 = smallTree();
        CalibratedBirthDeathModel bdSub = bdModel(tree2, 1.0, 2.0, 0.1, 20.0);
        CalibratedAgeDependentExtinctionModel erlangSub = erlangModel(tree2, 1.0, 2.0, 0.1,20.0);

        assertEquals(bdSub.calculateTreeLogLikelihood(tree2),
                     erlangSub.calculateTreeLogLikelihood(tree2),
                     1e-6, "Erlang(1) tree log likelihood ≠ BD (subcritical)");
    }

    // -------------------------------------------------------------------------
    // VIDE numerical solver (Exponential distribution) vs CalibratedBirthDeathModel
    // -------------------------------------------------------------------------

    @Test
    public void videExponentialMatchesBDLogDensityTest() {
        // VIDE tolerance is looser: implicit trapezoidal with h=20/2000=0.01 gives O(h^2)~1e-4 error
        double[] times = {0.5, 1.0, 2.0, 5.0, 10.0};

        Tree tree = smallTree();
        CalibratedBirthDeathModel bdSuper = bdModel(tree, 2.0, 1.0, 0.1, 20.0);
        CalibratedAgeDependentExtinctionModel videSuper = videModel(tree, 2.0, 1.0, 0.1,20.0);
        videSuper.calculateTreeLogLikelihood(tree); // solves VIDE on [0, 20]

        for (double t : times) {
            assertEquals(bdSuper.calculateLogNodeAgeDensity(t),
                         videSuper.calculateLogNodeAgeDensity(t),
                         1e-6, "VIDE (Exponential) density ≠ BD (supercritical) at t=" + t);
        }

        Tree tree2 = smallTree();
        CalibratedBirthDeathModel bdSub = bdModel(tree2, 1.0, 2.0, 0.1, 20.0);
        CalibratedAgeDependentExtinctionModel videSub = videModel(tree2, 1.0, 2.0, 0.1, 20.0);
        videSub.calculateTreeLogLikelihood(tree2);

        for (double t : times) {
            assertEquals(bdSub.calculateLogNodeAgeDensity(t),
                         videSub.calculateLogNodeAgeDensity(t),
                         1e-6, "VIDE (Exponential) density ≠ BD (subcritical) at t=" + t);
        }
    }

    @Test
    public void videExponentialMatchesBDLogCDFTest() {
        double[] times = {0.5, 1.0, 2.0, 5.0, 10.0};

        Tree tree = smallTree();
        CalibratedBirthDeathModel bdSuper = bdModel(tree, 2.0, 1.0, 0.1, 20.0);
        CalibratedAgeDependentExtinctionModel videSuper = videModel(tree, 2.0, 1.0, 0.1, 20.0);
        videSuper.calculateTreeLogLikelihood(tree);

        for (double t : times) {
            assertEquals(bdSuper.calculateLogNodeAgeCDF(t),
                         videSuper.calculateLogNodeAgeCDF(t),
                         1e-6, "VIDE (Exponential) CDF ≠ BD (supercritical) at t=" + t);
        }

        Tree tree2 = smallTree();
        CalibratedBirthDeathModel bdSub = bdModel(tree2, 1.0, 2.0, 0.1, 20.0);
        CalibratedAgeDependentExtinctionModel videSub = videModel(tree2, 1.0, 2.0, 0.1, 20.0);
        videSub.calculateTreeLogLikelihood(tree2);

        for (double t : times) {
            assertEquals(bdSub.calculateLogNodeAgeCDF(t),
                         videSub.calculateLogNodeAgeCDF(t),
                         1e-6, "VIDE (Exponential) CDF ≠ BD (subcritical) at t=" + t);
        }
    }

    @Test
    public void videExponentialMatchesBDTreeLogLikelihoodTest() {
        // Both use conditionOnRoot=true; VIDE solver runs on [0, rootAge=2]
        // All internal node ages fall within this range.
        Tree tree = smallTree();
        CalibratedBirthDeathModel bdSuper = bdModel(tree, 2.0, 1.0, 0.1, 20.0);

        Exponential expDist = new Exponential();
        expDist.initByName("mean", new RealScalarParam<>(1.0, PositiveReal.INSTANCE));
        CalibratedAgeDependentExtinctionModel videSuper = new CalibratedAgeDependentExtinctionModel();
        videSuper.initByName("tree", tree,
                "origin",              new RealScalarParam<>(20.0, PositiveReal.INSTANCE),
                "birthRate",           new RealScalarParam<>(2.0, PositiveReal.INSTANCE),
                "rho",                 new RealScalarParam<>(0.1, UnitInterval.INSTANCE),
                "lifetimeDistribution", expDist,
                "gridSize",            2000);

        assertEquals(bdSuper.calculateTreeLogLikelihood(tree),
                     videSuper.calculateTreeLogLikelihood(tree),
                     1e-6, "VIDE (Exponential) tree log likelihood ≠ BD (supercritical)");

        // Subcritical
        Tree tree2 = smallTree();
        CalibratedBirthDeathModel bdSub = bdModel(tree2, 1.0, 2.0, 0.1, 20.0);

        Exponential expDist2 = new Exponential();
        expDist2.initByName("mean", new RealScalarParam<>(0.5, PositiveReal.INSTANCE)); // mean = 1/mu = 1/2
        CalibratedAgeDependentExtinctionModel videSub = new CalibratedAgeDependentExtinctionModel();
        videSub.initByName("tree", tree2,
                "origin",              new RealScalarParam<>(20.0, PositiveReal.INSTANCE),
                "birthRate",           new RealScalarParam<>(1.0, PositiveReal.INSTANCE),
                "rho",                 new RealScalarParam<>(0.1, UnitInterval.INSTANCE),
                "lifetimeDistribution", expDist2,
                "gridSize",            20000);

        assertEquals(bdSub.calculateTreeLogLikelihood(tree2),
                     videSub.calculateTreeLogLikelihood(tree2),
                     1e-6, "VIDE (Exponential) tree log likelihood ≠ BD (subcritical)");
    }
}
