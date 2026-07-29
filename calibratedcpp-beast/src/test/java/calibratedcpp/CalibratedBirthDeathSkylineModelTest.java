package calibratedcpp;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeParser;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.RealScalar;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
//TODO import bdsky.evolution.speciation.BirthDeathSkylineModel;
import static org.junit.jupiter.api.Assertions.*;

public class CalibratedBirthDeathSkylineModelTest {

   private Tree tree;
   private RealScalar<PositiveReal> b;
   private RealScalar<NonNegativeReal> d;
   private RealScalarParam<UnitInterval> rho;
   private RealScalarParam<PositiveReal> origin;

   @BeforeEach
   public void setUp() {
      tree = new TreeParser();
      tree.initByName("newick", "((A:1,B:1):1,C:2);",
              "IsLabelledNewick", true,
              "adjustTipHeights", true);

      b = new RealScalarParam<>(2.0, PositiveReal.INSTANCE);
      d = new RealScalarParam<>(1.0, NonNegativeReal.INSTANCE);
      rho = new RealScalarParam<>(1.0, UnitInterval.INSTANCE);
      origin = new RealScalarParam<>(5.0, PositiveReal.INSTANCE);
   }

   private SkylineParameter createSkyline(double... rateValues) {
      RealVectorParam<NonNegativeReal> rates = new RealVectorParam<>(rateValues, NonNegativeReal.INSTANCE);
      SkylineParameter sp = new SkylineParameter();
      sp.initByName("values", rates);
      return sp;
   }

   private SkylineParameter createSkyline(double[] rateValues, double[] timeValues, boolean relative, boolean reverse) {
      RealVectorParam<NonNegativeReal> rates = new RealVectorParam<>(rateValues, NonNegativeReal.INSTANCE);
      SkylineParameter sp = new SkylineParameter();
      if (timeValues != null) {
         RealVectorParam<NonNegativeReal> times = new RealVectorParam<>(timeValues, NonNegativeReal.INSTANCE);
         sp.initByName("values", rates, "changeTimes", times, "timesAreRelative", relative, "timesAreAges", reverse);
      } else {
         sp.initByName("values", rates, "timesAreRelative", relative, "timesAreAges", reverse);
      }
      return sp;
   }

   @Test
   public void testOnlyOneRateThrows() {
      CalibratedBirthDeathSkylineModel model = new CalibratedBirthDeathSkylineModel();

       assertThrows(RuntimeException.class, () -> model.initByName(
               "birthRate", createSkyline(b.get()),
               "rho", rho,
               "tree", tree,
               "origin", origin,
               "conditionOnRoot", true
       ), "Should throw because exactly TWO rates must be specified.");
   }

   @Test
   public void testThreeRatesThrows() {
      CalibratedBirthDeathSkylineModel model = new CalibratedBirthDeathSkylineModel();

      assertThrows(RuntimeException.class, () -> model.initByName(
              "birthRate", createSkyline(b.get()),
              "deathRate", createSkyline(d.get()),
              "diversificationRate", createSkyline(0.5),
              "rho", rho,
              "tree", tree,
              "origin", origin,
              "conditionOnRoot", true
      ), "Should throw because mutual exclusion check failed (3 rates provided).");
   }

   @Test
   public void testRelativeOutOfBoundsThrows() {
      assertThrows(RuntimeException.class,
              () -> createSkyline(new double[]{2.0, 1.0}, new double[]{1.5}, true, false),
              "Should throw because relative time 1.5 is outside [0, 1].");
   }

   @Test
   public void testLikelihoodMatchesConstantBD() {
      CalibratedBirthDeathSkylineModel bdsky = new CalibratedBirthDeathSkylineModel();
      bdsky.initByName(
              "birthRate", createSkyline(b.get()),
              "deathRate", createSkyline(d.get()),
              "rho", rho,
              "tree", tree,
              "origin", origin,
              "conditionOnRoot", true
      );

      CalibratedBirthDeathModel constantBD = new CalibratedBirthDeathModel();
      constantBD.initByName("birthRate", b, "deathRate", d, "rho", rho,
              "tree", tree, "origin", origin, "conditionOnRoot", true);

      double logL_BDSKY = bdsky.calculateTreeLogLikelihood(tree);
      double logL_BD = constantBD.calculateTreeLogLikelihood(tree);

      assertEquals(logL_BD, logL_BDSKY, 1e-10,
              "BDSKY likelihood should match constant Birth-Death when no intervals are defined.");
   }

   @Test
   public void calculateLogNodeAgeDensityTest() {
      CalibratedBirthDeathSkylineModel model = new CalibratedBirthDeathSkylineModel();
      model.initByName(
              "birthRate", createSkyline(b.get()),
              "deathRate", createSkyline(d.get()),
              "rho", rho,
              "tree", tree,
              "origin", origin,
              "conditionOnRoot", true
      );

      model.calculateTreeLogLikelihood(tree);
      double density = model.calculateLogNodeAgeDensity(0.5);
      assertTrue(Double.isFinite(density), "Density should be a finite log value.");
   }

   @Test
   public void calculateLogNodeAgeCDFTest() {
      CalibratedBirthDeathSkylineModel model = new CalibratedBirthDeathSkylineModel();
      model.initByName(
              "birthRate", createSkyline(b.get()),
              "deathRate", createSkyline(d.get()),
              "rho", rho,
              "tree", tree,
              "origin", origin,
              "conditionOnRoot", true
      );

      model.calculateTreeLogLikelihood(tree);

      double cdfZero = model.calculateLogNodeAgeCDF(0.0);
      assertEquals(Double.NEGATIVE_INFINITY, cdfZero, 1e-10, "log(CDF) at t=0 should be -Infinity.");

      double cdfRoot = model.calculateLogNodeAgeCDF(tree.getRoot().getHeight());
      assertTrue(cdfRoot <= 0.0, "log(CDF) should never be positive.");
   }

//   @Test
//   public void calculateTreeLogLikelihood() {
//      // BEAST2 RealParameter objects for BirthDeathSkylineModel (bdsky package)
//      RealParameter birthRatesBeast2 = new RealParameter("2.0 1.0 3.0");
//      RealParameter deathRatesBeast2 = new RealParameter("1.1 2.0 0.5");
//      RealParameter rhoBeast2 = new RealParameter("0.5");
//      RealParameter birthChangesBeast2 = new RealParameter("0.0 1.0 1.5");
//      RealParameter deathChangesBeast2 = new RealParameter("0.0 0.5 1.25");
//      TreeInterface localTree = new TreeParser("((A:3,B:3):4,(C:6,D:6):1);");
//
//      BirthDeathSkylineModel BDSKY = new BirthDeathSkylineModel();
//      BDSKY.initByName("birthRate", birthRatesBeast2,
//              "deathRate", deathRatesBeast2,
//              "rho", rhoBeast2,
//              "birthRateChangeTimes", birthChangesBeast2,
//              "deathRateChangeTimes", deathChangesBeast2,
//              "samplingRate", new RealParameter("0.0"),
//              "origin", new RealParameter("8.0"),
//              "conditionOnSurvival", true,
//              "reverseTimeArrays", new BooleanParameter("true true true true true"),
//              "tree", localTree);
//
//      // BEAST3 types for CalibratedBirthDeathSkylineModel
//      SkylineParameter birthRateSkyline = new SkylineParameter();
//      SkylineParameter deathRateSkyline = new SkylineParameter();
//      birthRateSkyline.initByName("values", new RealVectorParam<>(new double[]{2.0, 1.0, 3.0}, NonNegativeReal.INSTANCE),
//              "changeTimes", new RealVectorParam<>(new double[]{1.0, 1.5}, NonNegativeReal.INSTANCE),
//              "timesAreAges", true);
//      deathRateSkyline.initByName("values", new RealVectorParam<>(new double[]{1.1, 2.0, 0.5}, NonNegativeReal.INSTANCE),
//              "changeTimes", new RealVectorParam<>(new double[]{0.5, 1.25}, NonNegativeReal.INSTANCE),
//              "timesAreAges", true);
//
//      CalibratedBirthDeathSkylineModel CalibratedBDSKY = new CalibratedBirthDeathSkylineModel();
//      CalibratedBDSKY.initByName("birthRate", birthRateSkyline,
//              "deathRate", deathRateSkyline,
//              "rho", new RealScalarParam<>(0.5, UnitInterval.INSTANCE),
//              "tree", localTree,
//              "origin", new RealScalarParam<>(8.0, PositiveReal.INSTANCE));
//
//      assertEquals(BDSKY.calculateTreeLogLikelihood(localTree) + 3.0 * Math.log(2.0) - Math.log(4.0) - Math.log(3.0) - Math.log(2.0),
//              CalibratedBDSKY.calculateTreeLogLikelihood(localTree), 1e-8,
//              "Likelihood of tree under BDSKY does not match likelihood under Calibrated BDSKY.");
//   }
}
