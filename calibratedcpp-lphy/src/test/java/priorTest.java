import calibratedcpp.lphy.prior.Calibration;
import calibratedcpp.lphy.prior.ConditionedMRCAPrior;
import lphy.core.model.Value;
import org.junit.jupiter.api.Test;

import static calibratedcpp.lphy.prior.ConditionedMRCAPrior.mapBetaNodes;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class priorTest {

    private static Calibration[] buildCalibrations(String[][] taxa, Double[] upper, Double[] lower) {
        Calibration[] cals = new Calibration[taxa.length];
        for (int i = 0; i < taxa.length; i++) {
            cals[i] = new Calibration(taxa[i], upper[i], lower[i]);
        }
        return cals;
    }

    @Test
    void testMapBetaArray_nonOverlap() {
        Double[] upperBounds = new Double[]{7.0, 4.5, 3.0, 3.0, 4.0};
        Double[] lowerBounds = new Double[]{5.0, 4.0, 2.0, 1.0, 2.0};
        int[] parent = new int[]{-1, 0, 1, 1, 0};
        Calibration[] cals = buildCalibrations(
                new String[][]{{"a"}, {"b"}, {"c"}, {"d"}, {"e"}},
                upperBounds, lowerBounds);

        boolean[] expect = new boolean[]{false, false, false, false, false};
        boolean[] actual = mapBetaNodes(cals, parent);

        for (int i = 0; i < actual.length; i++) assertEquals(expect[i], actual[i]);
    }

    @Test
    void testMapBetaArray_withOverlap() {
        Double[] upperBounds = new Double[]{7.0, 4.5, 3.0, 4.1, 4.0};
        Double[] lowerBounds = new Double[]{5.0, 4.0, 2.0, 1.0, 2.0};
        int[] parent = new int[]{-1, 0, 1, 1, 0};
        Calibration[] cals = buildCalibrations(
                new String[][]{{"a"}, {"b"}, {"c"}, {"d"}, {"e"}},
                upperBounds, lowerBounds);

        boolean[] expect = new boolean[]{false, false, false, true, false};
        boolean[] actual = mapBetaNodes(cals, parent);

        for (int i = 0; i < actual.length; i++) assertEquals(expect[i], actual[i]);
    }

    @Test
    void testExample() {
        // without root calibration
        Double[] upperBounds = new Double[]{5.0, 5.0, 5.0, 5.0, 6.0};
        Double[] lowerBounds = new Double[]{4.0, 4.0, 4.8, 4.8, 4.9};
        int[] parent = new int[]{2, 3, -1, 4, -1};
        Calibration[] cals = buildCalibrations(
                new String[][]{{"a"}, {"b"}, {"c"}, {"d"}, {"e"}},
                upperBounds, lowerBounds);

        boolean[] expect = new boolean[]{true, true, false, true, false};
        boolean[] actual = mapBetaNodes(cals, parent);

        for (int i = 0; i < actual.length; i++) assertEquals(expect[i], actual[i]);
    }

    @Test
    void testExample2() {
        // with root calibration
        Double[] upperBounds = new Double[]{6.0, 5.0, 5.0};
        Double[] lowerBounds = new Double[]{4.9, 4.0, 4.8};
        int[] parent = new int[]{-1, 2, 0};
        Calibration[] cals = buildCalibrations(
                new String[][]{{"a"}, {"b"}, {"c"}},
                upperBounds, lowerBounds);

        boolean[] expect = new boolean[]{false, true, true};
        boolean[] actual = mapBetaNodes(cals, parent);

        for (int i = 0; i < actual.length; i++) assertEquals(expect[i], actual[i]);
    }

    @Test
    void outputArray() {
        Calibration[] calibrationSpecs = new Calibration[]{
                new Calibration(new String[]{"1", "2", "3", "4", "5"}, 6.0, 5.0),
                new Calibration(new String[]{"1", "2", "3", "4"}, 4.0, 2.4),
                new Calibration(new String[]{"1", "2"}, 3.0, 2.0),
                new Calibration(new String[]{"3", "4"}, 2.5, 1.0)
        };
        ConditionedMRCAPrior conditionedMRCAPrior = new ConditionedMRCAPrior(
                new Value<>("", calibrationSpecs), null);
        Calibration[] observed = conditionedMRCAPrior.sample().value().getCalibrationArray();
        assertEquals(observed.length, calibrationSpecs.length);
        for (int i = 0; i < observed.length; i++) {
            if (i == 0) {
                assertEquals(5, observed[i].getTaxa().length);
            } else {
                assertEquals(observed[i].getTaxa().length, calibrationSpecs[i].getTaxa().length);
                assertTrue(observed[0].getAge() > observed[1].getAge());
                assertTrue(observed[1].getAge() > observed[3].getAge());
            }
        }
    }

    @Test
    void outputArrayNoRoot() {
        Calibration[] calibrationSpecs = new Calibration[]{
                new Calibration(new String[]{"1", "2", "3", "4"}, 4.0, 2.4),
                new Calibration(new String[]{"1", "2"}, 3.0, 2.0),
                new Calibration(new String[]{"3", "4"}, 2.5, 1.0)
        };
        ConditionedMRCAPrior conditionedMRCAPrior = new ConditionedMRCAPrior(
                new Value<>("", calibrationSpecs), null);
        Calibration[] observed = conditionedMRCAPrior.sample().value().getCalibrationArray();
        assertEquals(observed.length, calibrationSpecs.length);
        for (int i = 0; i < observed.length; i++) {
            assertEquals(observed[i].getTaxa().length, calibrationSpecs[i].getTaxa().length);
            assertTrue(observed[0].getAge() > observed[1].getAge());
            assertTrue(observed[0].getAge() > observed[2].getAge());
        }
    }

    @Test
    void inferenceTest() {
        Calibration[] calibrationSpecs = new Calibration[]{
                new Calibration(new String[]{"1", "2", "3"}, 1.5, 1.0),
                new Calibration(new String[]{"1", "2", "3", "4"}, 1.9, 1.5),
                new Calibration(new String[]{"6", "7"}, 1.5, 1.1)
        };
        ConditionedMRCAPrior conditionedMRCAPrior = new ConditionedMRCAPrior(
                new Value<>("", calibrationSpecs), null);
        Calibration[] observed = conditionedMRCAPrior.sample().value().getCalibrationArray();
        for (int i = 0; i < observed.length; i++) {
            assertEquals(observed[i].getTaxa().length, calibrationSpecs[i].getTaxa().length);
            assertTrue(observed[0].getAge() < observed[1].getAge());
        }
    }

    @Test
    void overlappingTest() {
        Double[] upperBounds = new Double[]{6.0, 5.0, 5.0};
        Double[] lowerBounds = new Double[]{4.9, 4.9, 4.0};
        int[] parent = new int[]{-1, 0, 1};
        Calibration[] cals = buildCalibrations(
                new String[][]{{"a"}, {"b"}, {"c"}},
                upperBounds, lowerBounds);

        boolean[] expect = new boolean[]{false, true, true};
        boolean[] actual = mapBetaNodes(cals, parent);

        for (int i = 0; i < actual.length; i++) assertEquals(expect[i], actual[i]);

        Calibration[] calibrationSpecs = new Calibration[]{
                new Calibration(new String[]{"1", "2", "3", "4"}, 6.0, 4.9),
                new Calibration(new String[]{"1", "2", "3"}, 5.0, 4.9),
                new Calibration(new String[]{"1", "2"}, 5.0, 4.0)
        };
        ConditionedMRCAPrior conditionedMRCAPrior = new ConditionedMRCAPrior(
                new Value<>("", calibrationSpecs), null);
        Calibration[] observed = conditionedMRCAPrior.sample().value().getCalibrationArray();
        for (int i = 0; i < observed.length; i++) {
            assertEquals(observed[i].getTaxa().length, calibrationSpecs[i].getTaxa().length);
            assertTrue(observed[0].getAge() > observed[1].getAge());
            assertTrue(observed[1].getAge() > observed[2].getAge());
        }
    }
}
