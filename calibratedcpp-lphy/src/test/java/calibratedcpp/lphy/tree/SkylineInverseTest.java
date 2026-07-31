package calibratedcpp.lphy.tree;

import lphy.core.model.Value;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

/**
 * Validates the skyline's closed-form {@code invertCDF} against {@code cdf} (the directly-ported BEAST
 * log-CDF) by round-tripping across a three-interval process: cdf(invertCDF(p)) must equal p, and the
 * inverse must be monotone and span the interval boundaries. Same-package so the protected hooks are
 * reachable.
 */
public class SkylineInverseTest {

    private static Value<Double[]> arr(double... v) {
        Double[] a = new Double[v.length];
        for (int i = 0; i < v.length; i++) a[i] = v[i];
        return new Value<>("", a);
    }

    @Test
    public void analyticInverseRoundTripsAcrossIntervals() {
        // three intervals with change times at ages 2 and 5; mixed super/sub-critical segments
        CalibratedBirthDeathSkylineTree sky = new CalibratedBirthDeathSkylineTree(
                arr(2.0, 0.8, 1.5), arr(0.5, 1.2, 0.9), null, null, null, arr(2.0, 5.0),
                new Value<>("", 1.0), new Value<>("", 5), null, null, null, new Value<>("", 8.0));
        sky.resolveRates();

        double prev = -1.0;
        for (double p = 0.005; p < 0.999; p += 0.005) {
            double t = sky.invertCDF(p);
            assertTrue(t > prev, "invertCDF must be monotone increasing in p (p=" + p + ")");
            prev = t;
            assertEquals(p, sky.cdf(t), 1e-9, "cdf(invertCDF(p)) must equal p at p=" + p + " (t=" + t + ")");
        }
        // the sweep must actually cross both change times, exercising all three interval inverses
        assertTrue(sky.invertCDF(0.005) < 2.0, "small p should land in the first interval");
        assertTrue(sky.invertCDF(0.99) > 5.0, "large p should land in the last interval");
    }
}
