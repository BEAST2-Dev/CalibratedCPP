package calibratedcpp.lphy.prior;

import lphy.core.model.Value;
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.special.Gamma;

public class ConditionedPriorUtils {

    /** Converts an array-literal Value (e.g. {@code [cal1, cal2, ...]}) to a Calibration[]. */
    public static Calibration[] extractCalibrations(Value<?> v) {
        Object raw = v.value();
        if (raw instanceof Calibration[] arr) return arr;
        Object[] objs = (Object[]) raw;
        Calibration[] result = new Calibration[objs.length];
        for (int i = 0; i < objs.length; i++) result[i] = (Calibration) objs[i];
        return result;
    }

    // methods copied from calibratedcpp-beast/calibrationprior/CalibrationPrior.java
    public static double computeLogTargetsMu(double tLo, double tHi, double p){
        // z-score for the coverage probability via Erf.erfInv
        double z = Math.sqrt(2.0) * Erf.erfInv(p);
        double sigma = (Math.log(tHi) - Math.log(tLo)) / (2 * z);
        return Math.log(tLo) + z * sigma;
    }

    public static double computeLogTargetsSigma2(double tLo, double tHi, double p){
        // z-score for the coverage probability via Erf.erfInv
        double z = Math.sqrt(2.0) * Erf.erfInv(p);
        double sigma = (Math.log(tHi) - Math.log(tLo)) / (2 * z);
        return sigma * sigma;
    }

    private static final double MIN_CONCENTRATION_FACTOR = 1.0;

    public static double[] invertLogMomentsToBeta(double m, double v) {
        double mu0 = Math.exp(m);
        mu0 = Math.min(1.0 - 1e-6, Math.max(1e-6, mu0));
        double n0 = Math.max(1e-3, (1.0 - mu0) / (mu0 * Math.max(v, 1e-15)));
        double a = Math.max(1e-6, mu0 * n0);
        double b = Math.max(1e-6, n0 - a);

        for (int iter = 0; iter < 200; iter++) {
            double n = a + b;
            double f1 = Gamma.digamma(a) - Gamma.digamma(n) - m;
            double f2 = Gamma.trigamma(a) - Gamma.trigamma(n) - v;
            if (Math.abs(f1) < 1e-10 && Math.abs(f2) < 1e-10) break;

            double triA = Gamma.trigamma(a);
            double triN = Gamma.trigamma(n);
            double hA = Math.max(a * 1e-5, 1e-9);
            double hN = Math.max(n * 1e-5, 1e-9);
            double tetA = (Gamma.trigamma(a + hA) - Gamma.trigamma(a - hA)) / (2 * hA);
            double tetN = (Gamma.trigamma(n + hN) - Gamma.trigamma(n - hN)) / (2 * hN);

            double j11 = triA - triN;
            double j12 = -triN;
            double j21 = tetA - tetN;
            double j22 = -tetN;
            double det = j11 * j22 - j12 * j21;
            if (Math.abs(det) < 1e-20) break;

            double da = (-j22 * f1 + j12 * f2) / det;
            double db = ( j21 * f1 - j11 * f2) / det;

            // Armijo line search: halve step until ||F||^2 decreases
            double res0 = f1 * f1 + f2 * f2;
            double step = 1.0;
            for (int ls = 0; ls < 50; ls++) {
                double an = Math.max(a + step * da, 1e-6);
                double bn = Math.max(b + step * db, 1e-6);
                double nn = an + bn;
                double g1 = Gamma.digamma(an) - Gamma.digamma(nn) - m;
                double g2 = Gamma.trigamma(an) - Gamma.trigamma(nn) - v;
                if (g1 * g1 + g2 * g2 < res0) break;
                step *= 0.5;
            }
            a = Math.max(a + step * da, 1e-6);
            b = Math.max(b + step * db, 1e-6);
        }

        double mu = a / (a + b);
        double minN = MIN_CONCENTRATION_FACTOR / (mu * (1.0 - mu));
        if (a + b < minN) {
            a = minN * mu;
            b = minN * (1.0 - mu);
        }

        return new double[]{a, b};
    }
}
