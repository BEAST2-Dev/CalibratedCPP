import calibratedcpp.lphy.prior.Calibration;
import calibratedcpp.lphy.prior.ConditionedMRCAPrior;
import lphy.core.model.Value;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Standalone runner: samples ConditionedMRCAPrior N_SIM times using the same
 * 10-node example tree as calibrationprior_simulation.r, then writes the node
 * ages (Wsim) and per-node coverage to a CSV file.
 *
 * Run this first, then run compare_distributions.R to visualise the overlap.
 *
 * Usage: run main() directly from IDE, or via maven exec:java.
 * Output: validation/calibrationprior/java_wsim.csv
 */
public class ConditionedMRCAPriorSampler {

    // Same parameters as calibrationprior_simulation.r
    private static final Double[] T_LO = {10.0, 9.4, 4.8, 4.0, 8.0, 6.8, 1.8, 2.0, 1.7, 1.6};
    private static final Double[] T_HI = {10.5, 9.6, 5.0, 5.0, 9.5, 8.0, 2.5, 3.2, 2.5, 2.5};
    private static final double P_COV = 0.90;
    private static final int N_SIM = 50000;
    private static final int N = T_LO.length;

    private static final String OUTPUT_CSV = "./calibratedcpp-lphy/src/test/java/javaResults.csv";

    public static void main(String[] args) throws IOException {
        Calibration[] specs = buildCalibrationSpecs();
        ConditionedMRCAPrior prior = new ConditionedMRCAPrior(
                new Value<>("", specs),
                new Value<>("", P_COV)
        );

        double[][] wsim = new double[N_SIM][N];
        for (int s = 0; s < N_SIM; s++) {
            Calibration[] cals = prior.sample().value().getCalibrationArray();
            for (int i = 0; i < N; i++) {
                wsim[s][i] = cals[i].getAge();
            }
        }

        // per-node coverage
        double[] coverage = new double[N];
        for (int i = 0; i < N; i++) {
            int count = 0;
            for (int s = 0; s < N_SIM; s++) {
                if (wsim[s][i] >= T_LO[i] && wsim[s][i] <= T_HI[i]) count++;
            }
            coverage[i] = (double) count / N_SIM;
        }

        System.out.println("Java per-node coverage:");
        for (int i = 0; i < N; i++) {
            System.out.printf("  Node %2d: %.4f%n", i + 1, coverage[i]);
        }

        writeCsv(wsim, OUTPUT_CSV);
        System.out.println("Written: " + OUTPUT_CSV);
    }

    private static void writeCsv(double[][] wsim, String path) throws IOException {
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(path))) {
            // header
            StringBuilder hdr = new StringBuilder();
            for (int i = 0; i < N; i++) {
                if (i > 0) hdr.append(',');
                hdr.append("node").append(i + 1);
            }
            bw.write(hdr.toString());
            bw.newLine();

            // rows
            StringBuilder row = new StringBuilder();
            for (int s = 0; s < N_SIM; s++) {
                row.setLength(0);
                for (int i = 0; i < N; i++) {
                    if (i > 0) row.append(',');
                    row.append(wsim[s][i]);
                }
                bw.write(row.toString());
                bw.newLine();
            }
        }
    }

    private static Calibration[] buildCalibrationSpecs() {
        String[][] taxa = {
                {"t1","t2","t3","t4","t5","t6","t7","t8","t9","t10","t11","t12","t13","t14"},
                {"t1","t2","t3","t4","t5","t6","t7","t8","t9","t10"},
                {"t11","t12","t13","t14"},
                {"t11","t12"},
                {"t1","t2"},
                {"t3","t4","t5","t6","t7","t8","t9","t10"},
                {"t3","t4"},
                {"t5","t6","t7","t8","t9","t10"},
                {"t5","t6","t7","t8","t9"},
                {"t8","t9"},
        };
        Calibration[] specs = new Calibration[N];
        for (int i = 0; i < N; i++) {
            specs[i] = new Calibration(taxa[i], T_HI[i], T_LO[i]);
        }
        return specs;
    }
}
