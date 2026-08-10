import calibratedcpp.lphy.prior.Calibration;
import calibratedcpp.lphy.prior.CalibrationArray;
import calibratedcpp.lphy.prior.ConditionedMRCAPrior;
import lphy.core.model.Value;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Samples ConditionedMRCAPrior N_SIM times using the calibrations from
 * validation/calibrationprior/test_calibration_prior.xml
 * and writes output as a TSV whose columns match the BEAST MCMC log
 * (mrca.clade1 … mrca.clade10).
 *
 * Run main() from IDE or with:
 *   mvn exec:java -pl calibratedcpp-lphy \
 *     -Dexec.mainClass=CalibrationPriorValidationSampler \
 *     -Dexec.classpathScope=test
 *
 * Then run calibration_prior_comparison.R to visualise.
 */
public class CalibrationPriorValidationSampler {

    private static final int N_SIM = 50000;

    // Taxa in the order passed to ConditionedMRCAPrior (root first, then by clade hierarchy).
    // Clade numbering matches calibrationprior_simulation.r (clade1 = root).
    private static final String[][] TAXA = {
        {"leaf_1","leaf_2","leaf_3","leaf_4","leaf_5","leaf_6","leaf_7","leaf_8","leaf_9","leaf_10","leaf_11"}, // clade1 (root)
        {"leaf_1","leaf_2","leaf_3","leaf_4","leaf_5","leaf_6","leaf_7","leaf_8"},  // clade2
        {"leaf_1","leaf_2","leaf_3","leaf_4","leaf_5","leaf_6"},                    // clade6
        {"leaf_1","leaf_2","leaf_3","leaf_4"},                                      // clade8
        {"leaf_1","leaf_2","leaf_3"},                                               // clade9
        {"leaf_1","leaf_2"},                                                        // clade10
        {"leaf_5","leaf_6"},                                                        // clade7
        {"leaf_7","leaf_8"},                                                        // clade5
        {"leaf_9","leaf_10","leaf_11"},                                             // clade3
        {"leaf_9","leaf_10"},                                                       // clade4
    };

    // Bounds in same order as TAXA
    private static final Double[] LOWER = {10.0, 9.4, 6.8, 2.0, 1.7, 1.6, 1.8, 8.0, 4.8, 4.0};
    private static final Double[] UPPER = {10.5, 9.6, 8.0, 3.2, 2.5, 2.5, 2.5, 9.5, 5.0, 5.0};

    // Maps LPhy output index → BEAST column name
    private static final String[] COL_NAMES = {
        "mrca.clade1", "mrca.clade2", "mrca.clade6", "mrca.clade8",
        "mrca.clade9", "mrca.clade10", "mrca.clade7", "mrca.clade5",
        "mrca.clade3", "mrca.clade4"
    };

    private static final String OUTPUT = "../validation/calibrationprior/lphy_wsim.tsv";

    public static void main(String[] args) throws IOException {
        int N = TAXA.length;
        Calibration[] specs = new Calibration[N];
        for (int i = 0; i < N; i++) {
            specs[i] = new Calibration(TAXA[i], UPPER[i], LOWER[i]);
        }
        ConditionedMRCAPrior prior = new ConditionedMRCAPrior(
                new Value<>("calibrations", specs),
                new Value<>("p", 0.90)
        );

        double[][] wsim = new double[N_SIM][N];

        for (int s = 0; s < N_SIM; s++) {
            CalibrationArray ca = prior.sample().value();
            Calibration[] cals = ca.getCalibrationArray();
            for (int i = 0; i < N; i++) {
                wsim[s][i] = cals[i].getAge();
            }
        }

        // per-node coverage
        System.out.println("LPhy per-node coverage:");
        for (int i = 0; i < N; i++) {
            int count = 0;
            for (int s = 0; s < N_SIM; s++) {
                if (wsim[s][i] >= LOWER[i] && wsim[s][i] <= UPPER[i]) count++;
            }
            System.out.printf("  %-15s: %.4f%n", COL_NAMES[i], (double) count / N_SIM);
        }

        writeTsv(wsim, N, OUTPUT);
        System.out.println("Written: " + OUTPUT);
    }

    private static void writeTsv(double[][] wsim, int N, String path) throws IOException {
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(path))) {
            bw.write(String.join("\t", COL_NAMES));
            bw.newLine();
            StringBuilder row = new StringBuilder();
            for (int s = 0; s < N_SIM; s++) {
                row.setLength(0);
                for (int i = 0; i < N; i++) {
                    if (i > 0) row.append('\t');
                    row.append(wsim[s][i]);
                }
                bw.write(row.toString());
                bw.newLine();
            }
        }
    }
}
