package calibratedcpp.lphybeast;

import lphybeast.LPhyBeastMain;

import java.util.ArrayList;
import java.util.List;

/**
 * Thin wrapper around {@link LPhyBeastMain} that adds two calibratedcpp-only flags to the
 * {@code convert} subcommand: {@code -MRCAPrior}/{@code --mrcaPrior} and
 * {@code -conditionOnCalibrations}/{@code --conditionOnCalibrations}. Neither flag is known to
 * lphybeast's own picocli parser (that CLI lives in a separate repo we don't own), so this class
 * strips them out of {@code args} before delegating, translating each into the matching System
 * property ({@code calibratedcppMRCAPrior}, {@code calibratedcppConditionOnCalibrations}) that
 * {@link calibratedcpp.lphybeast.tobeast.generators.MRCAPriorCalibrationUtils#isMrcaPriorMode()}
 * / {@code #getConditionOnCalibrationsOverride()} already read.
 *
 * <p>Usage (same {@code convert}/{@code run} subcommands as lphybeast itself, flags can appear
 * anywhere in the argument list):
 * <pre>
 *   convert -MRCAPrior -o out.xml script.lphy
 *   convert -conditionOnCalibrations false -o out.xml script.lphy
 * </pre>
 */
public class CalibratedCPPLPhyBeastMain {

    public static void main(String[] args) {
        List<String> remaining = new ArrayList<>();

        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            if (arg.equals("-MRCAPrior") || arg.equals("--mrcaPrior")) {
                System.setProperty("calibratedcppMRCAPrior", "true");
            } else if (arg.equals("-conditionOnCalibrations") || arg.equals("--conditionOnCalibrations")) {
                if (i + 1 >= args.length) {
                    throw new IllegalArgumentException(arg + " requires a value (true/false)");
                }
                System.setProperty("calibratedcppConditionOnCalibrations", args[++i]);
            } else {
                remaining.add(arg);
            }
        }

        LPhyBeastMain.main(remaining.toArray(new String[0]));
    }
}
