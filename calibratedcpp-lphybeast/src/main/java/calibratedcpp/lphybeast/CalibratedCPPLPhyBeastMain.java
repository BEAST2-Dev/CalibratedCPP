package calibratedcpp.lphybeast;

import lphybeast.LPhyBeastMain;

import java.util.ArrayList;
import java.util.List;

/**
 * Wraps {@link LPhyBeastMain} to add a calibratedcpp-only flag that its picocli parser (a
 * separate repo) does not know: it is stripped from {@code args} and translated into the
 * system property {@link calibratedcpp.lphybeast.tobeast.generators.MRCAPriorCalibrationUtils}
 * reads. The flag may appear anywhere in the argument list:
 * <pre>
 *   convert -conditionOnCalibrations false script.lphy
 * </pre>
 */
public class CalibratedCPPLPhyBeastMain {

    public static void main(String[] args) {
        List<String> remaining = new ArrayList<>();

        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            if (arg.equals("-conditionOnCalibrations") || arg.equals("--conditionOnCalibrations")) {
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
