package calibrationprior;

import java.io.PrintStream;

import beast.base.core.Description;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.Distribution;

/**
 * Wrapper prior for a calibrated tree that holds <em>either</em> a topology-consistent
 * {@link CalibrationPrior} <em>or</em> a set of independent
 * {@link beast.base.spec.evolution.tree.MRCAPrior}s as its children — the two mutually-exclusive
 * calibration modes offered by BEAUti's calibration panel.
 *
 * <p>It is a plain {@link CompoundDistribution}: it simply sums the log-densities of whatever
 * children it currently holds, so neither {@code CalibrationPrior} nor {@code MRCAPrior} needs to
 * know about the other. It exists as its own type only so BEAUti can attach the calibration editor
 * (mode toggle + table) to it while keeping its children off the top-level Priors list.</p>
 */
@Description("Container for a calibrated tree's clade-age prior: holds either a CalibrationPrior "
        + "(topology-consistent) or a set of independent MRCAPriors, and sums its children.")
public class CalibrationDistribution extends CompoundDistribution {

    @Override
    public void init(PrintStream out) {
        super.init(out);
        for (Distribution d : pDistributions.get()) d.init(out);
    }

    @Override
    public void log(long sample, PrintStream out) {
        super.log(sample, out);
        for (Distribution d : pDistributions.get()) d.log(sample, out);
    }

    @Override
    public void close(PrintStream out) {
        super.close(out);
        for (Distribution d : pDistributions.get()) d.close(out);
    }
}
