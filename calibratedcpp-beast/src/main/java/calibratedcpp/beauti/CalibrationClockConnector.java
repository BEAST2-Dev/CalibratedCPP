package calibratedcpp.beauti;

import java.util.ArrayList;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.spec.evolution.likelihood.GenericTreeLikelihood;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.Distribution;
import beast.base.inference.StateNode;
import beast.base.spec.evolution.tree.MRCAPrior;
import beastfx.app.inputeditor.BeautiDoc;
import calibrationprior.CalibrationDistribution;

/**
 * BEAUti "method connector" — invoked by the template via
 * {@code <connect method='calibratedcpp.beauti.CalibrationClockConnector.scrub'/>}. BEAUti calls the
 * named static method (with the {@link BeautiDoc}) on every {@code scrubAll}, as a side effect of
 * evaluating the connector; it has no srcID/targetID, so it connects nothing itself.
 *
 * <p>It does two jobs, both keyed off one condition — whether a partition's
 * {@link CalibrationDistribution} wrapper is currently in the top-level {@code prior} (i.e. a
 * calibrated tree prior is active):
 *
 * <p>It is registered at the {@code aux-partitiontemplate} merge point (see CalibratedCPP.xml), i.e.
 * inside the StandardPartitionTemplate whose {@code mainid='mcmc'} always exists — so it runs on
 * every {@code scrubAll}, first (before the tree-prior and clock-model subtemplates), regardless of
 * which tree prior is active.
 *
 * <ol>
 *   <li><b>Clock-rate flag (Remco's point 1).</b> When the wrapper <em>is</em> active, force the
 *       partition's clock rate to be estimated. Core {@code BeautiDoc.setClockRate} only recognises
 *       {@code MRCAPrior} as timing information, so in CalibrationPrior mode the flag would otherwise
 *       stay {@code false}. The flag alone is not enough: the clock rate is wired into {@code <state>}
 *       (plus operator and ClockPrior) by ClockModels.xml connectors gated on
 *       {@code clockRate.c:$(n)/estimate=true}. {@code setClockRate} resets the flag to false at the
 *       <em>top</em> of every scrub, so this method must set it true <em>before</em> those connectors
 *       run — which is why it lives in the partition template (processed before the clock subtemplate)
 *       rather than in a tree-prior subtemplate (processed after, leaving the flag set too late and
 *       the clock unwired). The rate input is {@code "clock.rate"} for both strict (clockRate) and
 *       Optimised Relaxed (ucldMean) clocks.</li>
 *   <li><b>Stale-logger side effect (Remco's point 3).</b> When the wrapper is <em>not</em> active
 *       (e.g. the user switched to Yule), detach every {@code MRCAPrior} that is no longer a live
 *       child of an active wrapper. They are created imperatively and no connector governs them, so
 *       they keep their {@code tree} reference and linger in {@code tree.getOutputs()} — exactly what
 *       {@code setClockRate} scans, wrongly keeping the flag set after a switch. Detaching discards
 *       only the derived MRCAPrior objects; the source of truth (TaxonSets + CalibrationCladePrior
 *       bounds) is untouched, so calibrations rebuild on demand when switching back.</li>
 * </ol>
 *
 * <p>Note the {@code instanceof} types must be the {@code beast.base.spec.*} classes core actually
 * instantiates ({@code spec…GenericTreeLikelihood}, {@code spec…MRCAPrior}), not their classic twins.
 */
public class CalibrationClockConnector {

    /** Template entry point: {@code <connect method='calibratedcpp.beauti.CalibrationClockConnector.scrub'/>}. */
    public static void scrub(BeautiDoc doc) {
        CompoundDistribution prior =
                (doc.pluginmap.get("prior") instanceof CompoundDistribution c) ? c : null;

        // Snapshot the pluginmap: detachFromTree mutates it via unregisterPlugin.
        for (BEASTInterface bi : new ArrayList<>(doc.pluginmap.values())) {
            // Point 1: estimate the clock whenever a calibration wrapper is active for the partition.
            if (bi instanceof CalibrationDistribution wrapper && isActive(prior, wrapper))
                setClockEstimated(doc, partitionOf(wrapper));

            // Point 3: detach any MRCAPrior that is not a live child of an active wrapper, so core
            // setClockRate stops seeing it via tree.getOutputs().
            if (bi instanceof MRCAPrior mrca && !isLive(prior, mrca))
                detachFromTree(doc, mrca);
        }
    }

    /** True when the wrapper is a current child of the top-level prior (a calibrated prior is active). */
    private static boolean isActive(CompoundDistribution prior, CalibrationDistribution wrapper) {
        return prior != null && prior.pDistributions.get().contains(wrapper);
    }

    /** An MRCAPrior is "live" iff it is a child of a wrapper that is itself active. */
    private static boolean isLive(CompoundDistribution prior, MRCAPrior mrca) {
        if (prior == null) return false;
        for (Distribution d : prior.pDistributions.get())
            if (d instanceof CalibrationDistribution wrapper && wrapper.pDistributions.get().contains(mrca))
                return true;
        return false;
    }

    /** Removes the MRCAPrior from its Tree's output set and the pluginmap, so setClockRate can't see it. */
    private static void detachFromTree(BeautiDoc doc, MRCAPrior mrca) {
        Object tree = inputValue(mrca, "tree");
        if (tree instanceof BEASTInterface t) t.getOutputs().remove(mrca);
        doc.unregisterPlugin(mrca);
    }

    /** Sets the partition's clock-rate StateNode to estimated, mirroring BeautiDoc.setClockRate's traversal. */
    private static void setClockEstimated(BeautiDoc doc, String partition) {
        if (!(doc.pluginmap.get("likelihood") instanceof CompoundDistribution likelihood)) return;
        // 'partition' already carries the "t:" prefix (e.g. "t:align"), so the tree id is "Tree." + it.
        Object tree = doc.pluginmap.get("Tree." + partition);
        for (Distribution d : likelihood.pDistributions.get()) {
            if (!(d instanceof GenericTreeLikelihood tl)) continue;
            if (tree != null && tl.treeInput.get() != tree) continue;
            BEASTInterface clock = tl.branchRateModelInput.get();
            if (clock != null && inputValue(clock, "clock.rate") instanceof StateNode rate)
                rate.isEstimatedInput.setValue(true, rate);
        }
    }

    /** Value of a named input, or {@code null} if the object has no such input. */
    private static Object inputValue(BEASTInterface o, String name) {
        try {
            Input<?> in = o.getInput(name);
            return in == null ? null : in.get();
        } catch (Exception e) {
            return null;
        }
    }

    /** Partition suffix from a wrapper ID (CalibrationDistribution.{p} or CalibrationDistributionADB.{p}). */
    private static String partitionOf(CalibrationDistribution wrapper) {
        return wrapper.getID().replaceFirst("^CalibrationDistribution(ADB)?\\.", "");
    }
}
