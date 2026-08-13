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
