package calibratedcpp.beauti;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beastfx.app.inputeditor.BeautiDoc;
import calibrationprior.CalibrationCladePrior;

@Description("Editor-only holder for a partially specified clade age bound.")
public class CladeBounds extends BEASTObject {

    public Double lower;
    public Double upper;

    @Override
    public void initAndValidate() { }

    static String idFor(String label, String partition) {
        return "CladeBounds." + label + "." + partition;
    }

    static CladeBounds find(BeautiDoc doc, String label, String partition) {
        return doc.pluginmap.get(idFor(label, partition)) instanceof CladeBounds cb ? cb : null;
    }

    static Double lowerOf(BeautiDoc doc, String label, String partition) {
        if (doc.pluginmap.get("CalibrationCladePrior." + label + "." + partition)
                instanceof CalibrationCladePrior c) return c.getLower();
        CladeBounds cb = find(doc, label, partition);
        return cb != null ? cb.lower : null;
    }

    static Double upperOf(BeautiDoc doc, String label, String partition) {
        if (doc.pluginmap.get("CalibrationCladePrior." + label + "." + partition)
                instanceof CalibrationCladePrior c) return c.getUpper();
        CladeBounds cb = find(doc, label, partition);
        return cb != null ? cb.upper : null;
    }


    static void store(BeautiDoc doc, String label, String partition, Double lower, Double upper) {
        String id = idFor(label, partition);
        if (lower == null && upper == null) { doc.pluginmap.remove(id); return; }
        CladeBounds cb = find(doc, label, partition);
        boolean isNew = cb == null;
        if (isNew) { cb = new CladeBounds(); cb.setID(id); }
        cb.lower = lower;
        cb.upper = upper;
        if (isNew) doc.addPlugin(cb);
    }

    static void remove(BeautiDoc doc, String label, String partition) {
        doc.pluginmap.remove(idFor(label, partition));
    }
}
