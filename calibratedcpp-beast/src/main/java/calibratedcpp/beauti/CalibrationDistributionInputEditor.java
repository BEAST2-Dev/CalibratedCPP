package calibratedcpp.beauti;

import java.util.*;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.inference.CompoundDistribution;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.Real;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.evolution.tree.MRCAPrior;
import beast.base.spec.inference.distribution.Beta;
import beast.base.spec.inference.distribution.Exponential;
import beast.base.spec.inference.distribution.Gamma;
import beast.base.spec.inference.distribution.InverseGamma;
import beast.base.spec.inference.distribution.Laplace;
import beast.base.spec.inference.distribution.LogNormal;
import beast.base.spec.inference.distribution.Normal;
import beast.base.spec.inference.distribution.OffsetReal;
import beast.base.spec.inference.distribution.ScalarDistribution;
import beast.base.spec.inference.distribution.Uniform;
import beast.base.spec.inference.parameter.RealScalarParam;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.InputEditor;
import beastfx.app.util.FXUtils;
import calibratedcpp.CalibratedCoalescentPointProcess;
import calibrationprior.CalibrationCladePrior;
import calibrationprior.CalibrationDistribution;
import calibrationprior.CalibrationPrior;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.*;
import javafx.scene.layout.*;

/**
 * BEAUti editor for {@link CalibrationDistribution} — the wrapper prior that holds either a
 * {@link CalibrationPrior} (topology-consistent clade prior) or a set of independent
 * {@link MRCAPrior}s as its children. Offers a mode toggle plus a single table; only the wrapper
 * appears in the Priors list, so its children never show as separate rows.
 */
public class CalibrationDistributionInputEditor extends InputEditor.Base {

    record DistDef(String name, List<ParamDef> params) {}
    record ParamDef(String key, String label, double defaultVal) {}

    static final List<DistDef> DISTS = List.of(
        new DistDef("Log-Normal", List.of(
            new ParamDef("M",      "M (log)", 1.0),
            new ParamDef("S",      "S (log)", 0.5),
            new ParamDef("offset", "Offset",  0.0))),
        new DistDef("Normal", List.of(
            new ParamDef("mean",   "Mean",   5.0),
            new ParamDef("sigma",  "Sigma",  1.0),
            new ParamDef("offset", "Offset", 0.0))),
        new DistDef("Exponential", List.of(
            new ParamDef("mean",   "Mean",   2.0),
            new ParamDef("offset", "Offset", 0.0))),
        new DistDef("Gamma", List.of(
            new ParamDef("alpha",  "Shape",  2.0),
            new ParamDef("lambda", "Rate",   0.5),
            new ParamDef("offset", "Offset", 0.0))),
        new DistDef("Inverse Gamma", List.of(
            new ParamDef("alpha",  "Alpha",  3.0),
            new ParamDef("beta",   "Beta",   2.0),
            new ParamDef("offset", "Offset", 0.0))),
        new DistDef("Beta", List.of(
            new ParamDef("alpha",  "Alpha",  1.0),
            new ParamDef("beta",   "Beta",   1.0))),
        new DistDef("Uniform", List.of(
            new ParamDef("lower",  "Lower",  0.0),
            new ParamDef("upper",  "Upper", 10.0))),
        new DistDef("Laplace", List.of(
            new ParamDef("mu",     "Mean",   5.0),
            new ParamDef("scale",  "Scale",  1.0),
            new ParamDef("offset", "Offset", 0.0)))
    );

    public CalibrationDistributionInputEditor(BeautiDoc doc) { super(doc); }
    public CalibrationDistributionInputEditor() { super(); }

    @Override
    public Class<?> type() { return CalibrationDistribution.class; }

    @Override
    public void init(Input<?> input, BEASTInterface beastObject, int listItemNr,
                     ExpandOption isExpandOption, boolean addButtons) {
        super.init(input, beastObject, listItemNr, isExpandOption, addButtons);
        pane.getChildren().removeIf(n -> !(n instanceof Label));
        VBox mainBox = new VBox(6);
        mainBox.setMaxWidth(Double.MAX_VALUE);
        HBox.setHgrow(mainBox, Priority.ALWAYS);
        pane.getChildren().add(mainBox);
        rebuild((CalibrationDistribution) beastObject, mainBox);
    }

    private void rebuild(CalibrationDistribution cd, VBox mainBox) {
        mainBox.getChildren().clear();
        String partition = partitionOf(cd);
        CalibratedCoalescentPointProcess cpp = getCPP(partition);

        // After a tree-prior template switch the active model is freshly built with an empty
        // calibrations list, but the TaxonSet.*/.partition plugins survive in the pluginmap.
        if (cpp != null && cpp.calibrationsInput.get().isEmpty())
            recoverCalibrationTaxonSets(partition, cpp.calibrationsInput.get());

        if (cpp == null || cpp.calibrationsInput.get().isEmpty()) {
            mainBox.getChildren().add(new Label("No calibration clades defined. Use 'Manage calibrations' to add them."));
            return;
        }

        List<TaxonSet> clades = new ArrayList<>(cpp.calibrationsInput.get());
        boolean mrcaMode = hasMRCAPriors(cd);

        ToggleGroup tg = new ToggleGroup();
        RadioButton calRb  = new RadioButton("CalibrationPrior");
        RadioButton mrcaRb = new RadioButton("MRCAPriors");
        calRb.setToggleGroup(tg);
        mrcaRb.setToggleGroup(tg);
        calRb.setSelected(!mrcaMode);
        mrcaRb.setSelected(mrcaMode);
        HBox modeRow = new HBox(12, new Label("Calibration Prior:"), calRb, mrcaRb);
        modeRow.setPadding(new Insets(4, 0, 8, 0));
        modeRow.setAlignment(Pos.CENTER_LEFT);
        mainBox.getChildren().add(modeRow);

        // Keep the CalibrationPrior child + clades consistent on every render (not just on click).
        if (!mrcaMode) syncCalibrationPriorClades(cd, clades, partition);

        VBox content = FXUtils.newVBox();
        content.setSpacing(4);
        mainBox.getChildren().add(content);
        renderRows(content, cd, clades, partition, !mrcaMode);

        calRb.setOnAction(e -> {
            switchToCalibrationPrior(cd, clades, partition);
            renderRows(content, cd, clades, partition, true);
        });
        mrcaRb.setOnAction(e -> {
            switchToMRCAPriors(cd, clades, partition);
            renderRows(content, cd, clades, partition, false);
        });
    }

    // ── Row rendering ─────────────────────────────────────────────────────────────

    private void renderRows(VBox box, CalibrationDistribution cd, List<TaxonSet> clades,
                             String partition, boolean calMode) {
        box.getChildren().clear();

        // CalibrationPrior is a pure clade-density prior: every clade must have both a lower and an
        // upper bound. Flag any clade still missing one so the requirement is visible rather than a
        // clade silently contributing nothing to the density.
        if (calMode) {
            List<String> missing = new ArrayList<>();
            for (TaxonSet ts : clades) {
                String label = labelOf(ts);
                if (CladeBounds.lowerOf(doc, label, partition) == null
                        || CladeBounds.upperOf(doc, label, partition) == null)
                    missing.add(label);
            }
            if (!missing.isEmpty()) {
                Label warn = new Label("⚠ CalibrationPrior requires a lower and an upper bound for "
                    + "every clade. Missing bounds: " + String.join(", ", missing) + ".");
                warn.setWrapText(true);
                warn.setStyle("-fx-text-fill: #b00020; -fx-font-size: 11px; -fx-font-weight: bold;");
                warn.setPadding(new Insets(0, 0, 4, 0));
                box.getChildren().add(warn);
            }
        }

        HBox hdr = new HBox(8);
        hdr.getChildren().add(fixedLabel("Clade", 160, true));
        if (calMode) {
            hdr.getChildren().addAll(fixedLabel("Lower bound", 100, true),
                                     fixedLabel("Upper bound", 100, true),
                                     fixedLabel("Confidence level", 100, true));
        } else {
            hdr.getChildren().addAll(fixedLabel("Distribution", 130, true),
                                     fixedLabel("Parameters", 300, true),
                                     fixedLabel("Monophyletic", 90, true));
        }
        box.getChildren().add(hdr);
        box.getChildren().add(new Separator());

        // Re-render after a cal-mode edit so the missing-bounds banner and the per-field required
        // highlighting reflect the value just entered.
        Runnable rerender = () -> renderRows(box, cd, clades, partition, calMode);

        for (TaxonSet ts : clades) {
            HBox row = new HBox(8);
            row.setPadding(new Insets(3, 0, 3, 0));
            row.setAlignment(Pos.CENTER_LEFT);
            row.getChildren().add(fixedLabel(labelOf(ts), 160, false));
            if (calMode) buildCalRow(row, cd, ts, partition, rerender);
            else         buildMRCARow(row, cd, ts, partition);
            box.getChildren().add(row);
        }
    }

    private void buildCalRow(HBox row, CalibrationDistribution cd, TaxonSet ts, String partition,
                             Runnable rerender) {
        String label = labelOf(ts);
        CalibrationCladePrior ccp = findCladePrior(ts, partition);
        Double lo = CladeBounds.lowerOf(doc, label, partition);
        Double hi = CladeBounds.upperOf(doc, label, partition);
        TextField loTf = numField(lo != null ? fmt(lo) : "");
        TextField hiTf = numField(hi != null ? fmt(hi) : "");
        TextField clTf = numField(ccp != null ? fmt(ccp.getCoverage()) : "0.9");
        loTf.setPrefWidth(100);
        hiTf.setPrefWidth(100);
        clTf.setPrefWidth(100); clTf.setPromptText("0–1");
        // Both bounds are required for a CalibrationPrior; mark the ones still empty.
        markRequired(loTf, lo == null);
        markRequired(hiTf, hi == null);
        clTf.setTooltip(new Tooltip("Probability mass within the bounds (between 0 and 1). Default 0.9."));
        boolean loWasEmpty = lo == null;
        boolean hiWasEmpty = hi == null;
        Runnable apply = () -> {
            applyCladePrior(cd, ts, partition, loTf.getText(), hiTf.getText(), clTf.getText());
            // Only rebuild when a bound has appeared or disappeared — that is when the banner and the
            // required-field highlighting change. Skipping it otherwise avoids stealing focus on blur.
            if (loTf.getText().trim().isEmpty() != loWasEmpty
                    || hiTf.getText().trim().isEmpty() != hiWasEmpty)
                rerender.run();
        };
        loTf.focusedProperty().addListener((obs, o, n) -> { if (!n) apply.run(); });
        hiTf.focusedProperty().addListener((obs, o, n) -> { if (!n) apply.run(); });
        clTf.focusedProperty().addListener((obs, o, n) -> { if (!n) apply.run(); });
        row.getChildren().addAll(loTf, hiTf, clTf);
    }

    /** Marks a bound field as required (empty) or satisfied, restyling its border and prompt. */
    private static void markRequired(TextField tf, boolean required) {
        if (required) {
            tf.setPromptText("required");
            tf.setStyle("-fx-font-size: 12px; -fx-border-color: #d9534f; -fx-border-width: 1;");
        } else {
            tf.setPromptText("");
            tf.setStyle("-fx-font-size: 12px;");
        }
    }

    private void buildMRCARow(HBox row, CalibrationDistribution cd, TaxonSet ts, String partition) {
        MRCAPrior mrca = findMRCAPrior(ts, partition);
        DistDef curDef = findDef(detectDistName(mrca));
        Double lowerBound = cladeLowerBound(ts, partition);

        ChoiceBox<String> distCb = new ChoiceBox<>();
        DISTS.forEach(d -> distCb.getItems().add(d.name()));
        distCb.setValue(curDef.name());
        distCb.setPrefWidth(130);

        HBox paramsBox = new HBox(6);
        paramsBox.setAlignment(Pos.CENTER_LEFT);
        renderParamFields(paramsBox, curDef, mrca, lowerBound);

        CheckBox monoCb = new CheckBox();
        monoCb.setSelected(mrca == null || Boolean.TRUE.equals(mrca.getInputValue("monophyletic")));

        Runnable apply = () -> applyMRCAPrior(cd, ts, partition, findDef(distCb.getValue()),
                                               readParams(paramsBox), monoCb.isSelected());
        attachParamListeners(paramsBox, apply);
        monoCb.setOnAction(e -> apply.run());

        distCb.setOnAction(e -> {
            DistDef nd = findDef(distCb.getValue());
            renderParamFields(paramsBox, nd, null, lowerBound);
            attachParamListeners(paramsBox, apply);
            applyMRCAPrior(cd, ts, partition, nd, defaultParams(nd, lowerBound), monoCb.isSelected());
        });

        row.getChildren().addAll(distCb, paramsBox, monoCb);
    }

    private void renderParamFields(HBox paramsBox, DistDef def, MRCAPrior existingMrca, Double lowerBound) {
        paramsBox.getChildren().clear();
        for (ParamDef pd : def.params()) {
            double dflt = pd.key().equals("offset") && lowerBound != null ? lowerBound : pd.defaultVal();
            double val = getDistParam(existingMrca, def.name(), pd.key(), dflt);
            Label lbl = new Label(pd.label() + ":");
            lbl.setStyle("-fx-font-size:11px");
            TextField tf = numField(fmt(val));
            tf.setPrefWidth(72);
            tf.setUserData(pd.key());
            paramsBox.getChildren().addAll(lbl, tf);
        }
    }

    private void attachParamListeners(HBox paramsBox, Runnable apply) {
        for (var node : paramsBox.getChildren()) {
            if (node instanceof TextField tf) {
                tf.focusedProperty().addListener((obs, o, n) -> { if (!n) apply.run(); });
            }
        }
    }

    private Map<String, Double> readParams(HBox paramsBox) {
        Map<String, Double> result = new LinkedHashMap<>();
        for (var node : paramsBox.getChildren()) {
            if (node instanceof TextField tf && tf.getUserData() instanceof String key) {
                try { result.put(key, Double.parseDouble(tf.getText().trim())); } catch (NumberFormatException ignored) {}
            }
        }
        return result;
    }

    private static Map<String, Double> defaultParams(DistDef def) {
        Map<String, Double> m = new LinkedHashMap<>();
        def.params().forEach(pd -> m.put(pd.key(), pd.defaultVal()));
        return m;
    }

    /** Default params, but with the offset seeded from the clade's lower bound when one exists. */
    private static Map<String, Double> defaultParams(DistDef def, Double lowerBound) {
        Map<String, Double> m = defaultParams(def);
        if (lowerBound != null && m.containsKey("offset")) m.put("offset", lowerBound);
        return m;
    }

    /**
     * The clade's lower-bound age — from its CalibrationCladePrior when both bounds are set, else
     * from a partial CladeBounds holder (a lone lower bound), or null if none is provided. Used to
     * seed the MRCAPrior offset.
     */
    private Double cladeLowerBound(TaxonSet ts, String partition) {
        Double lo = CladeBounds.lowerOf(doc, labelOf(ts), partition);
        return (lo == null || Double.isNaN(lo)) ? null : lo;
    }

    // ── Mode switching ────────────────────────────────────────────────────────────

    private void switchToCalibrationPrior(CalibrationDistribution cd, List<TaxonSet> clades, String partition) {
        removeMRCAPriors(cd, partition);
        syncCalibrationPriorClades(cd, clades, partition);
    }

    private void syncCalibrationPriorClades(CalibrationDistribution cd, List<TaxonSet> clades, String partition) {
        CalibrationPrior cp = calibrationPriorOf(cd, partition);
        cp.cladesInput.get().clear();
        for (TaxonSet ts : clades) {
            CalibrationCladePrior ccp = findCladePrior(ts, partition);
            if (ccp != null) cp.cladesInput.get().add(ccp);
        }
        if (!cd.pDistributions.get().contains(cp)) cd.pDistributions.get().add(cp);
        ensureWrapperConnected(cd);
    }

    private void switchToMRCAPriors(CalibrationDistribution cd, List<TaxonSet> clades, String partition) {
        // Drop the CalibrationPrior child (kept in the pluginmap so a toggle back restores it) and
        // seed one Uniform(lower, upper) MRCAPrior per clade from its constraint bounds.
        CalibrationPrior cp = calibrationPriorOf(cd, partition);
        cp.cladesInput.get().clear();
        cd.pDistributions.get().remove(cp);

        DistDef uniform = findDef("Uniform");
        for (TaxonSet ts : clades) {
            Map<String, Double> params = defaultParams(uniform);
            // Seed the Uniform from whichever bounds the clade has (full CCP or a partial holder);
            // a missing side keeps the Uniform default rather than being dropped.
            String label = labelOf(ts);
            Double lo = CladeBounds.lowerOf(doc, label, partition);
            Double hi = CladeBounds.upperOf(doc, label, partition);
            if (lo != null) params.put("lower", lo);
            if (hi != null) params.put("upper", hi);
            applyMRCAPrior(cd, ts, partition, uniform, params, true);
        }
        ensureWrapperConnected(cd);
    }

    private void removeMRCAPriors(CalibrationDistribution cd, String partition) {
        String pfx = "MRCAPrior.", sfx = "." + partition;
        for (String id : new ArrayList<>(doc.pluginmap.keySet())) {
            if (id.startsWith(pfx) && id.endsWith(sfx) && doc.pluginmap.get(id) instanceof MRCAPrior mrca) {
                cd.pDistributions.get().remove(mrca);
            }
        }
    }

    // ── Apply actions ─────────────────────────────────────────────────────────────

    private void applyCladePrior(CalibrationDistribution cd, TaxonSet ts, String partition,
                                  String loStr, String hiStr, String pcovStr) {
        String label = labelOf(ts);
        CalibrationPrior cp = calibrationPriorOf(cd, partition);
        String ccpId = "CalibrationCladePrior." + label + "." + partition;
        CalibrationCladePrior existing = (doc.pluginmap.get(ccpId) instanceof CalibrationCladePrior c) ? c : null;
        cp.cladesInput.get().remove(existing);

        loStr = loStr.trim(); hiStr = hiStr.trim();
        pcovStr = pcovStr.trim();
        if (loStr.isEmpty() || hiStr.isEmpty()) {
            // A CalibrationCladePrior needs both bounds; keep a lone bound in a CladeBounds holder so
            // it is not lost (it still feeds the MRCAPrior offset and pre-fills these fields).
            if (existing != null) doc.pluginmap.remove(ccpId);
            Double lone = null;
            try { if (!loStr.isEmpty()) lone = Double.parseDouble(loStr); } catch (NumberFormatException ignored) {}
            Double loneHi = null;
            try { if (!hiStr.isEmpty()) loneHi = Double.parseDouble(hiStr); } catch (NumberFormatException ignored) {}
            CladeBounds.store(doc, label, partition, lone, loneHi);
            return;
        }
        // Both bounds given: a full CCP supersedes any partial holder.
        CladeBounds.remove(doc, label, partition);
        try {
            double lo = Double.parseDouble(loStr);
            double hi = Double.parseDouble(hiStr);
            double pcov = Double.parseDouble(pcovStr);
            CalibrationCladePrior ccp;
            if (existing != null) {
                ccp = existing;
                ccp.lowerAgeInput.setValue(new RealScalarParam<>(lo, NonNegativeReal.INSTANCE), ccp);
                ccp.upperAgeInput.setValue(new RealScalarParam<>(hi, NonNegativeReal.INSTANCE), ccp);
                ccp.pCoverageInput.setValue(new RealScalarParam<>(pcov, UnitInterval.INSTANCE), ccp);
                try { ccp.initAndValidate(); } catch (Exception ignored) {}
            } else {
                ccp = new CalibrationCladePrior();
                ccp.initByName("taxa", ts,
                    "lowerAge", new RealScalarParam<>(lo, NonNegativeReal.INSTANCE),
                    "upperAge", new RealScalarParam<>(hi, NonNegativeReal.INSTANCE),
                        "confidenceLevel", new RealScalarParam<UnitInterval>(pcov, UnitInterval.INSTANCE));
                ccp.setID(ccpId);
                doc.addPlugin(ccp);
            }
            if (!cp.cladesInput.get().contains(ccp))
                cp.cladesInput.get().add(ccp);
        } catch (Exception ignored) {}
    }

    private void applyMRCAPrior(CalibrationDistribution cd, TaxonSet ts, String partition, DistDef def,
                                  Map<String, Double> params, boolean monophyletic) {
        CalibratedCoalescentPointProcess cpp = getCPP(partition);
        if (cpp == null) return;

        String mrcaId = "MRCAPrior." + labelOf(ts) + "." + partition;
        String distId = "MRCAPriorDist." + labelOf(ts) + "." + partition;

        if (doc.pluginmap.get(mrcaId) instanceof MRCAPrior old) {
            cd.pDistributions.get().remove(old);
        }
        doc.pluginmap.remove(mrcaId);
        doc.pluginmap.remove(distId);

        ScalarDistribution dist = buildDistribution(def, params);
        if (dist != null) { dist.setID(distId); doc.addPlugin(dist); }

        MRCAPrior mrca = new MRCAPrior();
        mrca.setInputValue("tree",        cpp.treeInput.get());
        mrca.setInputValue("taxonset",    ts);
        mrca.setInputValue("monophyletic", monophyletic);
        if (dist != null) mrca.setInputValue("distr", dist);
        try { mrca.initAndValidate(); } catch (Exception e) { e.printStackTrace(); }
        mrca.setID(mrcaId);
        doc.addPlugin(mrca);
        if (!cd.pDistributions.get().contains(mrca)) cd.pDistributions.get().add(mrca);
        ensureWrapperConnected(cd);
    }

    // ── Distribution builder ──────────────────────────────────────────────────────

    private ScalarDistribution buildDistribution(DistDef def, Map<String, Double> params) {
        try {
            // The spec distributions have NO "offset" input (unlike the legacy BEAST ones); an offset
            // is applied by wrapping the distribution in an OffsetReal. Build the base distribution
            // here (never setting a non-existent "offset" input, which would throw and silently drop
            // the whole distribution), then wrap it below if an offset was requested.
            ScalarDistribution base = switch (def.name()) {
                case "Log-Normal" -> {
                    LogNormal d = new LogNormal();
                    d.setInputValue("M", new RealScalarParam<>(params.getOrDefault("M", 1.0), Real.INSTANCE));
                    d.setInputValue("S", new RealScalarParam<>(params.getOrDefault("S", 0.5), PositiveReal.INSTANCE));
                    d.initAndValidate();
                    yield d;
                }
                case "Normal" -> {
                    Normal d = new Normal();
                    d.setInputValue("mean",  new RealScalarParam<>(params.getOrDefault("mean",  5.0), Real.INSTANCE));
                    d.setInputValue("sigma", new RealScalarParam<>(params.getOrDefault("sigma", 1.0), PositiveReal.INSTANCE));
                    d.initAndValidate();
                    yield d;
                }
                case "Exponential" -> {
                    Exponential d = new Exponential();
                    d.setInputValue("mean", new RealScalarParam<>(params.getOrDefault("mean", 2.0), PositiveReal.INSTANCE));
                    d.initAndValidate();
                    yield d;
                }
                case "Gamma" -> {
                    Gamma d = new Gamma();
                    // Gamma's rate parameter is "lambda" (Shape–Rate form); it has no "beta" input.
                    d.setInputValue("alpha",  new RealScalarParam<>(params.getOrDefault("alpha",  2.0), PositiveReal.INSTANCE));
                    d.setInputValue("lambda", new RealScalarParam<>(params.getOrDefault("lambda", 0.5), PositiveReal.INSTANCE));
                    d.initAndValidate();
                    yield d;
                }
                case "Inverse Gamma" -> {
                    InverseGamma d = new InverseGamma();
                    d.setInputValue("alpha", new RealScalarParam<>(params.getOrDefault("alpha", 3.0), PositiveReal.INSTANCE));
                    d.setInputValue("beta",  new RealScalarParam<>(params.getOrDefault("beta",  2.0), PositiveReal.INSTANCE));
                    d.initAndValidate();
                    yield d;
                }
                case "Beta" -> {
                    Beta d = new Beta();
                    d.setInputValue("alpha", new RealScalarParam<>(params.getOrDefault("alpha", 1.0), PositiveReal.INSTANCE));
                    d.setInputValue("beta",  new RealScalarParam<>(params.getOrDefault("beta",  1.0), PositiveReal.INSTANCE));
                    d.initAndValidate();
                    yield d;
                }
                case "Uniform" -> {
                    Uniform d = new Uniform();
                    d.setInputValue("lower", new RealScalarParam<>(params.getOrDefault("lower",  0.0), Real.INSTANCE));
                    d.setInputValue("upper", new RealScalarParam<>(params.getOrDefault("upper", 10.0), Real.INSTANCE));
                    d.initAndValidate();
                    yield d;
                }
                case "Laplace" -> {
                    Laplace d = new Laplace();
                    d.setInputValue("mu",    new RealScalarParam<>(params.getOrDefault("mu",    5.0), Real.INSTANCE));
                    d.setInputValue("scale", new RealScalarParam<>(params.getOrDefault("scale", 1.0), PositiveReal.INSTANCE));
                    d.initAndValidate();
                    yield d;
                }
                default -> null;
            };
            if (base == null) return null;

            double offset = params.getOrDefault("offset", 0.0);
            if (offset != 0.0 && defHasOffset(def)) {
                OffsetReal wrapped = new OffsetReal();
                wrapped.setInputValue("distribution", base);
                wrapped.setInputValue("offset", new RealScalarParam<>(offset, Real.INSTANCE));
                wrapped.initAndValidate();
                return wrapped;
            }
            return base;
        } catch (Exception e) { e.printStackTrace(); return null; }
    }

    private static boolean defHasOffset(DistDef def) {
        return def.params().stream().anyMatch(p -> p.key().equals("offset"));
    }

    /** Unwraps an OffsetReal to its inner distribution; returns the argument unchanged otherwise. */
    private static Object innerDist(Object distr) {
        if (distr instanceof OffsetReal or) {
            try { return or.getInputValue("distribution"); } catch (Exception ignored) {}
        }
        return distr;
    }

    // ── Lookups / connection helpers ───────────────────────────────────────────────

    /** The CalibrationPrior for this partition — the wrapper's current child, else the
     *  disconnected one kept in the pluginmap across toggles, else a freshly created one. */
    private CalibrationPrior calibrationPriorOf(CalibrationDistribution cd, String partition) {
        for (var d : cd.pDistributions.get())
            if (d instanceof CalibrationPrior c) return c;
        String sfx = "." + partition;
        for (BEASTInterface bi : doc.pluginmap.values())
            if (bi instanceof CalibrationPrior c && c.getID() != null && c.getID().endsWith(sfx))
                return c;
        CalibratedCoalescentPointProcess cpp = getCPP(partition);
        CalibrationPrior c = new CalibrationPrior();
        if (cpp != null) c.setInputValue("tree", cpp.treeInput.get());
        c.setID(cd.getID().replace("CalibrationDistribution", "CalibrationPrior"));
        try { c.initAndValidate(); } catch (Exception ignored) {}
        doc.addPlugin(c);
        return c;
    }

    private void ensureWrapperConnected(CalibrationDistribution cd) {
        CompoundDistribution prior = topLevelPrior();
        if (prior != null && !prior.pDistributions.get().contains(cd))
            prior.pDistributions.get().add(cd);
    }

    private CompoundDistribution topLevelPrior() {
        return (doc.pluginmap.get("prior") instanceof CompoundDistribution c) ? c : null;
    }

    private boolean hasMRCAPriors(CalibrationDistribution cd) {
        return cd.pDistributions.get().stream().anyMatch(d -> d instanceof MRCAPrior);
    }

    private CalibrationCladePrior findCladePrior(TaxonSet ts, String partition) {
        String id = "CalibrationCladePrior." + labelOf(ts) + "." + partition;
        return (doc.pluginmap.get(id) instanceof CalibrationCladePrior c) ? c : null;
    }

    private MRCAPrior findMRCAPrior(TaxonSet ts, String partition) {
        String id = "MRCAPrior." + labelOf(ts) + "." + partition;
        return (doc.pluginmap.get(id) instanceof MRCAPrior m) ? m : null;
    }

    private String detectDistName(MRCAPrior mrca) {
        if (mrca == null) return "Uniform";
        Object d = innerDist(mrca.getInputValue("distr"));
        if (d instanceof LogNormal)    return "Log-Normal";
        if (d instanceof Normal)       return "Normal";
        if (d instanceof Exponential)  return "Exponential";
        if (d instanceof Gamma)        return "Gamma";
        if (d instanceof InverseGamma) return "Inverse Gamma";
        if (d instanceof Beta)         return "Beta";
        if (d instanceof Uniform)      return "Uniform";
        if (d instanceof Laplace)      return "Laplace";
        return "Uniform";
    }

    private double getDistParam(MRCAPrior mrca, String distName, String key, double defaultVal) {
        if (mrca == null) return defaultVal;
        Object distr = mrca.getInputValue("distr");
        Object target = key.equals("offset") ? distr : innerDist(distr);
        if (!(target instanceof BEASTInterface bi)) return defaultVal;
        try {
            Object v = bi.getInputValue(key);
            if (v instanceof Number n)            return n.doubleValue();
            if (v instanceof RealScalarParam<?> r) return (double) r.get(0);
        } catch (Exception ignored) {}
        return defaultVal;
    }

    private DistDef findDef(String name) {
        return DISTS.stream().filter(d -> d.name().equals(name)).findFirst()
               .orElse(DISTS.stream().filter(d -> d.name().equals("Uniform")).findFirst().orElse(DISTS.get(0)));
    }

    private static String labelOf(TaxonSet ts) {
        String id = ts.getID();
        if (id == null) return "unknown";
        // ID format: TaxonSet.{label}.{partition}
        int first = id.indexOf('.');
        int second = first >= 0 ? id.indexOf('.', first + 1) : -1;
        return second >= 0 ? id.substring(first + 1, second) : id;
    }

    private static String partitionOf(CalibrationDistribution cd) {
        // Skyline template: CalibrationDistribution.{partition}; age-dependent: CalibrationDistributionADB.{partition}.
        return cd.getID().replaceFirst("^CalibrationDistribution(ADB)?\\.", "");
    }

    /** Finds the active calibrated tree-prior instance (of any concrete type) for this partition. */
    private CalibratedCoalescentPointProcess getCPP(String partition) {
        String suffix = "." + partition;
        if (doc.pluginmap.get("prior") instanceof CompoundDistribution cd) {
            for (var d : cd.pDistributions.get()) {
                if (d instanceof CalibratedCoalescentPointProcess m
                        && m.getID() != null && m.getID().endsWith(suffix))
                    return m;
            }
        }
        for (BEASTInterface bi : doc.pluginmap.values()) {
            if (bi instanceof CalibratedCoalescentPointProcess m && m.getID() != null && m.getID().endsWith(suffix))
                return m;
        }
        return null;
    }

    private void recoverCalibrationTaxonSets(String partition, List<TaxonSet> target) {
        String prefix = "TaxonSet.";
        String suffix = "." + partition;
        int minLen = prefix.length() + suffix.length() + 1;
        for (BEASTInterface bi : doc.pluginmap.values()) {
            if (bi instanceof TaxonSet ts && ts.getID() != null
                    && ts.getID().startsWith(prefix) && ts.getID().endsWith(suffix)
                    && ts.getID().length() >= minLen
                    && !target.contains(ts)) {
                target.add(ts);
            }
        }
    }

    // ── UI helpers ────────────────────────────────────────────────────────────────

    private static Label fixedLabel(String text, double width, boolean bold) {
        Label lbl = new Label(text);
        lbl.setMinWidth(width);
        lbl.setPrefWidth(width);
        if (bold) lbl.setStyle("-fx-font-weight: bold; -fx-font-size: 11px;");
        else      lbl.setStyle("-fx-font-size: 12px;");
        return lbl;
    }

    private static TextField numField(String text) {
        TextField tf = new TextField(text);
        tf.setStyle("-fx-font-size: 12px;");
        return tf;
    }

    private static String fmt(double v) {
        if (Double.isNaN(v)) return "";
        return v == Math.floor(v) && !Double.isInfinite(v) ? String.valueOf((long) v) : String.valueOf(v);
    }
}
