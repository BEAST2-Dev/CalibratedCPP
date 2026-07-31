package calibratedcpp.beauti;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.inference.StateNode;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.distribution.IID;
import beast.base.spec.inference.distribution.LogNormal;
import beast.base.spec.inference.distribution.Uniform;
import beast.base.spec.inference.operator.ScaleOperator;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.Tensor;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.util.FXUtils;
import calibratedcpp.CalibratedBirthDeathSkylineModel;
import calibratedcpp.CalibratedCoalescentPointProcess;
import calibratedcpp.SkylineParameter;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Node;
import javafx.scene.control.Button;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ChoiceBox;
import javafx.scene.control.Label;
import javafx.scene.control.Spinner;
import javafx.scene.control.TextField;
import javafx.scene.control.Tooltip;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Pane;
import javafx.scene.layout.VBox;

/**
 * BEAUti panel for {@link CalibratedBirthDeathSkylineModel}. Inherits the calibration-management
 * popup, live preview, conditioning row, and model persistence from
 * {@link CalibratedCPPInputEditor}; this class contributes only the skyline
 * rate-parameterization UI (a rate-pair picker plus per-rate epoch/change-time editors) and the
 * sampling-proportion row.
 */
public class CalibratedBirthDeathSkylineInputEditor extends CalibratedCPPInputEditor {

    // ── Parameterization tables ───────────────────────────────────────────────────

    private static final String[][] PARAM_PAIRS = {
        {"Diversification rate + Turnover",            "diversificationRate", "turnover"},
        {"Birth rate + Death rate",                    "birthRate",           "deathRate"},
        {"Birth rate + Turnover",                      "birthRate",           "turnover"},
        {"Birth rate + Reproductive number",           "birthRate",           "reproductiveNumber"},
        {"Diversification rate + Death rate",          "diversificationRate", "deathRate"},
        {"Diversification rate + Reproductive number", "diversificationRate", "reproductiveNumber"},
        {"Death rate + Turnover",                      "deathRate",           "turnover"},
        {"Death rate + Reproductive number",           "deathRate",           "reproductiveNumber"},
        {"Birth rate + Diversification rate",          "birthRate",           "diversificationRate"},
    };

    private static final String[] ALL_RATE_NAMES =
        {"birthRate", "deathRate", "diversificationRate", "reproductiveNumber", "turnover"};

    private VBox skylineEditorsBox;

    public CalibratedBirthDeathSkylineInputEditor(BeautiDoc doc) { super(doc); }
    public CalibratedBirthDeathSkylineInputEditor() { super(); }

    @Override
    public Class<?> type() { return CalibratedBirthDeathSkylineModel.class; }

    // ── Base-class hooks ────────────────────────────────────────────────────────────

    @Override
    protected String modelIdPrefix() { return "CalibratedCPP"; }

    @Override
    protected String originParamId(String partition) { return "originParam." + partition; }

    @Override
    protected void buildModelUI(Pane pane, CalibratedCoalescentPointProcess model) {
        CalibratedBirthDeathSkylineModel cpp = (CalibratedBirthDeathSkylineModel) model;

        HBox paramRow = FXUtils.newHBox();
        paramRow.setSpacing(8);
        paramRow.setPadding(new Insets(8, 0, 4, 0));
        paramRow.getChildren().add(new Label("Parameterization:"));
        ChoiceBox<String> paramChoice = new ChoiceBox<>();
        for (String[] pair : PARAM_PAIRS) paramChoice.getItems().add(pair[0]);
        int currentPairIdx = detectCurrentPairIndex(cpp);
        paramChoice.getSelectionModel().select(currentPairIdx);
        paramRow.getChildren().add(paramChoice);
        pane.getChildren().add(paramRow);

        // Ensure priors and operators exist for the currently active pair (bootstraps on first load).
        // resetEstimate=false: a plain panel load must respect each rate's current Estimate choice.
        String[] currentPair = PARAM_PAIRS[currentPairIdx];
        String partition = partitionOf(cpp);
        for (int i = 1; i <= 2; i++) activateSkylineParam(cpp, currentPair[i], partition, false);

        skylineEditorsBox = FXUtils.newVBox();
        skylineEditorsBox.setSpacing(10);
        pane.getChildren().add(skylineEditorsBox);
        refreshSkylineEditors(cpp);

        // Sampling proportion (rho) row — always shown with Estimate checkbox
        addSamplingProportionRow(pane, cpp);

        paramChoice.getSelectionModel().selectedIndexProperty().addListener((obs, oldIdx, newIdx) -> {
            if (newIdx.intValue() < 0) return;
            String[] pair = PARAM_PAIRS[newIdx.intValue()];
            try { switchParameterization(cpp, pair[1], pair[2]); } catch (Exception ex) { ex.printStackTrace(); }
            refreshSkylineEditors(cpp);
            sync();
        });
    }

    @Override
    protected void ensureOriginPriorAndOperator(String partition, RealScalarParam<?> scalar) {
        String priorId  = "originParam.prior." + partition;
        String scalerId = "originParamScaler." + partition;
        if (!doc.pluginmap.containsKey(priorId)) {
            try {
                LogNormal prior = new LogNormal();
                prior.setInputValue("M", new RealScalarParam<>(1.0, Real.INSTANCE));
                prior.setInputValue("S", new RealScalarParam<>(1.0, PositiveReal.INSTANCE));
                prior.setInputValue("param", scalar);
                prior.initAndValidate();
                pluginPut(priorId, prior);
            } catch (Exception e) { e.printStackTrace(); }
        }
        if (!doc.pluginmap.containsKey(scalerId)) {
            try {
                ScaleOperator scaler = new ScaleOperator();
                scaler.setInputValue("parameter", scalar);
                scaler.setInputValue("weight",    1.0);
                scaler.initAndValidate();
                pluginPut(scalerId, scaler);
            } catch (Exception e) { e.printStackTrace(); }
        }
    }

    // ── Parameterization ─────────────────────────────────────────────────────────

    private int detectCurrentPairIndex(CalibratedBirthDeathSkylineModel cpp) {
        List<String> active = new ArrayList<>();
        for (String name : ALL_RATE_NAMES)
            if (getSkylineInput(cpp, name).get() != null) active.add(name);
        if (active.size() == 2) {
            for (int i = 0; i < PARAM_PAIRS.length; i++) {
                if (PARAM_PAIRS[i][1].equals(active.get(0)) && PARAM_PAIRS[i][2].equals(active.get(1))) return i;
                if (PARAM_PAIRS[i][1].equals(active.get(1)) && PARAM_PAIRS[i][2].equals(active.get(0))) return i;
            }
        }
        return 0;
    }

    private void switchParameterization(CalibratedBirthDeathSkylineModel cpp, String p1, String p2) {
        String partition = partitionOf(cpp);
        for (String name : ALL_RATE_NAMES) {
            Input<SkylineParameter> in = getSkylineInput(cpp, name);
            if (in.get() != null) {
                // Always deactivate the named values StateNode (the connect conditions key off it)
                // so the deselected rate is dropped from the sampled state.
                BEASTInterface valuesNode = doc.pluginmap.get(rateValuesId(name, partition));
                if (valuesNode instanceof StateNode sn) setEstimated(sn, false);
                // The change times need the same treatment: their connectors key off this rate's
                // own changeTimes/estimate flag, so leaving it set would keep the prior, operator
                // and state node attached to a SkylineParameter the model no longer references —
                // an operator acting on a StateNode outside the posterior, which BEAST rejects.
                setEstimated(changeTimesParam(in.get()), false);
                in.setValue(null, cpp);
            }
        }
        // resetEstimate=true: (re)selecting a parameterization pair estimates both of its rates.
        activateSkylineParam(cpp, p1, partition, true);
        activateSkylineParam(cpp, p2, partition, true);
    }

    /**
     * ID of a skyline rate's values StateNode. Prefixed with {@code CBDS_} so it cannot collide with
     * standard BEAST tree-prior parameters that share these names — notably the Yule model's
     * {@code birthRate.t:$(n)} (and its Gamma prior / scaler). Without the prefix, whichever
     * tree-prior template BEAUti processes last wins the {@code birthRate.t:$(n)} slot in the
     * pluginmap, clobbering our vector with a scalar and leaving birthRate out of the state.
     */
    private static String rateValuesId(String paramName, String partition) {
        return "CBDS_" + paramName + "." + partition;
    }

    /**
     * Ensures the state node, SkylineParameter, prior and scaler for {@code paramName} exist and are
     * wired to the model. {@code resetEstimate} controls the Estimate default: pass {@code true} when
     * the user actively (re)selects this rate's parameterization pair (both rates should then be
     * sampled), and {@code false} on a plain panel reload so an existing rate keeps whatever its
     * Estimate checkbox last set. A freshly-created rate is always estimated regardless.
     */
    private void activateSkylineParam(CalibratedBirthDeathSkylineModel cpp, String paramName,
                                      String partition, boolean resetEstimate) {
        String paramId  = rateValuesId(paramName, partition);
        String spId     = paramName + "SP." + partition;
        String priorId  = paramName + ".prior." + partition;
        String scalerId = paramName + "Scaler." + partition;

        // The per-epoch rate values are the sampled StateNode. They are stored as a RealVectorParam
        // (dimension 1 for a single epoch) so that increasing the epoch count simply resizes this
        // same vector in place — the prior and operator below stay attached to the actual values.
        // Using a RealScalarParam here would break as soon as the values were promoted to a vector,
        // leaving the prior/operator/estimate flag on an orphaned node unrelated to the model.
        RealVectorParam<?> values;
        boolean created = false;
        BEASTInterface existing = doc.pluginmap.get(paramId);
        if (existing instanceof RealVectorParam<?> rvp) {
            values = rvp;
        } else {
            // The rate values must be a RealVectorParam (the template declares them as such). If a
            // RealScalarParam is found here it comes from a stale template/XML: replacing it in the
            // pluginmap via a plain addPlugin leaves BEAUti's state tracking pointing at the old
            // object, so the new vector never enters <state> and BEAST rejects the run
            // ("operator has a statenode ... missing from the state"). Unregister the stale node
            // first so the vector takes over the id cleanly.
            if (existing != null) {
                System.err.println("CalibratedCPP: replacing stale " + existing.getClass().getSimpleName()
                        + " '" + paramId + "' with a RealVectorParam — rebuild the package so its "
                        + "fxtemplates/CalibratedCPP.xml declares this rate as RealVectorParam.");
                doc.unregisterPlugin(existing);
            }
            double seed = (existing instanceof RealScalarParam<?> rsp) ? (double) rsp.get() : defaultValueFor(paramName);
            values = new RealVectorParam<>(new double[]{seed}, domainFor(paramName));
            pluginPut(paramId, values);
            created = true;
        }
        // Estimate by default when the rate is newly created or its pair is (re)selected; on a plain
        // reload leave the flag alone so an unchecked Estimate box is not silently undone.
        if (created || resetEstimate) setEstimated(values, true);

        SkylineParameter sp;
        BEASTInterface existingSP = doc.pluginmap.get(spId);
        if (existingSP instanceof SkylineParameter esp) {
            sp = esp;
        } else {
            sp = new SkylineParameter();
            pluginPut(spId, sp);
        }
        // Always bind the SkylineParameter's values to the canonical state node registered under
        // paramId. A reused SkylineParameter (persisted in the pluginmap across sessions, or loaded
        // from an older XML where the values were a promoted "<param>SP.values" vector) may still
        // reference a stale object; if that diverges from the state node the exported XML would
        // sample a parameter unrelated to the rate the model actually reads.
        if (sp.valuesInput.get() != values) {
            sp.valuesInput.setValue(values, sp);
            sp.initAndValidate();
        }
        // Only call setValue if not already set — avoids triggering spurious
        // initAndValidate() calls on the model during every init().
        if (getSkylineInput(cpp, paramName).get() != sp)
            getSkylineInput(cpp, paramName).setValue(sp, cpp);

        // Create prior and operator if not already present (only needed for non-default-pair params).
        // The prior is an IID of a LogNormal so it applies to every epoch value in the vector,
        // regardless of how many epochs the user later configures.
        if (!doc.pluginmap.containsKey(priorId)) {
            try {
                LogNormal base = new LogNormal();
                base.setInputValue("M", new RealScalarParam<>(0.0, Real.INSTANCE));
                base.setInputValue("S", new RealScalarParam<>(1.0, PositiveReal.INSTANCE));
                base.initAndValidate();

                IID prior = new IID();
                prior.setInputValue("distr", base);
                prior.setInputValue("param", values);
                prior.initAndValidate();
                pluginPut(priorId, prior);
            } catch (Exception e) { e.printStackTrace(); }
        }
        if (!doc.pluginmap.containsKey(scalerId)) {
            try {
                ScaleOperator scaler = new ScaleOperator();
                scaler.setInputValue("parameter", values);
                scaler.setInputValue("weight",    1.0);
                scaler.initAndValidate();
                pluginPut(scalerId, scaler);
            } catch (Exception e) { e.printStackTrace(); }
        }
    }

    @SuppressWarnings("unchecked")
    private Input<SkylineParameter> getSkylineInput(CalibratedBirthDeathSkylineModel cpp, String name) {
        return (Input<SkylineParameter>) switch (name) {
            case "birthRate"           -> cpp.birthRateInput;
            case "deathRate"           -> cpp.deathRateInput;
            case "diversificationRate" -> cpp.diversificationRateInput;
            case "reproductiveNumber"  -> cpp.reproductiveNumberInput;
            case "turnover"            -> cpp.turnoverInput;
            default -> throw new IllegalArgumentException("Unknown rate: " + name);
        };
    }

    // ── Sampling proportion ─────────────────────────────────────────────────────────

    private void addSamplingProportionRow(Pane parent, CalibratedBirthDeathSkylineModel cpp) {
        String partition = partitionOf(cpp);
        String spId = "samplingProportion." + partition;
        if (!(doc.pluginmap.get(spId) instanceof RealScalarParam<?> scalar)) return;

        VBox box = FXUtils.newVBox();
        box.setSpacing(4);
        box.setPadding(new Insets(6));
        box.setStyle("-fx-border-color: #b0b0b0; -fx-border-radius: 4;");

        HBox headerRow = FXUtils.newHBox();
        headerRow.setSpacing(8);
        headerRow.setAlignment(Pos.CENTER_LEFT);
        headerRow.getChildren().add(new Label("Sampling proportion (ρ)"));

        boolean estimated = scalar instanceof StateNode sn && sn.isEstimatedInput.get();
        CheckBox estimateCb = new CheckBox("Estimate");
        estimateCb.setSelected(estimated);
        headerRow.getChildren().add(estimateCb);
        box.getChildren().add(headerRow);

        double initVal = ((RealScalarParam<?>) scalar).valuesInput.get();
        TextField valueTf = compactField(initVal);
        valueTf.setPrefWidth(100);
        HBox valRow = FXUtils.newHBox();
        valRow.setSpacing(4);
        valRow.getChildren().addAll(new Label("Value:"), valueTf);
        box.getChildren().add(valRow);

        valueTf.textProperty().addListener((obs, o, nw) -> {
            try {
                double v = Double.parseDouble(nw);
                @SuppressWarnings({"unchecked","rawtypes"})
                RealScalarParam rsp = (RealScalarParam) scalar;
                rsp.valuesInput.setValue(v, rsp);
                rsp.initAndValidate();
            } catch (Exception ignored) {}
        });

        estimateCb.setOnAction(e -> {
            setEstimated(scalar, estimateCb.isSelected());
            if (estimateCb.isSelected()) ensureSamplingProportionPriorAndOperator(partition, scalar);
            sync();
        });

        parent.getChildren().add(box);
    }

    private void ensureSamplingProportionPriorAndOperator(String partition, RealScalarParam<?> scalar) {
        String priorId  = "samplingProportion.prior." + partition;
        String scalerId = "samplingProportionScaler." + partition;
        if (!doc.pluginmap.containsKey(priorId)) {
            try {
                Uniform prior = new Uniform();
                prior.setInputValue("lower", new RealScalarParam<>(0.0, Real.INSTANCE));
                prior.setInputValue("upper", new RealScalarParam<>(1.0, Real.INSTANCE));
                prior.setInputValue("param", scalar);
                prior.initAndValidate();
                pluginPut(priorId, prior);
            } catch (Exception e) { e.printStackTrace(); }
        }
        if (!doc.pluginmap.containsKey(scalerId)) {
            try {
                ScaleOperator scaler = new ScaleOperator();
                scaler.setInputValue("parameter", scalar);
                scaler.setInputValue("weight",    0.5);
                scaler.initAndValidate();
                pluginPut(scalerId, scaler);
            } catch (Exception e) { e.printStackTrace(); }
        }
    }

    // ── Skyline editors ───────────────────────────────────────────────────────────

    private void refreshSkylineEditors(CalibratedBirthDeathSkylineModel cpp) {
        skylineEditorsBox.getChildren().clear();
        for (String name : ALL_RATE_NAMES) {
            SkylineParameter sp = getSkylineInput(cpp, name).get();
            if (sp != null) skylineEditorsBox.getChildren().add(buildSkylineEditor(sp, name));
        }
    }

    private VBox buildSkylineEditor(SkylineParameter sp, String paramName) {
        VBox box = FXUtils.newVBox();
        box.setSpacing(4);
        box.setPadding(new Insets(6));
        box.setStyle("-fx-border-color: #b0b0b0; -fx-border-radius: 4;");

        // Header: parameter name + Estimate checkbox. Unchecking holds this rate fixed at its epoch
        // value(s): the values StateNode is dropped from the sampled state (the template's connectors
        // key its prior/operator/state off <values>/estimate=true), while the model still reads its
        // value as a constant. Both rates may be fixed for a fully-fixed birth-death process.
        HBox headerRow = FXUtils.newHBox();
        headerRow.setSpacing(8);
        headerRow.setAlignment(Pos.CENTER_LEFT);
        headerRow.getChildren().add(new Label(displayNameOf(paramName)));
        CheckBox estimateCb = new CheckBox("Estimate");
        estimateCb.setSelected(sp.valuesInput.get() instanceof StateNode sn && sn.isEstimatedInput.get());
        estimateCb.setTooltip(new Tooltip("Sample this rate. Uncheck to hold it fixed at its epoch value(s)."));
        estimateCb.setOnAction(e -> {
            setEstimated(sp.valuesInput.get(), estimateCb.isSelected());
            sync();
        });
        headerRow.getChildren().add(estimateCb);
        box.getChildren().add(headerRow);

        int nChanges = sp.changeTimesInput.get() == null ? 0 : sp.changeTimesInput.get().size();
        // Guard against corrupted state (e.g. from prior append-instead-of-replace bug).
        if (nChanges > 99) {
            nChanges = 0;
            sp.changeTimesInput.setValue(null, sp);
        }
        // Also sanitize the values vector — if absurdly large, shrink to match nChanges+1.
        if (sp.valuesInput.get() instanceof RealVectorParam<?> vp && vp.size() > 99
                && vp.size() != nChanges + 1) {
            double[] init = new double[nChanges + 1];
            Arrays.fill(init, 1.0);
            resetVector(vp, init);
        }

        // Mirror BDMM-Prime's ensureValuesConsistency(): proactively sync values dimension
        // to match change times count so the model is always consistent when the panel is built.
        ensureValues(sp, nChanges + 1);

        HBox epochRow = FXUtils.newHBox();
        epochRow.setSpacing(6);
        epochRow.getChildren().add(new Label("Epochs:"));
        Spinner<Integer> epochSpinner = new Spinner<>(1, 100, nChanges + 1);
        epochSpinner.setPrefWidth(70);
        epochRow.getChildren().add(epochSpinner);
        box.getChildren().add(epochRow);

        HBox flagsRow = FXUtils.newHBox();
        flagsRow.setSpacing(12);
        CheckBox agesCb = new CheckBox("Times are ages");
        agesCb.setSelected(sp.timesAreAgesInput.get());
        CheckBox relCb = new CheckBox("Relative to process length (0–1)");
        relCb.setSelected(sp.timesAreRelativeInput.get());
        relCb.setTooltip(new Tooltip("Change times as fractions of process length — must be between 0 and 1."));
        CheckBox estTimesCb = new CheckBox("Estimate change times");
        estTimesCb.setSelected(changeTimesParam(sp) instanceof StateNode sn && sn.isEstimatedInput.get());
        estTimesCb.setTooltip(new Tooltip("Sample the change times rather than holding them fixed. "
                + "Proposals come from a ChangeTimeOperator, which preserves their ordering."));
        estTimesCb.setDisable(nChanges == 0);
        Button distBtn = new Button("Distribute evenly");
        distBtn.setTooltip(new Tooltip("Space change times evenly (fractions if relative, else 1, 2, ...)."));
        flagsRow.getChildren().addAll(agesCb, relCb, estTimesCb, distBtn);
        box.getChildren().add(flagsRow);

        HBox ctRow = FXUtils.newHBox();
        ctRow.setSpacing(4);
        VBox ctBox = FXUtils.newVBox();
        ctBox.setSpacing(2);
        ctBox.getChildren().addAll(new Label("Change times:"), ctRow);
        ctBox.setVisible(nChanges > 0);
        ctBox.setManaged(nChanges > 0);
        box.getChildren().add(ctBox);

        HBox valRow = FXUtils.newHBox();
        valRow.setSpacing(4);
        VBox valBox = FXUtils.newVBox();
        valBox.setSpacing(2);
        valBox.getChildren().addAll(new Label("Epoch values:"), valRow);
        box.getChildren().add(valBox);

        rebuildSkylineRows(sp, nChanges, ctRow, valRow, agesCb, relCb);

        epochSpinner.valueProperty().addListener((obs, o, nv) -> {
            int nc = nv - 1;
            // Dropping to a single epoch detaches the change-times vector from the model, so clear
            // the estimate flag first — otherwise its prior and operator would stay connected to a
            // StateNode the model no longer reads, and BEAST would reject the run.
            if (nc == 0 && estTimesCb.isSelected()) estTimesCb.setSelected(false);
            ensureChangeTimes(sp, nc);
            ensureValues(sp, nv);
            ctBox.setVisible(nc > 0);
            ctBox.setManaged(nc > 0);
            estTimesCb.setDisable(nc == 0);
            // Adopt the flag of whatever vector now backs this rate, so the box never claims the
            // opposite of what will be written to the XML.
            estTimesCb.setSelected(changeTimesParam(sp) instanceof StateNode sn && sn.isEstimatedInput.get());
            rebuildSkylineRows(sp, nc, ctRow, valRow, agesCb, relCb);
        });
        estTimesCb.selectedProperty().addListener((obs, o, n) -> {
            // sync() re-parses the model, which runs SkylineParameter.initAndValidate() and its
            // strictly-increasing requirement. Normalise first so hand-entered times out of order
            // cannot abort the sync (and with it the connector run that attaches the prior).
            normaliseChangeTimes(sp);
            rebuildSkylineRows(sp, epochSpinner.getValue() - 1, ctRow, valRow, agesCb, relCb);
            setChangeTimesEstimated(doc, changeTimesParam(sp), sp.timesAreRelativeInput.get(), n);
            sync();
        });
        agesCb.selectedProperty().addListener((obs, o, n) -> {
            sp.timesAreAgesInput.setValue(n, sp);
            rebuildSkylineRows(sp, epochSpinner.getValue() - 1, ctRow, valRow, agesCb, relCb);
            revalidate(sp);
        });
        relCb.selectedProperty().addListener((obs, o, n) -> {
            sp.timesAreRelativeInput.setValue(n, sp);
            // The support changes with the scale (fractions vs absolute times), so refresh the
            // prior and operator bounds to match.
            if (estTimesCb.isSelected()) setChangeTimesEstimated(doc, changeTimesParam(sp), n, true);
            int nc = epochSpinner.getValue() - 1;
            ctBox.setVisible(nc > 0);
            ctBox.setManaged(nc > 0);
            rebuildSkylineRows(sp, nc, ctRow, valRow, agesCb, relCb);
            setTimeFieldValues(ctRow, evenlySpaced(nc, n));
            revalidate(sp);
        });
        distBtn.setOnAction(e -> {
            int nc = epochSpinner.getValue() - 1;
            setTimeFieldValues(ctRow, evenlySpaced(nc, relCb.isSelected()));
            revalidate(sp);
        });
        return box;
    }

    private void rebuildSkylineRows(SkylineParameter sp, int nChanges,
                                    HBox ctRow, HBox valRow, CheckBox agesCb, CheckBox relCb) {
        ctRow.getChildren().clear();
        valRow.getChildren().clear();
        boolean relative = relCb.isSelected();
        if (nChanges > 0) {
            RealVectorParam<?> ctParam = changeTimesParam(sp);
            double[] defaults = evenlySpaced(nChanges, relative);
            for (int i = 0; i < nChanges; i++) {
                double v = (ctParam != null && i < ctParam.size()) ? (double) ctParam.get(i) : defaults[i];
                ctRow.getChildren().add(new Label(agesCb.isSelected() ? ("Age " + (i + 1) + ":") : ("t" + (i + 1) + ":")));
                final int fi = i;
                TextField tf = compactField(v);
                String limits = relative ? "Relative change times must be between 0 and 1."
                                         : "Change times must not be negative.";
                tf.textProperty().addListener((obs, old, nw) -> {
                    boolean ok;
                    try { ok = setChangeTime(sp, fi, Double.parseDouble(nw)); }
                    catch (NumberFormatException e) { ok = false; }
                    // Mark rather than revert: reverting fights the user mid-keystroke, and the
                    // value is simply not committed to the parameter until it is valid.
                    tf.setStyle(ok ? "" : "-fx-border-color: #c0392b; -fx-border-width: 1;");
                    tf.setTooltip(ok ? null : new Tooltip(limits));
                });
                ctRow.getChildren().add(tf);
            }
        }
        int nEpochs = nChanges + 1;
        Tensor<?, Double> vParam = sp.valuesInput.get();
        for (int i = 0; i < nEpochs; i++) {
            double v = (vParam != null && i < vParam.size()) ? (double) vParam.get(i) : 1.0;
            valRow.getChildren().add(new Label("E" + (i + 1) + ":"));
            final int fi = i;
            TextField tf = compactField(v);
            tf.textProperty().addListener((obs, old, nw) -> {
                try { setEpochValue(sp, fi, Double.parseDouble(nw)); } catch (NumberFormatException ignored) {}
            });
            valRow.getChildren().add(tf);
        }
    }

    @SuppressWarnings("unchecked")
    private void ensureChangeTimes(SkylineParameter sp, int nChanges) {
        if (nChanges == 0) { sp.changeTimesInput.setValue(null, sp); return; }
        boolean relative = sp.timesAreRelativeInput.get();
        RealVectorParam<NonNegativeReal> ct = changeTimesParam(sp);
        if (ct == null) {
            // The vector may have been disconnected (setValue null) but still live in pluginmap —
            // reuse it to avoid a duplicate-ID dialog on the next addPlugin call.
            String spId = sp.getID();
            String ctId = spId != null ? spId + ".changeTimes" : null;
            if (ctId != null && doc.pluginmap.get(ctId) instanceof RealVectorParam<?> v)
                ct = (RealVectorParam<NonNegativeReal>) v;
            else {
                ct = new RealVectorParam<>(evenlySpaced(nChanges, relative), NonNegativeReal.INSTANCE);
                // StateNode's "estimate" input defaults to true, so a freshly created vector would
                // be sampled without the user ever ticking the box — and the unticked checkbox
                // would disagree with it. Change times are fixed unless explicitly estimated, which
                // is why BDMM-Prime's template spells out estimate="false" on its equivalents.
                setEstimated(ct, false);
                if (ctId != null) { ct.setID(ctId); doc.addPlugin(ct); }
            }
            sp.changeTimesInput.setValue(ct, sp);
        }
        double[] old = toDoubles(ct);
        if (old.length == nChanges) return;  // already correct size — skip resetVector and its cascade
        resetVector(ct, extendChangeTimes(old, nChanges, relative));
    }

    private static void setTimeFieldValues(HBox ctRow, double[] values) {
        int vi = 0;
        for (Node node : ctRow.getChildren()) {
            if (node instanceof TextField tf && vi < values.length)
                tf.setText(String.valueOf(values[vi++]));
        }
    }

    private static double[] evenlySpaced(int n, boolean relative) {
        double[] v = new double[n];
        for (int i = 0; i < n; i++) v[i] = relative ? (double)(i + 1) / (n + 1) : (double)(i + 1);
        return v;
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    private void ensureValues(SkylineParameter sp, int nEpochs) {
        Tensor<?, Double> existing = sp.valuesInput.get();
        Real dom = existing instanceof RealScalarParam<?> rsp ? rsp.domainTypeInput.get()
                 : existing instanceof RealVectorParam<?> vp  ? vp.domainTypeInput.get()
                 : Real.INSTANCE;
        if (existing instanceof RealVectorParam<?> vp) {
            double[] old = toDoubles(vp);
            // Guard: skip resetVector (and its initAndValidate) when size is already correct.
            // This prevents a rebuild cascade — resetVector fires BEAST2 change notifications
            // that can trigger another BEAUti sync(), which re-enters buildSkylineEditor.
            if (old.length == nEpochs) return;
            double[] upd = new double[nEpochs];
            for (int i = 0; i < nEpochs; i++) upd[i] = i < old.length ? old[i] : (old.length > 0 ? old[old.length - 1] : 1.0);
            resetVector(vp, upd);
        } else {
            // existing is a scalar (or null); a scalar is valid for 1 epoch — don't promote,
            // since doc.addPlugin() fires BEAUti events that can cascade into spurious syncs.
            if (nEpochs == 1) return;
            double scalar = (existing != null) ? (double) existing.get(0) : 1.0;
            double[] init = new double[nEpochs];
            for (int i = 0; i < nEpochs; i++) init[i] = scalar;
            String spId = sp.getID();
            String vpId = spId != null ? spId + ".values" : null;
            // Reuse an orphaned vector from pluginmap rather than adding a duplicate.
            RealVectorParam<?> vp;
            if (vpId != null && doc.pluginmap.get(vpId) instanceof RealVectorParam<?> v) {
                vp = v;
            } else {
                vp = new RealVectorParam(init, dom);
                if (vpId != null) { vp.setID(vpId); doc.addPlugin(vp); }
            }
            sp.valuesInput.setValue(vp, sp);
        }
    }

    /**
     * Resizes a change-time vector to {@code nChanges}, preserving existing entries and continuing
     * past the last one so the result stays strictly increasing.
     *
     * <p>The obvious implementation — filling new slots from {@link #evenlySpaced} — breaks as soon
     * as the user has typed a value larger than the defaults: extending {@code [66.0]} would give
     * {@code [66.0, 2.0]}, which {@link SkylineParameter} rejects, aborting the sync that renders
     * the panel and connects the priors.
     */
    static double[] extendChangeTimes(double[] old, int nChanges, boolean relative) {
        if (old.length == 0) return evenlySpaced(nChanges, relative);

        double[] upd = new double[nChanges];
        int keep = Math.min(old.length, nChanges);
        System.arraycopy(old, 0, upd, 0, keep);
        if (nChanges <= keep) return upd;

        double last = upd[keep - 1];
        int nNew = nChanges - keep;
        // Relative times must stay inside (0, 1), so spread the new entries over what is left of
        // the interval; absolute times continue at the existing spacing, or unit steps if there is
        // only one entry to infer it from.
        double step = relative
                ? (1.0 - last) / (nNew + 1)
                : (keep > 1 ? upd[keep - 1] - upd[keep - 2] : 1.0);
        // Guarantee a positive step even from a degenerate stored vector (equal entries, or a
        // relative time already at 1.0). Overshooting the relative range is preferable to emitting
        // duplicates: SkylineParameter then reports it clearly instead of the times silently tying.
        if (step <= 0) step = relative ? 1e-6 : 1.0;
        for (int i = keep; i < nChanges; i++) upd[i] = last + step * (i - keep + 1);
        return upd;
    }

    /**
     * Sorts the change times ascending, the ordering {@link SkylineParameter} requires. Panel input
     * is not MCMC state, so reordering it here is safe; it keeps the strictly-increasing check as a
     * diagnostic for hand-written XML rather than something a GUI edit can trip over.
     */
    private void normaliseChangeTimes(SkylineParameter sp) {
        RealVectorParam<?> ct = changeTimesParam(sp);
        if (ct == null || ct.size() < 2) return;
        double[] vals = toDoubles(ct);
        double[] sorted = vals.clone();
        Arrays.sort(sorted);
        if (!Arrays.equals(vals, sorted)) resetVector(ct, sorted);
    }

    /**
     * Re-runs validation after a panel edit. Failures are reported rather than swallowed: a throw
     * here means the panel has put the parameter into a state the model rejects, which previously
     * stayed invisible until the next sync() failed for no apparent reason.
     */
    private void revalidate(SkylineParameter sp) {
        try {
            sp.initAndValidate();
        } catch (Exception e) {
            System.err.println("CalibratedCPP: " + (sp.getID() != null ? sp.getID() : "skyline parameter")
                    + " is not in a valid state after this edit — " + e.getMessage());
        }
    }

    /**
     * Writes a change time, rejecting values the model would refuse.
     *
     * <p>Relative times are fractions of the process length and must lie in [0, 1]; note that
     * {@link #onOriginReconnected} switches every rate to relative when an origin is set, so a
     * panel showing an origin silently constrains these fields. Writing an out-of-range value would
     * make {@link SkylineParameter#initAndValidate} throw during the next sync, which aborts both
     * the panel refresh and the connector run that attaches the priors.
     */
    private boolean setChangeTime(SkylineParameter sp, int idx, double val) {
        RealVectorParam<?> ct = changeTimesParam(sp);
        if (ct == null) return false;
        if (val < 0.0) return false;
        if (sp.timesAreRelativeInput.get() && val > 1.0) return false;
        double[] vals = toDoubles(ct);
        if (idx >= vals.length) return false;
        vals[idx] = val;
        resetVector(ct, vals);
        // Widen the prior if the new time has outgrown its support; otherwise the prior becomes
        // unbuildable and silently vanishes from the Priors panel on the next sync.
        refreshChangeTimePrior(doc, ct, sp.timesAreRelativeInput.get());
        return true;
    }

    private void setEpochValue(SkylineParameter sp, int idx, double val) {
        if (!(sp.valuesInput.get() instanceof RealVectorParam<?> vp)) return;
        double[] vals = toDoubles(vp);
        if (idx < vals.length) { vals[idx] = val; resetVector(vp, vals); }
    }

    @SuppressWarnings("unchecked")
    private static RealVectorParam<NonNegativeReal> changeTimesParam(SkylineParameter sp) {
        Tensor<?, ?> t = sp.changeTimesInput.get();
        return (t instanceof RealVectorParam<?>) ? (RealVectorParam<NonNegativeReal>) t : null;
    }

    private static double[] toDoubles(RealVectorParam<?> vp) {
        double[] out = new double[vp.size()];
        for (int i = 0; i < out.length; i++) out[i] = (double) vp.get(i);
        return out;
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    private static void resetVector(RealVectorParam<?> vp, double[] vals) {
        // BEAST2's Input.setValue(List, ...) APPENDS to the existing list rather than replacing it.
        // We must mutate the existing list in-place to avoid unbounded growth.
        List<Double> existing = ((RealVectorParam<Real>) vp).valuesInput.get();
        existing.clear();
        for (double v : vals) existing.add(v);
        // initAndValidate() computes dim = max(dimensionInput, list.size()), so if dimensionInput
        // still holds a previously larger value, it expands the array back to that size — undoing
        // any shrink. Reset dimensionInput first so the new size takes effect.
        vp.dimensionInput.setValue(vals.length, vp);
        vp.initAndValidate();
    }

    private static String displayNameOf(String paramName) {
        return switch (paramName) {
            case "birthRate"           -> "Birth rate (λ)";
            case "deathRate"           -> "Death rate (μ)";
            case "diversificationRate" -> "Diversification rate (λ − μ)";
            case "turnover"            -> "Turnover (μ/λ)";
            case "reproductiveNumber"  -> "Reproductive number (λ/μ)";
            default -> paramName;
        };
    }

    private static double defaultValueFor(String paramName) {
        return switch (paramName) {
            case "birthRate"           -> 0.5;
            case "deathRate"           -> 0.25;
            case "diversificationRate" -> 0.1;
            case "turnover"            -> 0.5;
            case "reproductiveNumber"  -> 2.0;
            default -> 1.0;
        };
    }

    private static Real domainFor(String paramName) {
        return switch (paramName) {
            case "birthRate", "deathRate", "turnover", "reproductiveNumber" -> PositiveReal.INSTANCE;
            default -> Real.INSTANCE;
        };
    }
}
