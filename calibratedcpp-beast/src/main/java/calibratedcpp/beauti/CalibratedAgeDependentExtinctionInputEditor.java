package calibratedcpp.beauti;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.inference.MCMC;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;
import beast.base.spec.domain.PositiveInt;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.Real;
import beast.base.spec.inference.distribution.Dirichlet;
import beast.base.spec.inference.distribution.Exponential;
import beast.base.spec.inference.distribution.Gamma;
import beast.base.spec.inference.distribution.LogNormal;
import beast.base.spec.inference.distribution.Normal;
import beast.base.spec.inference.distribution.Uniform;
import beast.base.spec.inference.operator.DeltaExchangeOperator;
import beast.base.spec.inference.operator.RealRandomWalkOperator;
import beast.base.spec.inference.operator.ScaleOperator;
import beast.base.spec.inference.parameter.IntScalarParam;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.distribution.ScalarDistribution;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.inference.parameter.SimplexParam;
import beast.base.spec.type.RealScalar;
import beastfx.app.inputeditor.BeautiConnector;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.inputeditor.BeautiSubTemplate;
import beastfx.app.util.FXUtils;
import calibratedcpp.CalibratedAgeDependentExtinctionModel;
import calibratedcpp.CalibratedCoalescentPointProcess;
import calibratedcpp.distribution.Erlang;
import calibratedcpp.distribution.ScalarMixtureDistribution;
import calibratedcpp.distribution.Weibull;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Node;
import javafx.scene.control.CheckBox;
import javafx.scene.control.ChoiceBox;
import javafx.scene.control.Label;
import javafx.scene.control.Spinner;
import javafx.scene.control.TextField;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Pane;
import javafx.scene.layout.VBox;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * BEAUti panel for {@link CalibratedAgeDependentExtinctionModel}. Inherits the
 * calibration-management popup, live preview, conditioning row, and model persistence from
 * {@link CalibratedCPPInputEditor}; this class contributes only the model-specific
 * parameter UI: birth rate, extant sampling probability, and a lifetime-distribution picker
 * (exponential / gamma / Erlang / log-normal / Weibull / finite mixture) with its associated
 * parameter editors.
 */
public class CalibratedAgeDependentExtinctionInputEditor extends CalibratedCPPInputEditor {

    private enum LifetimeKind {
        EXPONENTIAL("Exponential"), GAMMA("Gamma"), ERLANG("Erlang"), LOGNORMAL("Log-Normal"),
        WEIBULL("Weibull"), MIXTURE("Mixture");
        final String display;
        LifetimeKind(String display) { this.display = display; }
    }

    /** Families offered for a mixture component. A component may not itself be a mixture. */
    private static final LifetimeKind[] COMPONENT_KINDS = {
        LifetimeKind.EXPONENTIAL, LifetimeKind.GAMMA, LifetimeKind.ERLANG,
        LifetimeKind.LOGNORMAL, LifetimeKind.WEIBULL
    };

    private static final int MIN_MIXTURE_COMPONENTS = 2;
    private static final int MAX_MIXTURE_COMPONENTS = 8;

    /**
     * Slot for the lifetime distribution when it is the model's own, rather than a mixture
     * component. Empty so that existing single-distribution IDs are unchanged.
     */
    private static final String TOP_LEVEL_SLOT = "";

    /**
     * Scale slot used when the mixture's components share one scale parameter. Distinct from any
     * {@code "mix" + i}, so every component's scale resolves to the same ID and hence the same
     * {@code RealScalarParam}.
     */
    private static final String SHARED_SLOT = "shared";

    private enum ParamKind { POSITIVE, UNIT_INTERVAL, REAL }

    private VBox lifetimeParamsBox;

    /**
     * Un-normalised weights as typed into the mixture's weight fields. The simplex parameter
     * only ever holds the normalised values; keeping the raw entries here lets the user type a
     * weight of 2 against a weight of 1 without the field fighting the caret on every keystroke.
     */
    private double[] mixtureRawWeights = new double[0];

    public CalibratedAgeDependentExtinctionInputEditor(BeautiDoc doc) { super(doc); }
    public CalibratedAgeDependentExtinctionInputEditor() { super(); }

    @Override
    public Class<?> type() { return CalibratedAgeDependentExtinctionModel.class; }

    // ── Base-class hooks ────────────────────────────────────────────────────────────

    @Override
    protected String modelIdPrefix() { return "CalibratedAgeDependentExtinction"; }

    @Override
    protected String originParamId(String partition) { return "adOriginParam." + partition; }

    @Override
    protected void ensureOriginPriorAndOperator(String partition, RealScalarParam<?> scalar) {
        ensurePriorAndOperator(ParamKind.POSITIVE, "adOriginParam", partition, scalar);
    }

    @Override
    protected void buildModelUI(Pane pane, CalibratedCoalescentPointProcess baseModel) {
        CalibratedAgeDependentExtinctionModel model = (CalibratedAgeDependentExtinctionModel) baseModel;
        String partition = partitionOf(model);

        // Birth rate
        if (doc.pluginmap.get("adBirthRate." + partition) instanceof RealScalarParam<?> birthRate) {
            bootstrapIfEstimated("adBirthRate", partition, birthRate, ParamKind.POSITIVE);
            pane.getChildren().add(scalarRow("Birth rate (λ)", "adBirthRate", partition, birthRate, ParamKind.POSITIVE));
        }

        // Lifetime distribution picker
        HBox distRow = FXUtils.newHBox();
        distRow.setSpacing(8);
        distRow.setPadding(new Insets(8, 0, 4, 0));
        distRow.setAlignment(Pos.CENTER_LEFT);
        distRow.getChildren().add(new Label("Lifetime distribution:"));
        ChoiceBox<String> distChoice = new ChoiceBox<>();
        for (LifetimeKind k : LifetimeKind.values()) distChoice.getItems().add(k.display);
        LifetimeKind currentKind = detectCurrentLifetimeKind(model);
        distChoice.getSelectionModel().select(currentKind.ordinal());
        distRow.getChildren().add(distChoice);
        pane.getChildren().add(distRow);

        // Bootstrap: ensure the currently-selected distribution's params/prior/operator exist
        // and that the model actually points at them (handles a freshly-applied template).
        // Guarded so a failure here cannot take down the whole panel (e.g. the birth-rate row).
        try {
            Object activeDist = activateLifetimeDistribution(currentKind, TOP_LEVEL_SLOT, TOP_LEVEL_SLOT, 1.0, partition);
            if (model.lifetimeDistributionInput.get() != activeDist) setLifetimeDistribution(model, activeDist);
        } catch (RuntimeException ex) {
            ex.printStackTrace();
        }

        lifetimeParamsBox = FXUtils.newVBox();
        lifetimeParamsBox.setSpacing(10);
        pane.getChildren().add(lifetimeParamsBox);
        refreshLifetimeParamsBox(model);

        // Extant sampling probability (rho)
        if (doc.pluginmap.get("adRho." + partition) instanceof RealScalarParam<?> rho) {
            bootstrapIfEstimated("adRho", partition, rho, ParamKind.UNIT_INTERVAL);
            pane.getChildren().add(scalarRow("Extant sampling probability (ρ)", "adRho", partition, rho, ParamKind.UNIT_INTERVAL));
        }

        distChoice.getSelectionModel().selectedIndexProperty().addListener((obs, oldIdx, newIdx) -> guard(() -> {
            if (newIdx.intValue() < 0) return;
            LifetimeKind kind = LifetimeKind.values()[newIdx.intValue()];
            switchLifetimeDistribution(model, kind);
            refreshLifetimeParamsBox(model);
            sync();
        }));
    }

    // ── Scalar parameter rows (birth rate / rho / lifetime-distribution params) ────

    private VBox scalarRow(String label, String name, String partition, RealScalarParam<?> scalar, ParamKind kind) {
        return scalarRow(label, name, TOP_LEVEL_SLOT, partition, scalar, kind);
    }

    /**
     * @param name base parameter name, before the mixture slot is applied
     * @param slot {@link #TOP_LEVEL_SLOT} for a parameter the fxtemplate declares connectors for,
     *             or a mixture slot, whose parameters exist only at runtime and so must be
     *             connected to the state, prior, operators and trace log by hand
     */
    private VBox scalarRow(String label, String name, String slot, String partition,
                           RealScalarParam<?> scalar, ParamKind kind) {
        final String id = slotted(name, slot);
        VBox box = FXUtils.newVBox();
        box.setSpacing(4);
        box.setPadding(new Insets(6));
        box.setStyle("-fx-border-color: #b0b0b0; -fx-border-radius: 4;");

        HBox headerRow = FXUtils.newHBox();
        headerRow.setSpacing(8);
        headerRow.setAlignment(Pos.CENTER_LEFT);
        headerRow.getChildren().add(new Label(label));

        boolean estimated = scalar instanceof StateNode sn && sn.isEstimatedInput.get();
        CheckBox estimateCb = new CheckBox("Estimate");
        estimateCb.setSelected(estimated);
        headerRow.getChildren().add(estimateCb);
        box.getChildren().add(headerRow);

        double initVal = scalar.valuesInput.get();
        TextField valueTf = compactField(initVal);
        HBox valRow = FXUtils.newHBox();
        valRow.setSpacing(4);
        valRow.getChildren().addAll(new Label("Value:"), valueTf);
        box.getChildren().add(valRow);

        valueTf.textProperty().addListener((obs, o, nw) -> {
            try {
                double v = Double.parseDouble(nw);
                @SuppressWarnings({"unchecked", "rawtypes"})
                RealScalarParam rsp = (RealScalarParam) scalar;
                rsp.valuesInput.setValue(v, rsp);
                rsp.initAndValidate();
            } catch (Exception ignored) {}
        });

        estimateCb.setOnAction(e -> guard(() -> {
            boolean estimate = estimateCb.isSelected();
            setEstimated(scalar, estimate);
            if (estimate) ensurePriorAndOperator(kind, id, partition, scalar);
            if (!slot.isEmpty()) wireSlottedParam(kind, id, partition, scalar, estimate);
            sync();
        }));

        return box;
    }

    private void bootstrapIfEstimated(String name, String partition, RealScalarParam<?> scalar, ParamKind kind) {
        if (scalar instanceof StateNode sn && sn.isEstimatedInput.get())
            ensurePriorAndOperator(kind, name, partition, scalar);
    }

    // ── Connecting mixture-slot objects by hand ─────────────────────────────────────

    /**
     * Connects (or disconnects) a mixture component's parameter, its prior, its operator and its
     * trace-log entry.
     *
     * <p>{@code BeautiDoc.scrubAll} wires objects up purely from the {@code connect} elements
     * declared in the fxtemplate, matched by literal ID -- it does not walk the model graph. The
     * template can only name a fixed set of IDs, so it covers the model's own lifetime parameters
     * but cannot cover mixture components, whose number and families are chosen at runtime.
     * Those have to be connected here or they are never added to the state, however they are
     * marked in the UI.
     *
     * <p>Safe to call repeatedly: {@code connect} ignores duplicates and {@code disconnect}
     * ignores absent entries.
     */
    private void wireSlottedParam(ParamKind kind, String id, String partition,
                                  RealScalarParam<?> scalar, boolean estimate) {
        BEASTInterface prior = doc.pluginmap.get(id + ".prior." + partition);
        BEASTInterface op = doc.pluginmap.get(operatorId(kind, id, partition));
        // A connector lets BEAUti's scrub manage the prior connection by the estimate flag, which
        // is what survives the save-time scrub; the manual wiring below is for immediate feedback.
        ensurePriorConnector(id);
        wire(scalar, "state", "stateNode", estimate);
        wire(scalar, "tracelog", "log", estimate);
        wire(prior, "prior", "distribution", estimate);
        wireOperator(op, estimate);
        // Evict by id too, in case a scrub left an orphaned instance in the compound's list that
        // the identity-based disconnect above cannot see. The prior/operator are kept in the plugin
        // map and reused on re-estimate -- removing them there makes ensurePriorAndOperator create
        // duplicates with numbered ids.
        if (!estimate) removeFromCompoundById(id + ".prior." + partition);
    }

    /** Removes every child of the "prior" compound whose id matches, regardless of instance. */
    private void removeFromCompoundById(String priorId) {
        BEASTInterface comp = doc.pluginmap.get("prior");
        if (comp == null || !(comp.getInput("distribution").get() instanceof List<?> raw)) return;
        @SuppressWarnings("unchecked")
        List<Object> list = (List<Object>) raw;
        list.removeIf(o -> o instanceof BEASTInterface bi && priorId.equals(bi.getID()));
    }

    /**
     * Registers a BEAUti connector that ties {@code <slottedId>.prior} to the "prior" compound,
     * conditional on the parameter being estimated. This makes BEAUti's own {@code scrubAll} --
     * which runs on every save -- disconnect the prior whenever {@code estimate=false}, the same
     * way it manages the model's declared parameters. Without it, a scrub re-adds the prior for a
     * parameter that is reachable in the model, and manual disconnects here do not survive the
     * save-time scrub. Idempotent: skips if a connector for this source already exists.
     */
    private void ensurePriorConnector(String slottedId) {
        BeautiSubTemplate sub = ageDependentSubtemplate();
        if (sub == null) return;
        String src = slottedId + ".prior.t:$(n)";
        for (BeautiConnector c : sub.connectorsInput.get())
            if (src.equals(c.sourceIDInput.get()) && "distribution".equals(c.inputNameInput.get()))
                return;
        String condition = "inposterior(CalibratedAgeDependentExtinction.t:$(n)) and "
                + slottedId + ".t:$(n)/estimate=true";
        try {
            sub.connectorsInput.get().add(new BeautiConnector(src, "prior", "distribution", condition));
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private BeautiSubTemplate ageDependentSubtemplate() {
        if (doc.beautiConfig == null || doc.beautiConfig.subTemplates == null) return null;
        for (BeautiSubTemplate t : doc.beautiConfig.subTemplates)
            if (t.getMainID() != null && t.getMainID().startsWith("CalibratedAgeDependentExtinction"))
                return t;
        return null;
    }


    private void wire(BEASTInterface source, String targetId, String inputName, boolean connect) {
        if (source == null) return;
        if (connect) doc.connect(source, targetId, inputName);
        else doc.disconnect(source, targetId, inputName);
    }

    /**
     * Attaches an operator to the MCMC object directly rather than through
     * {@code BeautiDoc.connect}, which resolves its target by looking the ID up in the plugin
     * map and returns silently when that misses. Template connectors do not use the literal
     * "mcmc" either -- {@code BeautiDoc.connect(BeautiConnector, ...)} rewrites it to the
     * runnable's actual ID -- so the literal is not a reliable key here.
     */
    private void wireOperator(BEASTInterface operator, boolean connect) {
        if (operator == null || !(doc.mcmc.get() instanceof MCMC runnable)) return;
        List<Operator> operators = runnable.operatorsInput.get();
        if (connect) {
            if (!operators.contains(operator))
                runnable.operatorsInput.setValue(operator, runnable);
        } else if (operators.remove(operator)) {
            operator.getOutputs().remove(runnable);
        }
    }

    /** Operator IDs follow the naming used by the {@code ensure*} methods below. */
    private static String operatorId(ParamKind kind, String id, String partition) {
        return switch (kind) {
            case POSITIVE, UNIT_INTERVAL -> id + "Scaler." + partition;
            case REAL -> id + "RW." + partition;
        };
    }


    private void ensurePriorAndOperator(ParamKind kind, String name, String partition, RealScalarParam<?> scalar) {
        switch (kind) {
            case POSITIVE -> ensureLogNormalPriorAndScaler(name, partition, scalar);
            case UNIT_INTERVAL -> ensureUniformPriorAndScaler(name, partition, scalar);
            case REAL -> ensureNormalPriorAndRW(name, partition, scalar);
        }
    }

    private void ensureLogNormalPriorAndScaler(String name, String partition, RealScalarParam<?> scalar) {
        String priorId = name + ".prior." + partition;
        String scalerId = name + "Scaler." + partition;
        if (!doc.pluginmap.containsKey(priorId)) {
            try {
                LogNormal prior = new LogNormal();
                prior.setInputValue("M", new RealScalarParam<>(0.0, Real.INSTANCE));
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
                scaler.setInputValue("weight", 1.0);
                scaler.initAndValidate();
                pluginPut(scalerId, scaler);
            } catch (Exception e) { e.printStackTrace(); }
        }
    }

    private void ensureUniformPriorAndScaler(String name, String partition, RealScalarParam<?> scalar) {
        String priorId = name + ".prior." + partition;
        String scalerId = name + "Scaler." + partition;
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
                scaler.setInputValue("weight", 0.5);
                scaler.initAndValidate();
                pluginPut(scalerId, scaler);
            } catch (Exception e) { e.printStackTrace(); }
        }
    }

    private void ensureNormalPriorAndRW(String name, String partition, RealScalarParam<?> scalar) {
        String priorId = name + ".prior." + partition;
        String rwId = name + "RW." + partition;
        if (!doc.pluginmap.containsKey(priorId)) {
            try {
                Normal prior = new Normal();
                prior.setInputValue("mean", new RealScalarParam<>(0.0, Real.INSTANCE));
                prior.setInputValue("sigma", new RealScalarParam<>(1.0, PositiveReal.INSTANCE));
                prior.setInputValue("param", scalar);
                prior.initAndValidate();
                pluginPut(priorId, prior);
            } catch (Exception e) { e.printStackTrace(); }
        }
        if (!doc.pluginmap.containsKey(rwId)) {
            try {
                RealRandomWalkOperator rw = new RealRandomWalkOperator();
                rw.setInputValue("scalar", scalar);
                rw.setInputValue("windowSize", 1.0);
                rw.setInputValue("weight", 1.0);
                rw.initAndValidate();
                pluginPut(rwId, rw);
            } catch (Exception e) { e.printStackTrace(); }
        }
    }

    // ── Lifetime distribution ────────────────────────────────────────────────────

    private LifetimeKind detectCurrentLifetimeKind(CalibratedAgeDependentExtinctionModel model) {
        Object d = model.lifetimeDistributionInput.get();
        if (d instanceof ScalarMixtureDistribution) return LifetimeKind.MIXTURE;
        if (d instanceof Erlang) return LifetimeKind.ERLANG;
        if (d instanceof Gamma) return LifetimeKind.GAMMA;
        if (d instanceof LogNormal) return LifetimeKind.LOGNORMAL;
        if (d instanceof Weibull) return LifetimeKind.WEIBULL;
        return LifetimeKind.EXPONENTIAL;
    }

    /**
     * Creates (or reuses) the distribution of the given family, along with its parameters,
     * priors and operators.
     *
     * @param slot      distinguishes the model's own lifetime distribution ({@link #TOP_LEVEL_SLOT})
     *                  from mixture component {@code i} ({@code "mix" + i}), so that two components
     *                  of the same family get separate parameters rather than sharing one set
     * @param scaleHint starting value for whichever parameter sets the distribution's scale. Mixture
     *                  components are spread over 1, 2, 3, ... so that a fresh mixture does not start
     *                  with identical components, which is a ridge the sampler cannot leave.
     * @param scaleSlot slot for the distribution's scale parameter alone. Equal to {@code slot}
     *                  except when the mixture links its components' scales, where it is
     *                  {@link #SHARED_SLOT} so every component resolves to one scale parameter.
     */
    private Object activateLifetimeDistribution(LifetimeKind kind, String slot, String scaleSlot,
                                                double scaleHint, String partition) {
        return switch (kind) {
            case EXPONENTIAL -> activateExponential(scaleSlot, scaleHint, partition);
            case GAMMA -> activateGamma(slot, scaleSlot, scaleHint, partition);
            case ERLANG -> activateErlang(slot, scaleSlot, scaleHint, partition);
            case LOGNORMAL -> activateLogNormal(slot, scaleHint, partition);
            case WEIBULL -> activateWeibull(slot, scaleSlot, scaleHint, partition);
            case MIXTURE -> activateMixture(partition);
        };
    }

    private void switchLifetimeDistribution(CalibratedAgeDependentExtinctionModel model, LifetimeKind kind) {
        String partition = partitionOf(model);
        // Partition passed so that switching away from a mixture disconnects its components.
        deactivateCurrentLifetimeParams(model.lifetimeDistributionInput.get(), TOP_LEVEL_SLOT, partition);
        setLifetimeDistribution(model, activateLifetimeDistribution(kind, TOP_LEVEL_SLOT, TOP_LEVEL_SLOT, 1.0, partition));
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    private void setLifetimeDistribution(CalibratedAgeDependentExtinctionModel model, Object dist) {
        ((Input) model.lifetimeDistributionInput).setValue(dist, model);
    }

    private void deactivateCurrentLifetimeParams(Object dist) {
        deactivateCurrentLifetimeParams(dist, TOP_LEVEL_SLOT, null);
    }

    /**
     * Marks a distribution's parameters as not estimated, and for mixture components also
     * disconnects them, since no template connector will do it on the next scrub. Leaving them
     * connected would put parameters in the state that nothing in the posterior refers to.
     */
    private void deactivateCurrentLifetimeParams(Object dist, String slot, String partition) {
        if (dist instanceof Exponential e) {
            release("lifetimeMean", slot, partition, scalarOf(e.meanInput), ParamKind.POSITIVE);
        } else if (dist instanceof Gamma g) {
            release("lifetimeGammaShape", slot, partition, scalarOf(g.alphaInput), ParamKind.POSITIVE);
            release("lifetimeGammaScale", slot, partition, scalarOf(g.thetaInput), ParamKind.POSITIVE);
        } else if (dist instanceof Erlang e) {
            release("lifetimeErlangScale", slot, partition, scalarOf(e.scaleInput), ParamKind.POSITIVE);
        } else if (dist instanceof LogNormal ln) {
            release("lifetimeLogNormalM", slot, partition, scalarOf(ln.MParameterInput), ParamKind.REAL);
            release("lifetimeLogNormalS", slot, partition, scalarOf(ln.SParameterInput), ParamKind.POSITIVE);
        } else if (dist instanceof Weibull weibull) {
            release("lifetimeWeibullScale", slot, partition, scalarOf(weibull.scaleInput), ParamKind.POSITIVE);
            release("lifetimeWeibullShape", slot, partition, scalarOf(weibull.shapeInput), ParamKind.POSITIVE);
        } else if (dist instanceof ScalarMixtureDistribution<?, ?> mix) {
            List<? extends ScalarDistribution<?, ?>> components = mix.distributionsInput.get();
            for (int i = 0; i < components.size(); i++)
                deactivateCurrentLifetimeParams(components.get(i), "mix" + i, partition);
            SimplexParam weights = (SimplexParam) mix.weightsInput.get();
            setEstimated(weights, false);
            if (partition != null) wireMixtureWeights(partition, weights, false);
        }
    }

    private void release(String name, String slot, String partition, RealScalarParam<?> scalar, ParamKind kind) {
        setEstimated(scalar, false);
        if (!slot.isEmpty() && partition != null && scalar != null)
            wireSlottedParam(kind, slotted(name, slot), partition, scalar, false);
    }

    /** Qualifies a parameter name with its mixture slot; top-level names are left unchanged. */
    private static String slotted(String name, String slot) {
        return slot.isEmpty() ? name : name + "." + slot;
    }

    private static String distId(String family, String slot, String partition) {
        return "lifetimeDistribution." + family + (slot.isEmpty() ? "" : "." + slot) + "." + partition;
    }

    /**
     * Points {@code dist}'s named input at {@code value}, but only if it is not already there.
     * Returns whether anything changed, so the caller can skip {@code initAndValidate()} on a
     * plain rebuild (where nothing changed). Re-pointing every build is what regressed the panel:
     * it ran {@code initAndValidate()} on a live model object on every refresh.
     */
    private static boolean repoint(BEASTInterface dist, String inputName, BEASTInterface value) {
        if (dist.getInput(inputName).get() == value) return false;
        dist.setInputValue(inputName, value);
        return true;
    }

    private Exponential activateExponential(String slot, double scaleHint, String partition) {
        RealScalarParam<?> mean = ensureScalar("lifetimeMean", slot, partition, scaleHint, PositiveReal.INSTANCE, ParamKind.POSITIVE);
        String distId = distId("exponential", slot, partition);
        boolean isNew = !(doc.pluginmap.get(distId) instanceof Exponential);
        Exponential dist = isNew ? new Exponential() : (Exponential) doc.pluginmap.get(distId);
        boolean changed = repoint(dist, "mean", mean);
        if (isNew || changed) dist.initAndValidate();
        if (isNew) pluginPut(distId, dist);
        return dist;
    }

    private Gamma activateGamma(String slot, String scaleSlot, double scaleHint, String partition) {
        RealScalarParam<?> shape = ensureScalar("lifetimeGammaShape", slot, partition, 2.0, PositiveReal.INSTANCE, ParamKind.POSITIVE);
        RealScalarParam<?> scale = ensureScalar("lifetimeGammaScale", scaleSlot, partition, scaleHint, PositiveReal.INSTANCE, ParamKind.POSITIVE);
        String distId = distId("gamma", slot, partition);
        boolean isNew = !(doc.pluginmap.get(distId) instanceof Gamma);
        Gamma dist = isNew ? new Gamma() : (Gamma) doc.pluginmap.get(distId);
        boolean changed = repoint(dist, "alpha", shape) | repoint(dist, "theta", scale);
        if (isNew || changed) dist.initAndValidate();
        if (isNew) pluginPut(distId, dist);
        return dist;
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    private Erlang activateErlang(String slot, String scaleSlot, double scaleHint, String partition) {
        String shapeId = slotted("lifetimeErlangShape", slot) + "." + partition;
        IntScalarParam shape;
        if (doc.pluginmap.get(shapeId) instanceof IntScalarParam existing) {
            shape = existing;
        } else {
            shape = new IntScalarParam(2, PositiveInt.INSTANCE);
            pluginPut(shapeId, shape);
        }
        RealScalarParam<?> scale = ensureScalar("lifetimeErlangScale", scaleSlot, partition, scaleHint, PositiveReal.INSTANCE, ParamKind.POSITIVE);
        String distId = distId("erlang", slot, partition);
        boolean isNew = !(doc.pluginmap.get(distId) instanceof Erlang);
        Erlang dist = isNew ? new Erlang() : (Erlang) doc.pluginmap.get(distId);
        boolean changed = repoint(dist, "shape", shape) | repoint(dist, "scale", scale);
        if (isNew || changed) dist.initAndValidate();
        if (isNew) pluginPut(distId, dist);
        return dist;
    }

    private LogNormal activateLogNormal(String slot, double scaleHint, String partition) {
        RealScalarParam<?> m = ensureScalar("lifetimeLogNormalM", slot, partition, Math.log(scaleHint), Real.INSTANCE, ParamKind.REAL);
        RealScalarParam<?> s = ensureScalar("lifetimeLogNormalS", slot, partition, 1.0, PositiveReal.INSTANCE, ParamKind.POSITIVE);
        String distId = distId("lognormal", slot, partition);
        boolean isNew = !(doc.pluginmap.get(distId) instanceof LogNormal);
        LogNormal dist = isNew ? new LogNormal() : (LogNormal) doc.pluginmap.get(distId);
        boolean changed = repoint(dist, "M", m) | repoint(dist, "S", s);
        if (isNew || changed) dist.initAndValidate();
        if (isNew) pluginPut(distId, dist);
        return dist;
    }

    private Weibull activateWeibull(String slot, String scaleSlot, double scaleHint, String partition) {
        RealScalarParam<?> scale = ensureScalar("lifetimeWeibullScale", scaleSlot, partition, scaleHint, PositiveReal.INSTANCE, ParamKind.POSITIVE);
        RealScalarParam<?> shape = ensureScalar("lifetimeWeibullShape", slot, partition, 1.0, PositiveReal.INSTANCE, ParamKind.POSITIVE);
        String distId = distId("weibull", slot, partition);
        boolean isNew = !(doc.pluginmap.get(distId) instanceof Weibull);
        Weibull dist = isNew ? new Weibull() : (Weibull) doc.pluginmap.get(distId);
        boolean changed = repoint(dist, "scale", scale) | repoint(dist, "shape", shape);
        if (isNew || changed) dist.initAndValidate();
        if (isNew) pluginPut(distId, dist);
        return dist;
    }

    // ── Mixture ─────────────────────────────────────────────────────────────────────

    /**
     * Creates (or reuses) the partition's mixture. An existing mixture keeps the component
     * families it already has; a fresh one starts as two exponentials with different means.
     */
    private ScalarMixtureDistribution<?, ?> activateMixture(String partition) {
        String distId = distId("mixture", TOP_LEVEL_SLOT, partition);
        ScalarMixtureDistribution<?, ?> mix;
        List<LifetimeKind> kinds;

        boolean scaleLinked;
        if (doc.pluginmap.get(distId) instanceof ScalarMixtureDistribution<?, ?> existing) {
            mix = existing;
            kinds = componentKindsOf(mix);
            scaleLinked = detectScaleLinked(mix);
        } else {
            mix = new ScalarMixtureDistribution<RealScalar<PositiveReal>, Double>();
            mix.setID(distId);
            kinds = new ArrayList<>(List.of(LifetimeKind.EXPONENTIAL, LifetimeKind.EXPONENTIAL));
            scaleLinked = false;
        }

        applyMixture(mix, kinds, partition, scaleLinked);
        if (doc.pluginmap.get(distId) != mix) pluginPut(distId, mix);
        return mix;
    }

    /**
     * Rebuilds the mixture's components from {@code kinds}, resizes the weights to match, and
     * revalidates. The weights must be resized before the mixture is validated, since
     * {@code ScalarMixtureDistribution.refresh()} rejects a weight vector of the wrong length.
     *
     * @param scaleLinked when true (and {@link #canLinkScale} allows it), every component's scale
     *                    resolves to one shared parameter rather than a per-component one
     */
    @SuppressWarnings({"unchecked", "rawtypes"})
    private void applyMixture(ScalarMixtureDistribution<?, ?> mix, List<LifetimeKind> kinds,
                             String partition, boolean scaleLinked) {
        boolean linked = scaleLinked && canLinkScale(kinds);
        List components = new ArrayList<>();
        for (int i = 0; i < kinds.size(); i++) {
            String scaleSlot = linked ? SHARED_SLOT : "mix" + i;
            components.add(activateLifetimeDistribution(kinds.get(i), "mix" + i, scaleSlot, i + 1.0, partition));
        }

        // Disconnect the scale parameters this configuration no longer references, so a scale left
        // over from the other linking state does not linger in the state and posterior. Disconnect
        // only -- never write estimate flags here, since this runs on every rebuild.
        if (linked)
            for (int i = 0; i < kinds.size(); i++) disconnectScaleSlot(kinds.get(i), "mix" + i, partition);
        else
            for (LifetimeKind fam : COMPONENT_KINDS) disconnectScaleSlot(fam, SHARED_SLOT, partition);

        SimplexParam weights = ensureMixtureWeights(partition, kinds.size());

        // Raw: BEAST's own distributions disagree on the domain type parameter (Exponential is
        // NonNegativeReal, the rest PositiveReal), so no single instantiation holds them all.
        // Erased at XML level anyway -- see ScalarMixtureDistribution for what actually constrains
        // the components.
        ((ScalarMixtureDistribution) mix).setComponents(components);
        if (mix.weightsInput.get() != weights)
            mix.weightsInput.setValue(weights, mix);
        mix.initAndValidate();
    }

    private List<LifetimeKind> componentKindsOf(ScalarMixtureDistribution<?, ?> mix) {
        List<LifetimeKind> kinds = new ArrayList<>();
        for (ScalarDistribution<?, ?> component : mix.distributionsInput.get()) {
            if (component instanceof Erlang) kinds.add(LifetimeKind.ERLANG);
            else if (component instanceof Gamma) kinds.add(LifetimeKind.GAMMA);
            else if (component instanceof LogNormal) kinds.add(LifetimeKind.LOGNORMAL);
            else if (component instanceof Weibull) kinds.add(LifetimeKind.WEIBULL);
            else kinds.add(LifetimeKind.EXPONENTIAL);
        }
        if (kinds.size() < MIN_MIXTURE_COMPONENTS)
            while (kinds.size() < MIN_MIXTURE_COMPONENTS) kinds.add(LifetimeKind.EXPONENTIAL);
        return kinds;
    }

    // ── Scale linking ───────────────────────────────────────────────────────────────

    /**
     * Base name of the family's scale parameter, or {@code null} for families with no scale to
     * link. Exponential is excluded: its mean is its only parameter, so linking it would make the
     * components identical. Log-normal is excluded: neither of its parameters is unambiguously the
     * scale (its scale is {@code exp(M)}, not a parameter in its own right).
     */
    private static String scaleParamName(LifetimeKind kind) {
        return switch (kind) {
            case GAMMA -> "lifetimeGammaScale";
            case ERLANG -> "lifetimeErlangScale";
            case WEIBULL -> "lifetimeWeibullScale";
            case EXPONENTIAL, LOGNORMAL, MIXTURE -> null;
        };
    }

    /** The scale {@link Input} of a component, or {@code null} if the family has no linkable scale. */
    private static Input<?> scaleInputOf(Object d) {
        if (d instanceof Weibull w) return w.scaleInput;
        if (d instanceof Gamma g) return g.thetaInput;
        if (d instanceof Erlang e) return e.scaleInput;
        return null;
    }

    /** Scale linking is offered only when every component is the same scale-bearing family. */
    private boolean canLinkScale(List<LifetimeKind> kinds) {
        if (kinds.size() < 2 || scaleParamName(kinds.get(0)) == null) return false;
        for (LifetimeKind k : kinds)
            if (k != kinds.get(0)) return false;
        return true;
    }

    /** True when the components currently resolve their scale to a single shared parameter. */
    private boolean detectScaleLinked(ScalarMixtureDistribution<?, ?> mix) {
        List<? extends ScalarDistribution<?, ?>> comps = mix.distributionsInput.get();
        if (comps.size() < 2) return false;
        Input<?> first = scaleInputOf(comps.get(0));
        if (first == null) return false;
        Object shared = first.get();
        for (ScalarDistribution<?, ?> c : comps) {
            Input<?> s = scaleInputOf(c);
            if (s == null || s.get() != shared) return false;
        }
        return true;
    }

    /**
     * Disconnects the scale parameter held in {@code scaleSlot} from the state/prior/operators,
     * if it exists. Must NOT touch the estimate flag: this runs on every rebuild, and writing
     * estimate flags on the per-build path corrupts the user's checkbox state.
     */
    private void disconnectScaleSlot(LifetimeKind kind, String scaleSlot, String partition) {
        String name = scaleParamName(kind);
        if (name == null) return;
        String id = slotted(name, scaleSlot);
        if (doc.pluginmap.get(id + "." + partition) instanceof RealScalarParam<?> scalar)
            wireSlottedParam(ParamKind.POSITIVE, id, partition, scalar, false);
    }

    /** Creates (or reuses and resizes) the partition's mixture-weight simplex. */
    private SimplexParam ensureMixtureWeights(String partition, int n) {
        String id = "lifetimeMixtureWeights." + partition;
        SimplexParam weights;
        if (doc.pluginmap.get(id) instanceof SimplexParam existing) {
            weights = existing;
        } else {
            weights = new SimplexParam(equalWeights(n));
            pluginPut(id, weights);
            setEstimated(weights, false);
        }
        if (weights.size() != n)
            setVectorValues(weights, normalised(equalWeights(n)));
        if (weights.isEstimatedInput.get()) {
            ensureWeightsPriorAndOperator(partition, weights);
            wireMixtureWeights(partition, weights, true);
        } else {
            // Symmetric with the estimated branch: actively disconnect, so a prior/operator left
            // over from a previous estimated state -- or loaded from an edited XML -- does not
            // linger in the posterior when the weights are fixed.
            wireMixtureWeights(partition, weights, false);
        }
        return weights;
    }

    /** As {@link #wireSlottedParam}: the weights exist only at runtime, so nothing else connects them. */
    private void wireMixtureWeights(String partition, SimplexParam weights, boolean estimate) {
        BEASTInterface prior = doc.pluginmap.get("lifetimeMixtureWeights.prior." + partition);
        BEASTInterface op = doc.pluginmap.get("lifetimeMixtureWeightsExchange." + partition);
        ensurePriorConnector("lifetimeMixtureWeights");
        wire(weights, "state", "stateNode", estimate);
        wire(weights, "tracelog", "log", estimate);
        wire(prior, "prior", "distribution", estimate);
        wireOperator(op, estimate);
        // As in wireSlottedParam: also evict any orphaned instance by id.
        if (!estimate) removeFromCompoundById("lifetimeMixtureWeights.prior." + partition);
    }

    /**
     * A Dirichlet(1,...,1) prior and a delta-exchange operator for the weights. Delta exchange
     * adds to one weight and subtracts the same amount from another, which is what keeps the
     * vector on the simplex; a scale operator would not.
     */
    private void ensureWeightsPriorAndOperator(String partition, SimplexParam weights) {
        String alphaId = "lifetimeMixtureWeightsAlpha." + partition;
        String priorId = "lifetimeMixtureWeights.prior." + partition;
        String operatorId = "lifetimeMixtureWeightsExchange." + partition;
        int n = weights.size();

        RealVectorParam<PositiveReal> alpha;
        if (doc.pluginmap.get(alphaId) instanceof RealVectorParam<?> existing) {
            @SuppressWarnings("unchecked")
            RealVectorParam<PositiveReal> reused = (RealVectorParam<PositiveReal>) existing;
            alpha = reused;
        } else {
            alpha = new RealVectorParam<>(onesArray(n), PositiveReal.INSTANCE);
            pluginPut(alphaId, alpha);
            setEstimated(alpha, false);
        }
        if (alpha.size() != n) setVectorValues(alpha, onesArray(n));

        try {
            if (doc.pluginmap.get(priorId) instanceof Dirichlet existing) {
                existing.initAndValidate(); // the number of components may have changed
            } else {
                Dirichlet prior = new Dirichlet();
                prior.setInputValue("alpha", alpha);
                prior.setInputValue("param", weights);
                prior.initAndValidate();
                pluginPut(priorId, prior);
            }
        } catch (Exception e) { e.printStackTrace(); }

        if (!doc.pluginmap.containsKey(operatorId)) {
            try {
                DeltaExchangeOperator op = new DeltaExchangeOperator();
                op.setInputValue("rvparameter", weights);
                op.setInputValue("delta", 0.1);
                op.setInputValue("weight", 1.0);
                op.initAndValidate();
                pluginPut(operatorId, op);
            } catch (Exception e) { e.printStackTrace(); }
        }
    }

    private static double[] equalWeights(int n) {
        double[] w = new double[n];
        Arrays.fill(w, 1.0 / n);
        return w;
    }

    private static double[] onesArray(int n) {
        double[] a = new double[n];
        Arrays.fill(a, 1.0);
        return a;
    }

    /** Scales to sum 1, falling back to equal weights if the input carries no mass. */
    private static double[] normalised(double[] raw) {
        double sum = 0.0;
        for (double v : raw) sum += Math.max(v, 0.0);
        if (!(sum > 0.0)) return equalWeights(raw.length);
        double[] w = new double[raw.length];
        for (int i = 0; i < raw.length; i++) w[i] = Math.max(raw[i], 0.0) / sum;
        return w;
    }

    /**
     * Replaces a vector parameter's values, resizing it if needed.
     *
     * <p>{@code dimension} must be set explicitly: {@code RealVectorParam.initAndValidate} takes
     * {@code max(dimension, values.length)}, so a stale larger dimension would silently prevent
     * the parameter from shrinking.
     */
    private static void setVectorValues(RealVectorParam<?> param, double[] values) {
        param.valuesInput.get().clear();
        for (double v : values) param.valuesInput.setValue(v, param);
        param.dimensionInput.setValue(values.length, param);
        param.initAndValidate();
    }

    /** Creates (or reuses) a scalar param, marking it estimated and ensuring its prior/operator exist. */
    private RealScalarParam<?> ensureScalar(String name, String partition, double defaultVal, Real domain, ParamKind kind) {
        return ensureScalar(name, TOP_LEVEL_SLOT, partition, defaultVal, domain, kind);
    }

    private RealScalarParam<?> ensureScalar(String name, String slot, String partition,
                                            double defaultVal, Real domain, ParamKind kind) {
        String id = slotted(name, slot);
        String pluginId = id + "." + partition;
        RealScalarParam<?> scalar;
        if (doc.pluginmap.get(pluginId) instanceof RealScalarParam<?> existing) {
            scalar = existing;
        } else {
            RealScalarParam<Real> created = new RealScalarParam<>(defaultVal, domain);
            pluginPut(pluginId, created);
            // Mixture components start fixed. Selecting "Mixture" would otherwise create up to
            // 16 estimated parameters at once, each with a prior the user never asked for.
            // Top-level parameters keep their existing default of being estimated.
            setEstimated(created, slot.isEmpty());
            scalar = created;
        }
        if (scalar instanceof StateNode sn && sn.isEstimatedInput.get()) {
            ensurePriorAndOperator(kind, id, partition, scalar);
            // Connect straight away; the template has no connectors that would do it on the next
            // scrub for a mixture-slot parameter.
            if (!slot.isEmpty()) wireSlottedParam(kind, id, partition, scalar, true);
        } else if (!slot.isEmpty()) {
            // Symmetric disconnect, for the same reason as the weights: clear any prior/operator
            // wiring left from an earlier estimated state or an edited XML.
            wireSlottedParam(kind, id, partition, scalar, false);
        }
        return scalar;
    }

    private void refreshLifetimeParamsBox(CalibratedAgeDependentExtinctionModel model) {
        lifetimeParamsBox.getChildren().clear();
        String partition = partitionOf(model);
        Object d = model.lifetimeDistributionInput.get();
        try {
            if (d instanceof ScalarMixtureDistribution<?, ?> mix)
                lifetimeParamsBox.getChildren().add(mixtureBox(model, mix, partition));
            else
                lifetimeParamsBox.getChildren().addAll(lifetimeParamRows(d, TOP_LEVEL_SLOT, partition, true));
        } catch (RuntimeException ex) {
            // Surface the failure in the terminal and keep the panel usable, rather than letting
            // it break the whole editor.
            ex.printStackTrace();
            Label err = new Label("Could not render the lifetime distribution: " + ex);
            err.setWrapText(true);
            lifetimeParamsBox.getChildren().add(err);
        }
    }

    /**
     * Parameter rows for one distribution, whether it is the model's own or a mixture component.
     *
     * @param renderScale when false, the scale row is omitted -- used for a mixture component
     *                    whose scale is linked and shown once, shared, elsewhere
     */
    private List<Node> lifetimeParamRows(Object d, String slot, String partition, boolean renderScale) {
        List<Node> rows = new ArrayList<>();
        if (d instanceof Exponential e) {
            rows.add(scalarRow("Mean lifetime", "lifetimeMean", slot, partition,
                scalarOf(e.meanInput), ParamKind.POSITIVE));
        } else if (d instanceof Gamma g) {
            rows.add(scalarRow("Shape (α)", "lifetimeGammaShape", slot, partition,
                scalarOf(g.alphaInput), ParamKind.POSITIVE));
            if (renderScale)
                rows.add(scalarRow("Scale (θ)", "lifetimeGammaScale", slot, partition,
                    scalarOf(g.thetaInput), ParamKind.POSITIVE));
        } else if (d instanceof Erlang e) {
            rows.add(erlangShapeRow(e));
            if (renderScale)
                rows.add(scalarRow("Scale (θ)", "lifetimeErlangScale", slot, partition,
                    scalarOf(e.scaleInput), ParamKind.POSITIVE));
        } else if (d instanceof LogNormal ln) {
            rows.add(scalarRow("M (log mean)", "lifetimeLogNormalM", slot, partition,
                scalarOf(ln.MParameterInput), ParamKind.REAL));
            rows.add(scalarRow("S (log sd)", "lifetimeLogNormalS", slot, partition,
                scalarOf(ln.SParameterInput), ParamKind.POSITIVE));
        } else if (d instanceof Weibull weibull) {
            rows.add(scalarRow("Shape (k)", "lifetimeWeibullShape", slot, partition,
                scalarOf(weibull.shapeInput), ParamKind.POSITIVE));
            if (renderScale)
                rows.add(scalarRow("Scale (θ)", "lifetimeWeibullScale", slot, partition,
                    scalarOf(weibull.scaleInput), ParamKind.POSITIVE));
        }
        return rows;
    }

    // ── Mixture UI ──────────────────────────────────────────────────────────────────

    /**
     * Runs a UI event action, printing any exception rather than letting JavaFX swallow it (which
     * leaves a control looking dead). The action still aborts, but the stack trace reaches the
     * terminal that launched BEAUti.
     */
    private static void guard(Runnable action) {
        try {
            action.run();
        } catch (RuntimeException ex) {
            ex.printStackTrace();
        }
    }

    private VBox mixtureBox(CalibratedAgeDependentExtinctionModel model,
                            ScalarMixtureDistribution<?, ?> mix, String partition) {
        VBox box = FXUtils.newVBox();
        box.setSpacing(8);

        List<? extends ScalarDistribution<?, ?>> components = new ArrayList<>(mix.distributionsInput.get());
        SimplexParam weights = (SimplexParam) mix.weightsInput.get();

        // Seed the raw weights from the current (normalised) values whenever the mixture is redrawn.
        mixtureRawWeights = new double[components.size()];
        for (int i = 0; i < components.size(); i++)
            mixtureRawWeights[i] = weights.get(i);

        List<LifetimeKind> currentKinds = componentKindsOf(mix);
        boolean linked = detectScaleLinked(mix);

        HBox countRow = FXUtils.newHBox();
        countRow.setSpacing(6);
        countRow.setAlignment(Pos.CENTER_LEFT);
        countRow.getChildren().add(new Label("Number of components:"));
        Spinner<Integer> countSpinner =
            new Spinner<>(MIN_MIXTURE_COMPONENTS, MAX_MIXTURE_COMPONENTS, components.size());
        countSpinner.setPrefWidth(80);
        countSpinner.valueProperty().addListener((obs, o, n) -> guard(() -> {
            List<LifetimeKind> kinds = componentKindsOf(mix);
            // Release the components being dropped, or their parameters stay in the state with
            // nothing in the posterior referring to them.
            for (int i = n; i < mix.distributionsInput.get().size(); i++)
                deactivateCurrentLifetimeParams(mix.distributionsInput.get().get(i), "mix" + i, partition);
            while (kinds.size() > n) kinds.remove(kinds.size() - 1);
            while (kinds.size() < n) kinds.add(LifetimeKind.EXPONENTIAL);
            applyMixture(mix, kinds, partition, detectScaleLinked(mix));
            refreshLifetimeParamsBox(model);
            sync();
        }));
        countRow.getChildren().add(countSpinner);
        box.getChildren().add(countRow);

        // Link-scale control, offered only when every component is the same scale-bearing family.
        if (canLinkScale(currentKinds) || linked) {
            CheckBox linkScale = new CheckBox("Link scale across components");
            linkScale.setSelected(linked);
            linkScale.setDisable(!canLinkScale(currentKinds));
            linkScale.setOnAction(e -> guard(() -> {
                applyMixture(mix, componentKindsOf(mix), partition, linkScale.isSelected());
                refreshLifetimeParamsBox(model);
                sync();
            }));
            HBox linkRow = FXUtils.newHBox();
            linkRow.setSpacing(6);
            linkRow.setAlignment(Pos.CENTER_LEFT);
            linkRow.getChildren().add(linkScale);
            box.getChildren().add(linkRow);

            // The one shared scale parameter, shown once rather than per component.
            if (linked)
                box.getChildren().add(sharedScaleRow(currentKinds.get(0), components.get(0), partition));
        }

        List<Label> normalisedLabels = new ArrayList<>();
        for (int i = 0; i < components.size(); i++)
            box.getChildren().add(componentBox(model, mix, i, components.get(i), partition, !linked, normalisedLabels));

        CheckBox estimateWeights = new CheckBox("Estimate mixture weights");
        estimateWeights.setSelected(weights.isEstimatedInput.get());
        estimateWeights.setOnAction(e -> guard(() -> {
            boolean estimate = estimateWeights.isSelected();
            setEstimated(weights, estimate);
            if (estimate) ensureWeightsPriorAndOperator(partition, weights);
            wireMixtureWeights(partition, weights, estimate);
            sync();
        }));
        HBox weightsRow = FXUtils.newHBox();
        weightsRow.setSpacing(6);
        weightsRow.setAlignment(Pos.CENTER_LEFT);
        weightsRow.setPadding(new Insets(4, 0, 0, 0));
        weightsRow.getChildren().add(estimateWeights);
        box.getChildren().add(weightsRow);

        return box;
    }

    private VBox componentBox(CalibratedAgeDependentExtinctionModel model, ScalarMixtureDistribution<?, ?> mix,
                              int index, Object component, String partition, boolean renderScale,
                              List<Label> normalisedLabels) {
        VBox box = FXUtils.newVBox();
        box.setSpacing(4);
        box.setPadding(new Insets(6));
        box.setStyle("-fx-border-color: #b0b0b0; -fx-border-radius: 4;");

        HBox header = FXUtils.newHBox();
        header.setSpacing(8);
        header.setAlignment(Pos.CENTER_LEFT);
        header.getChildren().add(new Label("Component " + (index + 1) + ":"));

        ChoiceBox<String> familyChoice = new ChoiceBox<>();
        for (LifetimeKind k : COMPONENT_KINDS) familyChoice.getItems().add(k.display);
        List<LifetimeKind> currentKinds = componentKindsOf(mix);
        familyChoice.getSelectionModel().select(indexOfKind(currentKinds.get(index)));
        familyChoice.getSelectionModel().selectedIndexProperty().addListener((obs, o, n) -> guard(() -> {
            if (n.intValue() < 0) return;
            List<LifetimeKind> kinds = componentKindsOf(mix);
            deactivateCurrentLifetimeParams(mix.distributionsInput.get().get(index), "mix" + index, partition);
            kinds.set(index, COMPONENT_KINDS[n.intValue()]);
            // Keep the scales linked if the new family set still allows it; applyMixture drops
            // linking when the components are no longer homogeneous.
            applyMixture(mix, kinds, partition, detectScaleLinked(mix));
            refreshLifetimeParamsBox(model);
            sync();
        }));
        header.getChildren().add(familyChoice);

        header.getChildren().add(new Label("Weight:"));
        TextField weightField = compactField(mixtureRawWeights[index]);
        header.getChildren().add(weightField);

        SimplexParam weights = (SimplexParam) mix.weightsInput.get();
        Label normalisedLabel = new Label(formatNormalised(weights.get(index)));
        normalisedLabels.add(normalisedLabel);
        header.getChildren().add(normalisedLabel);

        weightField.textProperty().addListener((obs, o, text) -> {
            try {
                double v = Double.parseDouble(text);
                if (v < 0.0) return;
                mixtureRawWeights[index] = v;
                applyWeights(mix, normalisedLabels);
            } catch (NumberFormatException ignored) {}
        });

        box.getChildren().add(header);
        box.getChildren().addAll(lifetimeParamRows(component, "mix" + index, partition, renderScale));
        return box;
    }

    /** The single scale row shown for a scale-linked mixture, editing the shared parameter. */
    private VBox sharedScaleRow(LifetimeKind kind, Object component, String partition) {
        Input<?> scaleInput = scaleInputOf(component);
        return scalarRow("Scale (θ, shared across components)", scaleParamName(kind), SHARED_SLOT,
                partition, scalarOf(scaleInput), ParamKind.POSITIVE);
    }

    /** Pushes the normalised raw weights into the simplex parameter and updates the readouts. */
    private void applyWeights(ScalarMixtureDistribution<?, ?> mix, List<Label> normalisedLabels) {
        SimplexParam weights = (SimplexParam) mix.weightsInput.get();
        double[] normalised = normalised(mixtureRawWeights);
        setVectorValues(weights, normalised);
        for (int i = 0; i < normalisedLabels.size() && i < normalised.length; i++)
            normalisedLabels.get(i).setText(formatNormalised(normalised[i]));
    }

    private static String formatNormalised(double w) {
        return String.format("(= %.3f)", w);
    }

    private static int indexOfKind(LifetimeKind kind) {
        for (int i = 0; i < COMPONENT_KINDS.length; i++)
            if (COMPONENT_KINDS[i] == kind) return i;
        return 0;
    }

    @SuppressWarnings({"unchecked", "rawtypes"})
    private HBox erlangShapeRow(Erlang e) {
        HBox row = FXUtils.newHBox();
        row.setSpacing(6);
        row.setAlignment(Pos.CENTER_LEFT);
        row.getChildren().add(new Label("Shape (k, positive integer):"));
        IntScalarParam shape = (IntScalarParam) e.shapeInput.get();
        Spinner<Integer> spinner = new Spinner<>(1, 1000, (Integer) shape.valuesInput.get());
        spinner.setEditable(true);
        spinner.setPrefWidth(80);
        spinner.valueProperty().addListener((obs, o, n) -> {
            shape.valuesInput.setValue(n, shape);
            try { shape.initAndValidate(); } catch (Exception ignored) {}
        });
        row.getChildren().add(spinner);
        return row;
    }

    @SuppressWarnings("unchecked")
    private static RealScalarParam<?> scalarOf(Input<?> input) {
        return (input.get() instanceof RealScalarParam<?> r) ? r : null;
    }
}
