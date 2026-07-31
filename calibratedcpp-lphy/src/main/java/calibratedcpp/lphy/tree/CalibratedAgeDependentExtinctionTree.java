package calibratedcpp.lphy.tree;

import calibratedcpp.lphy.prior.Calibration;
import calibratedcpp.lphy.prior.CalibrationArray;
import lphy.base.distribution.DistributionConstants;
import lphy.base.evolution.birthdeath.BirthDeathConstants;
import lphy.base.evolution.tree.TimeTree;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;
import org.apache.commons.math3.analysis.interpolation.SplineInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

import java.util.Map;

/**
 * Calibrated CPP where individuals have an age-dependent extinction hazard: lifetimes follow a
 * {@link LifetimeDistribution} and births happen at a constant rate. The node-age law Q(t) is derived
 * from the lifetime distribution via the Volterra integro-differential equation
 * <pre>  G'(t) = birthRate · ( G(t) - ∫₀ᵗ G(t-s) g(s) ds + S(t) ),   G(0) = 0  </pre>
 * with Q(t) = ρ·G(t) / (1 + ρ·G(t)). This solves the VIDE numerically (implicit-trapezoidal integrator,
 * product-integration convolution) — an implementation independent of the BEAST {@code CalibratedAgeDependentExtinctionModel},
 * so validation against it (and against the exponential-lifetime → constant-rate limit) is meaningful.
 * Only {@link #logCDF} is supplied; the generic numerical inverse-transform sampling in
 * {@link AbstractCalibratedCPPTree} handles the rest.
 *
 * <p>An explicit origin ({@code stemAge} or {@code rootAge}) is recommended: it sizes the solver grid
 * exactly. Without one, the grid is sized by a heuristic (mean lifetime and the calibration ages) and
 * ages beyond the horizon are clamped.
 */
public class CalibratedAgeDependentExtinctionTree extends AbstractCalibratedCPPTree {

    public static final String lifetimeDistName = "lifetimeDist";
    public static final String reproductiveNumberName = "reproductiveNumber";

    Value<Number> birthRate;
    Value<Number> reproductiveNumber;
    Value<LifetimeDistribution> lifetimeDist;

    // resolved per sample() (cached across samples while the inputs are unchanged)
    private double rhoVal, birthRateVal, horizon;
    private PolynomialSplineFunction gSpline;
    private LifetimeDistribution lastLifetime;
    private double cacheTarget = Double.NaN;
    private boolean cacheOriginGiven;

    public CalibratedAgeDependentExtinctionTree(
            @ParameterInfo(name = BirthDeathConstants.lambdaParamName, description = "per-lineage birth rate (alternative to reproductiveNumber).", optional = true) Value<Number> birthRate,
            @ParameterInfo(name = lifetimeDistName, description = "individual lifetime distribution (e.g. from weibullLifetime/gammaLifetime).") Value<LifetimeDistribution> lifetimeDist,
            @ParameterInfo(name = BirthDeathConstants.rhoParamName, description = "sampling probability.") Value<Number> rho,
            @ParameterInfo(name = DistributionConstants.nParamName, description = "the total number of taxa; omit for a random number of tips (uncalibrated only).", optional = true) Value<Integer> n,
            @ParameterInfo(name = calibrationsName, description = "an array of calibrations generated from a MRCA prior.", optional = true) Value<CalibrationArray> calibrations,
            @ParameterInfo(name = otherTaxaNames, description = "a string array of taxa names for non-calibrated tips.", optional = true) Value<String[]> otherNames,
            @ParameterInfo(name = stemAgeName, description = "the stem age working as condition time.", optional = true) Value<Number> stemAge,
            @ParameterInfo(name = rootAgeName, description = "the root age to condition on when no calibrations are provided.", optional = true) Value<Number> rootAge,
            @ParameterInfo(name = reproductiveNumberName, description = "reproductive number R0 = birthRate * mean(lifetime); R0 > 1 is supercritical. Alternative to birthRate.", optional = true) Value<Number> reproductiveNumber) {
        super(n, rho, calibrations, otherNames, stemAge, rootAge);
        if (lifetimeDist == null) throw new IllegalArgumentException("lifetimeDist must be provided.");
        if ((birthRate == null) == (reproductiveNumber == null)) {
            throw new IllegalArgumentException("Provide exactly one of birthRate or reproductiveNumber.");
        }
        this.birthRate = birthRate;
        this.reproductiveNumber = reproductiveNumber;
        this.lifetimeDist = lifetimeDist;
    }

    @GeneratorInfo(name = "CalibratedAgeDependentExtinctionTree", examples = {},
            description = "The Calibrated Coalescent Point Process with age-dependent individual lifetimes; the "
                    + "node-age law is derived from the lifetime distribution via a Volterra IDE. An exponential "
                    + "lifetime reduces to the constant-rate CalibratedCPP.")
    @Override
    public RandomVariable<TimeTree> sample() {
        return super.sample();
    }

    @Override
    protected void resolveRates() {
        double rho = getSamplingProb().value().doubleValue();
        LifetimeDistribution life = getLifetime().value();
        double birth = resolveBirthRate(life);
        boolean originGiven = getStemAge() != null || getRootAge() != null;

        // Grid extent. With an explicit origin the grid covers exactly [0, origin]; the sampler never
        // queries beyond it (node ages are truncated to the origin and no stem is sampled), so the
        // clamp can never trigger. Without one, the origin is sampled from Q(t)^n, so the grid must be
        // grown until Q saturates — otherwise a sampled origin past the grid would be silently clamped
        // (and the inverse-CDF search would run off to garbage), biasing the trees.
        double target;
        if (getStemAge() != null) {
            target = getStemAge().value().doubleValue();
        } else if (getRootAge() != null) {
            target = getRootAge().value().doubleValue();
        } else {
            target = 20.0 * life.mean();
            if (getCalibrations() != null) {
                for (Calibration c : getCalibrations().value().getCalibrationArray()) {
                    target = Math.max(target, c.getAge() * 2.0);
                }
            }
        }

        // cache on the inputs (the solution depends only on rho, birthRate, lifetime, extent, origin-given)
        if (gSpline != null && rho == rhoVal && birth == birthRateVal && life == lastLifetime
                && target == cacheTarget && originGiven == cacheOriginGiven) {
            return;
        }
        rhoVal = rho;
        birthRateVal = birth;
        lastLifetime = life;
        cacheTarget = target;
        cacheOriginGiven = originGiven;

        horizon = originGiven ? target * 1.001 : target;   // tiny margin so the origin is strictly interior
        solveVIDE(life, birth, gridSizeFor(horizon));

        if (!originGiven) {
            int attempts = 0;                              // grow until Q saturates -> sampled origin is inside the grid
            while (!isSaturated() && attempts++ < 25) {
                horizon *= 1.7;
                solveVIDE(life, birth, gridSizeFor(horizon));
            }
            if (!isSaturated()) {
                throw new RuntimeException("Age-dependent extinction process does not saturate "
                        + "(near-critical, subcritical, or numerically stiff): the max-of-n origin has no "
                        + "finite value. Provide stemAge or rootAge to bound the simulation.");
            }
        }
    }

    /** Q(horizon) has effectively reached 1. Returns false for NaN (a diverged/stiff VIDE), so the
     *  grow loop and the "does not saturate" guard treat that as unsaturated rather than falling through. */
    private boolean isSaturated() {
        return qAtHorizon() >= 1.0 - 1e-8;
    }

    /** Grid step ~0.01 (Richardson makes this O(h^4)); floored/capped and forced even for Richardson. */
    private int gridSizeFor(double H) {
        int N = Math.min(20000, Math.max(400, (int) Math.ceil(H * 100)));
        return (N % 2 == 0) ? N : N + 1;
    }

    private double qAtHorizon() {
        double G = gSpline.value(horizon);
        return rhoVal * G / (1.0 + rhoVal * G);
    }

    /**
     * Solve the VIDE on [0, horizon] to O(h^4): the implicit-trapezoidal scheme is run at step h (N
     * points) and step 2h (N/2 points), then Richardson-extrapolated, G4 = (4·G_h - G_2h)/3, at the
     * coarse-grid points, and a natural cubic spline is fitted through them. Mirrors the BEAST
     * CalibratedAgeDependentExtinctionModel; the direct O(N^2) product-integration convolution replaces
     * its FFT divide-and-conquer (Richardson keeps N small, so this is ample for simulation).
     */
    private void solveVIDE(LifetimeDistribution life, double birth, int N) {
        double[] gFine = computeGridG(N, life, birth);
        double[] gCoarse = computeGridG(N / 2, life, birth);

        int M = N / 2;
        double hc = horizon / M;
        double[] coarseT = new double[M + 1];
        double[] gRich = new double[M + 1];
        for (int j = 0; j <= M; j++) {
            coarseT[j] = j * hc;
            gRich[j] = (4.0 * gFine[2 * j] - gCoarse[j]) / 3.0;   // Richardson: cancels the O(h^2) term
        }
        coarseT[M] = horizon;
        gSpline = new SplineInterpolator().interpolate(coarseT, gRich);
    }

    /**
     * Implicit-trapezoidal solve of G'(t) = birth·(G - conv + S) on a uniform grid of N steps, where
     * conv is the convolution ∫₀^{t} G(t-s) g(s) ds. The convolution uses <em>product integration</em>
     * rather than sampling the density: each panel contributes G at its centre times the panel's exact
     * probability mass w[k], drawn from the (always finite, smooth) survival function. This is essential
     * for a Weibull shape &lt; 1, whose density g(s) ~ s^(k-1) diverges at s=0 — sampling it as h·g[k]
     * loses an order of accuracy (~30% error in G at a step of 0.01 for k=0.4) and forces the diverging
     * g(0) to be special-cased. The mass weights capture the singular mass exactly and restore the O(h^2)
     * convergence that the Richardson step in {@link #solveVIDE} relies on. K accumulates the past-G part
     * directly (O(N^2) overall); the G(t_{i+1}) endpoint is folded into the implicit denominator.
     */
    private double[] computeGridG(int N, LifetimeDistribution life, double birth) {
        double h = horizon / N;
        // S[j] = survival at the grid node; w[k] = probability mass of the panel centred on t_k (the
        // half-panel [0, h/2] for k=0), so the convolution never evaluates the (possibly infinite) density.
        double[] S = new double[N + 1], w = new double[N + 1];
        for (int j = 0; j <= N; j++) S[j] = life.survival(j * h);
        w[0] = 1.0 - life.survival(0.5 * h);
        for (int k = 1; k <= N; k++) w[k] = life.survival((k - 0.5) * h) - life.survival((k + 0.5) * h);

        double[] G = new double[N + 1], Gp = new double[N + 1];
        G[0] = 0.0;
        Gp[0] = birth;
        double denom = 1.0 - 0.5 * h * birth * (1.0 - w[0]);

        for (int i = 0; i < N; i++) {
            double K = 0.0;                                  // Σ_{k=1}^{i} w[k]·G[i+1-k]: past-G part of the convolution
            for (int k = 1; k <= i; k++) K += w[k] * G[i + 1 - k];
            G[i + 1] = (G[i] + 0.5 * h * Gp[i] + 0.5 * h * birth * (S[i + 1] - K)) / denom;
            double conv = w[0] * G[i + 1] + K;
            Gp[i + 1] = birth * (G[i + 1] - conv + S[i + 1]);
        }
        return G;
    }

    @Override
    protected double logCDF(double t) {
        double G = interpolateG(t);
        double fpMinusOne = rhoVal * G;                       // = F_p - 1
        double fp = 1.0 + fpMinusOne;
        if (fpMinusOne <= 0.0 || fp <= 0.0) return Double.NEGATIVE_INFINITY;
        return Math.log(fpMinusOne) - Math.log(fp);           // log( ρG / (1 + ρG) ) = log Q(t)
    }

    private double interpolateG(double t) {
        if (t <= 0.0) return 0.0;
        if (t >= horizon) return gSpline.value(horizon);   // clamp beyond the grid
        return gSpline.value(t);                            // cubic spline over the Richardson-extrapolated grid
    }
    // Clamping beyond the grid is safe by construction: for an explicit origin the sampler only needs
    // t <= origin <= horizon (a larger t only appears while invertCDF brackets, where the clamped flat
    // tail still yields the correct bracket); for a sampled origin resolveRates grows the grid until Q
    // saturates, so the clamp region is Q ~ 1.

    @Override
    protected AbstractCalibratedCPPTree newSubClade(int nTaxa, CalibrationArray subCalibrations) {
        return new CalibratedAgeDependentExtinctionTree(birthRate, lifetimeDist, getSamplingProb(),
                new Value<>("n", nTaxa), new Value<>("", subCalibrations), null, null, null, reproductiveNumber);
    }

    /** Birth rate, either given directly or derived from the reproductive number: lambda = R0 / mean(lifetime). */
    private double resolveBirthRate(LifetimeDistribution life) {
        if (getBirthRate() != null) return getBirthRate().value().doubleValue();
        return getReproductiveNumber().value().doubleValue() / life.mean();
    }

    @Override
    public Map<String, Value> getParams() {
        Map<String, Value> map = super.getParams();
        if (birthRate != null) map.put(BirthDeathConstants.lambdaParamName, birthRate);
        if (reproductiveNumber != null) map.put(reproductiveNumberName, reproductiveNumber);
        map.put(lifetimeDistName, lifetimeDist);
        return map;
    }

    @Override
    public void setParam(String paramName, Value value) {
        if (paramName.equals(BirthDeathConstants.lambdaParamName)) birthRate = value;
        else if (paramName.equals(reproductiveNumberName)) reproductiveNumber = value;
        else if (paramName.equals(lifetimeDistName)) lifetimeDist = value;
        else super.setParam(paramName, value);
    }

    public Value<Number> getBirthRate() {
        return getParams().get(BirthDeathConstants.lambdaParamName);
    }

    public Value<Number> getReproductiveNumber() {
        return getParams().get(reproductiveNumberName);
    }

    public Value<LifetimeDistribution> getLifetime() {
        return getParams().get(lifetimeDistName);
    }
}
