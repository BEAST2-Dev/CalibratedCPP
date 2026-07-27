package calibratedcpp;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.TreeInterface;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.inference.distribution.ScalarDistribution;
import beast.base.spec.type.RealScalar;
import java.util.Arrays;

import calibratedcpp.distribution.Erlang;
import org.hipparchus.analysis.integration.gauss.GaussIntegrator;
import org.hipparchus.analysis.integration.gauss.GaussIntegratorFactory;
import org.hipparchus.analysis.interpolation.SplineInterpolator;
import org.hipparchus.analysis.polynomials.PolynomialSplineFunction;
import org.hipparchus.complex.Complex;
import org.hipparchus.analysis.solvers.LaguerreSolver;
import org.hipparchus.exception.MathIllegalStateException;

/**
 * @author Marcus Overwater
 */
@Description("Implementation of the Calibrated Coalescent Point Process where individual lifetimes follow a " +
        "user-specified distribution and births happen at a constant rate. Erlang (integer-shape Gamma) lifetimes " +
        "are handled via closed-form partial fractions; all other distributions use a numerical Volterra IDE solver.")
public class CalibratedAgeDependentExtinctionModel extends CalibratedCoalescentPointProcess {

    public Input<ScalarDistribution<RealScalar<PositiveReal>, Double>> lifetimeDistributionInput = new Input<>("lifetimeDistribution",
            "Distribution of the lifetime of an individual.");
    public Input<RealScalar<UnitInterval>> rhoInput = new Input<>("rho", "Extant sampling probability.");
    public Input<RealScalar<PositiveReal>> birthRateInput = new Input<>("birthRate", "The birth rate.");
    public Input<RealScalar<PositiveReal>> reproductiveNumberInput = new Input<>("reproductiveNumber", "The reproductive number birthrate*mean lifetime.");
    public Input<Integer> gridSizeInput = new Input<>("gridSize",
            "Number of grid points for the numerical Volterra IDE solver (used for non-Erlang lifetime distributions).", 1000);

    protected boolean lifetimesAreErlang;
    protected boolean useNumericalSolver;
    protected double birthRate;
    protected double rho;

    // --- Erlang closed-form fields ---
    protected Erlang erlangDistribution;
    protected int n;        // Erlang shape (positive integer)
    protected double theta; // Erlang rate (1/scale)
    protected Complex[] roots;
    protected Complex[] alphas;
    protected double gammaConst;
    protected boolean erlangValid = true; // false when root-finding fails → return -Inf likelihood

    // --- Numerical VIDE fields ---
    // gSpline stores G(t) = F(t) - 1, so fp - 1 = rho*G(t) without cancellation
    protected PolynomialSplineFunction gSpline;
    protected PolynomialSplineFunction lifetimePdfSpline;  // g(s): lifetime PDF
    protected PolynomialSplineFunction survivalSpline;     // S(t) = 1 - CDF_g(t)
    private double coarsePanel;                            // first coarse-grid step (maxTime / M)
    private boolean densitySingularAt0;                   // g(0) = +inf (Weibull shape < 1)

    private static final GaussIntegratorFactory GAUSS_FACTORY = new GaussIntegratorFactory();

    @Override
    public void initAndValidate() {
        lifetimesAreErlang = lifetimeDistributionInput.get() instanceof Erlang;
        useNumericalSolver  = !lifetimesAreErlang && lifetimeDistributionInput.get() != null;
        super.initAndValidate();
    }

    @Override
    public void updateModel() {
        super.updateModel();
        preCalc();
    }

    public void preCalc() {
        boolean hasBirthRate = birthRateInput.get() != null;
        boolean hasReproductiveNumber = reproductiveNumberInput.get() != null;
        if (hasBirthRate && hasReproductiveNumber) {
            throw new IllegalArgumentException("Specify exactly one of birthRate or reproductiveNumber, not both.");
        }
        if (!hasBirthRate && !hasReproductiveNumber) {
            throw new IllegalArgumentException("Exactly one of birthRate or reproductiveNumber must be specified.");
        }
        // R0 = birthRate * mean(lifetime)  =>  birthRate = R0 / mean(lifetime).
        birthRate = hasBirthRate
                ? birthRateInput.get().get()
                : reproductiveNumberInput.get().get() / lifetimeDistributionInput.get().getMean();
        rho = rhoInput.get().get();

        if (lifetimesAreErlang) {
            preCalcErlang();
        } else if (useNumericalSolver) {
            solveVIDE();
        }
    }

    // -------------------------------------------------------------------------
    // Erlang closed-form path
    // -------------------------------------------------------------------------

    private void preCalcErlang() {
        erlangDistribution = (Erlang) lifetimeDistributionInput.get();
        double shapeParam = erlangDistribution.shapeInput.get().get();
        n = (int) Math.round(shapeParam);
        if (Math.abs(n - shapeParam) > 1e-10){
            throw new IllegalArgumentException("Shape parameter must be an integer.");
        }
        // thetaInput → input name "theta" (scale); betaInput → input name "lambda" (rate)
        theta = 1.0 / erlangDistribution.scaleInput.get().get();  // rate = 1/scale

        double[] coeffs = buildRnCoefficients(n, theta, birthRate);
        try {
            roots = new LaguerreSolver(1e-12).solveAllComplex(coeffs, 1000, 1.0);
            gammaConst = Math.pow(theta, n) / coeffs[0];
            alphas = computeResidues(roots, n, theta, coeffs);
            erlangValid = true;
        } catch (MathIllegalStateException e) {
            erlangValid = false;
        }
    }

    /**
     * Coefficients of R_n(x) in ascending power order (c[0] = constant, c[n] = leading).
     * Expanding R_n(x) = (x+theta)^n - lambda*[(x+theta)^n - theta^n]/x and collecting powers:
     * c[n] = 1; for m < n: c[m] = C(n,n-m)*theta^(n-m) - lambda*C(n,n-m-1)*theta^(n-m-1).
     * For n=2 this gives Q(x) from Lambert & Stadler (2013) Proposition 6.
     */
    private double[] buildRnCoefficients(int n, double theta, double lambda) {
        double[] c = new double[n + 1];
        c[n] = 1.0;
        for (int m = 0; m < n; m++) {
            int k = n - m;
            c[m] = binomial(n, k) * Math.pow(theta, k)
                    - lambda * binomial(n, k - 1) * Math.pow(theta, k - 1);
        }
        return c;
    }

    private Complex[] computeResidues(Complex[] roots, int n, double theta, double[] coeffs) {
        Complex[] result = new Complex[n];
        for (int j = 0; j < n; j++) {
            Complex xj = roots[j];
            Complex numerator = xj.add(theta).pow(n);
            Complex rPrime = evaluateDerivative(coeffs, xj);
            result[j] = numerator.divide(xj.multiply(rPrime));
        }
        return result;
    }

    private Complex evaluateDerivative(double[] coeffs, Complex x) {
        Complex result = Complex.ZERO;
        Complex xPow = Complex.ONE;
        for (int m = 1; m < coeffs.length; m++) {
            result = result.add(xPow.multiply(m * coeffs[m]));
            xPow = xPow.multiply(x);
        }
        return result;
    }

    /**
     * Computes scaled exponential sums to prevent overflow for large t.
     *
     * <p>Let {@code maxExp = max_j Re(roots[j]) * t}. Every term is divided by
     * {@code exp(maxExp)} before summing, so the values stay in range regardless of t.
     * The returned array contains:
     * <ul>
     *   <li>[0] maxExp</li>
     *   <li>[1] scaledF  = Re[sum_j alphas[j]          * exp(roots[j]*t - maxExp)]</li>
     *   <li>[2] scaledFP = Re[sum_j alphas[j]*roots[j]  * exp(roots[j]*t - maxExp)]</li>
     * </ul>
     * Then F_p(t) = exp(maxExp) * innerFp  where innerFp = constTerm*exp(-maxExp) + rho*scaledF,
     * and  F_p'(t) = rho * exp(maxExp) * scaledFP.
     * The exp(maxExp) factors cancel in the log density and log CDF formulas.</p>
     */
    private double[] computeScaledSums(double t) {
        double maxExp = 0.0;
        for (int j = 0; j < n; j++) {
            maxExp = Math.max(maxExp, roots[j].getReal() * t);
        }
        double scaledF = 0.0, scaledFP = 0.0;
        for (int j = 0; j < n; j++) {
            Complex expScaled = roots[j].multiply(t).subtract(maxExp).exp();
            scaledF  += alphas[j].multiply(expScaled).getReal();
            scaledFP += alphas[j].multiply(roots[j]).multiply(expScaled).getReal();
        }
        return new double[]{maxExp, scaledF, scaledFP};
    }

    // -------------------------------------------------------------------------
    // Numerical Volterra IDE solver
    // -------------------------------------------------------------------------

    /**
     * Solves the VIDE G'(t) = lambda*(G(t) - integral_0^t G(t-s)*g(s)ds + S(t))
     * on [0, maxTime] and builds cubic splines for G, the lifetime PDF, and the survival function.
     *
     * <p>Uses Richardson extrapolation over two implicit-trapezoidal solves (step h and step 2h)
     * to cancel the O(h^2) leading error, giving O(h^4) global accuracy in G(t) at ~25% extra cost.
     * The combined G values are at the coarse-grid points (N/2 + 1 points), which is sufficient
     * resolution for the cubic spline.</p>
     */
    private void solveVIDE() {
        int N = gridSizeInput.get();
        if (N % 2 != 0) N++;  // Richardson requires an even step count

        ScalarDistribution<RealScalar<PositiveReal>, Double> dist = lifetimeDistributionInput.get();
        double[] gridG_fine = computeGridG(N, dist);
        double[] gridG_coarse = computeGridG(N / 2, dist);

        // Richardson extrapolation: G_4th = (4*G_h - G_2h) / 3 at coarse-grid points.
        // Splines for the lifetime PDF and survival function are exact (no Richardson needed).
        int M = N / 2;
        double hc = maxTime / M;
        double[] coarseT  = new double[M + 1];
        double[] gridG    = new double[M + 1];
        double[] gValues  = new double[M + 1];
        double[] sValues  = new double[M + 1];

        try {
            for (int j = 0; j <= M; j++) {
                coarseT[j] = j * hc;
                gridG[j] = (4.0 * gridG_fine[2 * j] - gridG_coarse[j]) / 3.0;
                gValues[j] = dist.density(coarseT[j]);
                sValues[j] = 1.0 - dist.cumulativeProbability(coarseT[j]);
            }
            // density(0) is +inf for a Weibull shape < 1; a spline node cannot hold it. Only then, use the
            // mean density over the first panel, CDF(hc)/hc, so lifetimePdfSpline stays finite (the s=0 point
            // has measure zero in the G'(t) convolution, which is integrated over the open interval). Finite
            // densities (shape >= 1, exponential, ...) keep their exact g(0), leaving the smooth case intact.
            densitySingularAt0 = !Double.isFinite(gValues[0]);
            if (densitySingularAt0) gValues[0] = dist.cumulativeProbability(hc) / hc;
            coarsePanel = hc;
            coarseT[M] = maxTime;
        } catch (MathIllegalStateException e){
            throw new RuntimeException("Failed to evaluate lifetime distribution", e);
        }

        SplineInterpolator interp = new SplineInterpolator();
        gSpline           = interp.interpolate(coarseT, gridG);
        lifetimePdfSpline = interp.interpolate(coarseT, gValues);
        survivalSpline    = interp.interpolate(coarseT, sValues);
    }

    private double[] computeGridG(int N, ScalarDistribution<RealScalar<PositiveReal>, Double> dist) {
        double h = maxTime / N;
        double[] gridG      = new double[N + 1];
        double[] gridGPrime = new double[N + 1];
        double[] wValues    = new double[N + 1];   // per-panel probability mass (product-integration weights)
        double[] sValues    = new double[N + 1];

        // Product-integration convolution weights instead of pointwise density samples. Sampling the
        // density as h*g[k] loses an order of accuracy when g(s) ~ s^(k-1) diverges at s=0 (Weibull shape
        // < 1: ~30% error in G at a grid step of 0.01 for k=0.4) and needs the infinite g(0) special-cased.
        // w[k] is instead the exact probability mass of the panel centred on t_k, taken from the finite,
        // smooth survival function: w[k] = S((k-1/2)h) - S((k+1/2)h); w[0] the half-panel [0, h/2]. This
        // captures the singular mass at 0 exactly and restores the O(h^2) order the Richardson step needs.
        try {
            for (int j = 0; j <= N; j++) sValues[j] = 1.0 - dist.cumulativeProbability(j * h);
            wValues[0] = dist.cumulativeProbability(0.5 * h);
            for (int k = 1; k <= N; k++)
                wValues[k] = dist.cumulativeProbability((k + 0.5) * h)
                           - dist.cumulativeProbability((k - 0.5) * h);
        } catch (MathIllegalStateException e){
            throw new RuntimeException("Failed to evaluate lifetime distribution at a grid point", e);
        }

        gridG[0]      = 0.0;
        gridGPrime[0] = birthRate;
        double w0    = wValues[0];
        double denom = 1.0 - 0.5 * h * birthRate * (1.0 - w0);
        double[] K   = new double[N + 1];

        videRecurse(gridG, gridGPrime, wValues, sValues, K, w0, denom, h, 0, N);
        return gridG;
    }

    private void videRecurse(double[] gridG, double[] gridGPrime,
                             double[] wValues, double[] sValues,
                             double[] K, double w0, double denom, double h,
                             int l, int r) {
        if (r - l == 1) {
            gridG[l + 1]  = (gridG[l] + 0.5 * h * gridGPrime[l]
                    + 0.5 * h * birthRate * (sValues[l + 1] - K[l + 1])) / denom;
            double conv   = w0 * gridG[l + 1] + K[l + 1];   // w0 = mass of the s=0 half-panel
            gridGPrime[l + 1] = birthRate * (gridG[l + 1] - conv + sValues[l + 1]);
            return;
        }
        int m = (l + r) / 2;
        videRecurse(gridG, gridGPrime, wValues, sValues, K, w0, denom, h, l, m);
        addCrossContributions(gridG, wValues, K, l, m, r);
        videRecurse(gridG, gridGPrime, wValues, sValues, K, w0, denom, h, m, r);
    }

    /**
     * Adds the contribution of G[l+1..m] to K[m+1..r] via a single FFT-based linear convolution.
     *
     * <p>The required sum is K[m+1+k'] += Σ_{i=0}^{lenG-1} G[l+1+i] * w[lenG+k'-i],
     * which equals c[lenG-1+k'] where c = linearConvolve(G[l+1..m], w[1..lenG+lenK-1]) and w is the
     * product-integration mass-weight vector (already carrying the panel width, so no extra h factor).</p>
     */
    private void addCrossContributions(double[] gridG, double[] wValues, double[] K,
                                       int l, int m, int r) {
        int lenG   = m - l;
        int lenK   = r - m;
        double[] Gsub   = Arrays.copyOfRange(gridG, l + 1, m + 1);
        double[] wSlice = new double[lenG + lenK - 1];
        System.arraycopy(wValues, 1, wSlice, 0, wSlice.length);
        double[] c = linearConvolve(Gsub, wSlice);
        for (int kp = 0; kp < lenK; kp++) K[m + 1 + kp] += c[lenG - 1 + kp];
    }

    /** Zero-pad-and-FFT linear convolution of two real sequences. */
    private static double[] linearConvolve(double[] a, double[] b) {
        int n = a.length + b.length - 1;
        int m = 1;
        while (m < n) m <<= 1;
        // JTransforms uses interleaved real/imaginary in a single double[] of length 2*m
        double[] fa = new double[2 * m];
        double[] fb = new double[2 * m];
        for (int i = 0; i < a.length; i++) fa[2 * i] = a[i];
        for (int i = 0; i < b.length; i++) fb[2 * i] = b[i];
        org.jtransforms.fft.DoubleFFT_1D fft = new org.jtransforms.fft.DoubleFFT_1D(m);
        fft.complexForward(fa);
        fft.complexForward(fb);
        for (int i = 0; i < m; i++) {
            double re = fa[2*i] * fb[2*i] - fa[2*i+1] * fb[2*i+1];
            double im = fa[2*i] * fb[2*i+1] + fa[2*i+1] * fb[2*i];
            fa[2*i]   = re;
            fa[2*i+1] = im;
        }
        fft.complexInverse(fa, true); // true = scale by 1/m
        double[] result = new double[n];
        for (int i = 0; i < n; i++) result[i] = fa[2 * i];
        return result;
    }

    /**
     * Evaluates G'(t) = lambda*(G(t) - integral_0^t G(t-s)*g(s)ds + S(t)) on demand.
     *
     * <p>Rather than interpolating a stored G'(t) spline (which accumulates O(h^2) absolute
     * error and gives large relative error when G'(t) is exponentially small in the subcritical
     * case), we re-evaluate the VIDE formula directly using 32-point Gauss-Legendre quadrature
     * on the G spline. The O(h^2) equilibrium error in G(t) cancels with the corresponding
     * error in the convolution integral, leaving a remainder of order O(h^2 * S(t)) which
     * decays exponentially — giving stable relative accuracy at all t.</p>
     */
    private double evaluateGPrime(double t) {
        if (t <= 0.0) return birthRate;
        double integral;
        if (densitySingularAt0 && t > coarsePanel) {
            // g(s) has an integrable singularity at s=0 (Weibull shape < 1) that the pdf spline cannot
            // resolve. Split the convolution: on the first panel [0, δ] approximate G(t-s) ≈ G(t) and use
            // the exact mass ∫₀^δ g = 1 - S(δ) from the (smooth) survival spline; integrate the finite
            // remainder [δ, t] with Gauss-Legendre as usual. Reduces to the plain quadrature as δ → 0.
            double delta = coarsePanel;
            double singular = gSpline.value(t) * (1.0 - survivalSpline.value(delta));
            GaussIntegrator gl = GAUSS_FACTORY.legendre(32, delta, t);
            integral = singular + gl.integrate(s -> gSpline.value(t - s) * lifetimePdfSpline.value(s));
        } else {
            GaussIntegrator gl = GAUSS_FACTORY.legendre(32, 0.0, t);
            integral = gl.integrate(s -> gSpline.value(t - s) * lifetimePdfSpline.value(s));
        }
        return birthRate * (gSpline.value(t) - integral + survivalSpline.value(t));
    }

    // -------------------------------------------------------------------------
    // Density and CDF
    // -------------------------------------------------------------------------

    private static final double LOG_ARG_FLOOR = Double.MIN_NORMAL;

    private static double logFloored(double x) {
        return Math.log(Math.max(x, LOG_ARG_FLOOR));
    }

    @Override
    public double calculateLogNodeAgeDensity(double time) {
        time = Math.min(time, maxTime);
        if (lifetimesAreErlang) {
            double[] s       = computeScaledSums(time);
            double maxExp    = s[0], scaledF = s[1], scaledFP = s[2];
            double constTerm = (1.0 - rho) + rho * gammaConst;
            double innerFp = constTerm * Math.exp(-maxExp) + rho * scaledF;
            if (innerFp <= 0.0) return Double.NEGATIVE_INFINITY;
            return Math.log(rho) + logFloored(scaledFP) - maxExp - 2.0 * Math.log(innerFp);
        } else if (useNumericalSolver) {
            double g       = gSpline.value(time);
            double fp      = 1.0 + rho * g;
            if (fp <= 0.0) return Double.NEGATIVE_INFINITY; // G(t) < -1/rho: outside the model's support
            double fpPrime = rho * evaluateGPrime(time);
            return logFloored(fpPrime) - 2.0 * Math.log(fp);
        }
        return 0.0;
    }

    @Override
    public double calculateLogNodeAgeCDF(double time) {
        time = Math.min(time, maxTime);
        if (lifetimesAreErlang) {
            double[] s       = computeScaledSums(time);
            double maxExp    = s[0], scaledF = s[1];
            double constTerm = (1.0 - rho) + rho * gammaConst;
            double innerFp = constTerm * Math.exp(-maxExp) + rho * scaledF;
            if (innerFp <= 0.0) return Double.NEGATIVE_INFINITY;
            // Q(t) = 1 - exp(-maxExp)/innerFp; clamp the ratio into [0,1] so log1p stays finite.
            double oneMinusQ = Math.min(1.0, Math.exp(-maxExp) / innerFp);
            return Math.log1p(-oneMinusQ);
        } else if (useNumericalSolver) {
            double g          = gSpline.value(time); // G(t) = F(t) - 1
            double fp         = 1.0 + rho * g;
            if (fp <= 0.0) return Double.NEGATIVE_INFINITY;
            double fpMinusOne = rho * g;              // fp - 1 = rho*G, no cancellation
            return logFloored(fpMinusOne) - Math.log(fp);
        }
        return 0.0;
    }

    /**
     * Computes log(1 - Q(time)) directly rather than via the CDF.
     *
     * <p>Both paths have an exact closed form. For Erlang lifetimes the survival is the
     * {@code oneMinusQ} ratio the CDF already forms internally, kept in log space as
     * {@code -maxExp - log(innerFp)} so the linear decay in {@code maxExp} is carried
     * explicitly. For the numerical solver, F_p = 1 + rho*G and Q = rho*G / F_p, so
     * 1 - Q = 1 / F_p exactly.</p>
     *
     * <p>Deriving either from the CDF instead would lose all precision once Q rounds to 1.0.</p>
     */
    @Override
    public double calculateLogNodeAgeSurvival(double time) {
        time = Math.min(time, maxTime);
        if (lifetimesAreErlang) {
            double[] s       = computeScaledSums(time);
            double maxExp    = s[0], scaledF = s[1];
            double constTerm = (1.0 - rho) + rho * gammaConst;
            double innerFp = constTerm * Math.exp(-maxExp) + rho * scaledF;
            // Mirrors the CDF's guard: innerFp <= 0 is treated as Q = 0, hence 1 - Q = 1.
            if (innerFp <= 0.0) return 0.0;
            // min(0, .) is the log-space form of the CDF's min(1.0, .) clamp on the ratio.
            return Math.min(0.0, -maxExp - Math.log(innerFp));
        } else if (useNumericalSolver) {
            double g  = gSpline.value(time); // G(t) = F(t) - 1
            double fp = 1.0 + rho * g;
            if (fp <= 0.0) return 0.0;
            return -Math.log(fp);
        }
        // Degenerate branch: the CDF returns log Q = 0, i.e. Q = 1, so the survival is 0.
        return Double.NEGATIVE_INFINITY;
    }

    @Override
    public double calculateTreeLogLikelihood(TreeInterface tree) {
        updateModel();
        if (!erlangValid) {
            return Double.NEGATIVE_INFINITY;
        }
        return super.calculateTreeLogLikelihood(tree);
    }

    // -------------------------------------------------------------------------
    // Utilities
    // -------------------------------------------------------------------------

    private static long binomial(int n, int k) {
        if (k < 0 || k > n) return 0;
        if (k == 0 || k == n) return 1;
        k = Math.min(k, n - k);
        long result = 1;
        for (int i = 0; i < k; i++) {
            result = result * (n - i) / (i + 1);
        }
        return result;
    }
}
