package calibratedcpp.distribution;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.spec.inference.distribution.ScalarDistribution;
import beast.base.spec.type.Scalar;
import beast.base.spec.type.Simplex;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marcus Overwater
 */

@Description("Finite mixture f(x) = sum_i w_i f_i(x) over scalar component distributions.")
public class ScalarMixtureDistribution<S extends Scalar<?, T>, T> extends ScalarDistribution<S, T> {

    public Input<List<ScalarDistribution<S, T>>> distributionsInput =
            new Input<>("distribution", "Component distributions.", new ArrayList<>(), Input.Validate.REQUIRED);
    public Input<Simplex> weightsInput =
            new Input<>("weights", "Mixture weight for each component.", Input.Validate.REQUIRED);

    private List<ScalarDistribution<S, T>> distributions;
    private Simplex weights;

    @Override
    public void initAndValidate() {
        super.initAndValidate();
        refresh();
    }

    @Override
    public void refresh() {
        distributions = distributionsInput.get();
        weights = weightsInput.get();
        if (distributions.size() != weights.size())
            throw new IllegalArgumentException("ScalarMixtureDistribution: " + distributions.size() +
                    " distributions but " + weights.size() + " weights.");
        if (!weights.isValid())
            throw new IllegalArgumentException("ScalarMixtureDistribution: weights are not a valid simplex (sum = "
                    + weights.sum() + ").");
        for (ScalarDistribution<S, T> d : distributions)
            d.refresh();
    }

    /** log sum_i w_i f_i(x), by log-sum-exp over log w_i + log f_i(x). */
    public double logDensity(double x) {
        final int n = distributions.size();
        double max = Double.NEGATIVE_INFINITY;
        double[] terms = new double[n];
        for (int i = 0; i < n; i++) {
            double w = weights.get(i);
            terms[i] = w > 0.0 ? Math.log(w) + distributions.get(i).logDensity(x)
                    : Double.NEGATIVE_INFINITY;
            if (terms[i] > max) max = terms[i];
        }
        if (max == Double.NEGATIVE_INFINITY) return Double.NEGATIVE_INFINITY;
        double sum = 0.0;
        for (int i = 0; i < n; i++)
            sum += Math.exp(terms[i] - max);
        return max + Math.log(sum);
    }

    @Override
    public double density(double x) {
        double d = 0.0;
        for (int i = 0; i < distributions.size(); i++)
            d += weights.get(i) * distributions.get(i).density(x);
        return d;
    }

    @Override
    public double cumulativeProbability(double x) {
        double p = 0.0;
        for (int i = 0; i < distributions.size(); i++)
            p += weights.get(i) * distributions.get(i).cumulativeProbability(x);
        return p;
    }

    @Override
    protected double calcLogP(T value) {
        return logDensity(((Number) value).doubleValue());
    }

    @Override
    public double calculateLogP() {
        refresh();
        logP = calcLogP(param.get());
        return logP;
    }

    @Override
    public List<T> sample() {
        double u = rng.nextDouble();
        double cumulative = 0.0;
        int idx = distributions.size() - 1;   // guards against float drift in the last weight
        for (int i = 0; i < distributions.size(); i++) {
            cumulative += weights.get(i);
            if (u < cumulative) { idx = i; break; }
        }
        return distributions.get(idx).sample();
    }

    @Override
    public T getLowerBoundOfParameter() {
        double lo = Double.POSITIVE_INFINITY;
        for (ScalarDistribution<S, T> d : distributions)
            lo = Math.min(lo, ((Number) d.getLowerBoundOfParameter()).doubleValue());
        return (T) Double.valueOf(lo);
    }

    @Override
    public T getUpperBoundOfParameter() {
        double hi = Double.NEGATIVE_INFINITY;
        for (ScalarDistribution<S, T> d : distributions)
            hi = Math.max(hi, ((Number) d.getUpperBoundOfParameter()).doubleValue());
        return (T) Double.valueOf(hi);
    }

    @Override
    public Double getLower() { return ((Number) getLowerBoundOfParameter()).doubleValue(); }

    @Override
    public Double getUpper() { return ((Number) getUpperBoundOfParameter()).doubleValue(); }

    /**
     * Replaces the component list wholesale, keeping the {@code outputs} bookkeeping of the
     * discarded components consistent.
     *
     * <p>For BEAUti, which edits one long-lived mixture in place rather than building a new one
     * each time the user changes the components. {@link #initAndValidate()} must be called
     * afterwards, once the weights have been resized to match.
     */
    public void setComponents(List<ScalarDistribution<S, T>> components) {
        List<ScalarDistribution<S, T>> current = distributionsInput.get();
        for (ScalarDistribution<S, T> old : current)
            old.getOutputs().remove(this);
        current.clear();
        for (ScalarDistribution<S, T> c : components)
            distributionsInput.setValue(c, this);
    }
}