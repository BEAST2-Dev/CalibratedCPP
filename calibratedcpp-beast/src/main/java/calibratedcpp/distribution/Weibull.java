package calibratedcpp.distribution;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.distribution.ScalarDistribution;
import beast.base.spec.type.RealScalar;
import org.apache.commons.statistics.distribution.ContinuousDistribution;
import org.apache.commons.statistics.distribution.WeibullDistribution;

import java.util.List;

@Description("Weibull distribution: with scale θ and shape k. Density: f(x)=k/θ * (x/θ)^(k-1) * exp((-x/θ)^k)")
public class Weibull extends ScalarDistribution<RealScalar<PositiveReal>,Double> {
    public Input<RealScalar<PositiveReal>> shapeInput = new Input<>("shape", "Shape parameter (k).");
    public Input<RealScalar<PositiveReal>> scaleInput = new Input<>("scale", "Scale parameter (θ).");

    private WeibullDistribution dist = WeibullDistribution.of(1.0, 1.0);
    private ContinuousDistribution.Sampler sampler;

    public Weibull() {
        this.shapeInput = new Input<>("shape", "Shape parameter (k).");
        this.scaleInput = new Input<>("scale", "Scale parameter (k).");
        this.dist = WeibullDistribution.of(1.0, 1.0);
    }

    @Override
    public void initAndValidate() {
        refresh();
        super.initAndValidate();
    }

    @Override
    public void refresh() {
        double k = shapeInput.get().get();
        double theta = scaleInput.get().get();
        if (dist.getShape() != k || dist.getScale() != theta) {
            dist = WeibullDistribution.of(k, theta);
            sampler = null;
        }
    }

    @Override
    protected WeibullDistribution getApacheDistribution() {
        refresh();
        return dist;
    }

    @Override
    public List sample() {
        if (sampler == null)
            sampler = dist.createSampler(rng);
        return List.of(sampler.sample());
    }
}
