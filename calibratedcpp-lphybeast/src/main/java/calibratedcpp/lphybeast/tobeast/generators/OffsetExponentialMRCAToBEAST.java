package calibratedcpp.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeInterface;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.evolution.tree.MRCAPrior;
import beast.base.spec.inference.distribution.Exponential;
import beast.base.spec.inference.distribution.OffsetReal;
import beast.base.spec.inference.parameter.RealScalarParam;
import calibratedcpp.lphy.prior.OffsetExponentialMRCA;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

import static lphybeast.tobeast.TaxaUtils.getTaxonSet;

/**
 * Converts one {@link OffsetExponentialMRCA} calibration into a BEAST
 * {@code MRCAPrior(monophyletic=true, distr=OffsetReal(offset, Exponential(mean)))}.
 */
public class OffsetExponentialMRCAToBEAST implements GeneratorToBEAST<OffsetExponentialMRCA, MRCAPrior> {

    @Override
    public MRCAPrior generatorToBEAST(OffsetExponentialMRCA generator, BEASTInterface treeValue, BEASTContext context) {
        TaxonSet taxonSet = getTaxonSet((TreeInterface) treeValue, generator.getTaxa().value());
        return generatorToBEAST(generator, treeValue, taxonSet, context);
    }

    /** Variant that reuses an already-built TaxonSet, so callers can share it with other clades. */
    public MRCAPrior generatorToBEAST(OffsetExponentialMRCA generator, BEASTInterface treeValue, TaxonSet taxonSet, BEASTContext context) {
        Exponential exponential = new Exponential();
        exponential.setInputValue("mean", new RealScalarParam<>(generator.getMean().value().doubleValue(), PositiveReal.INSTANCE));
        exponential.initAndValidate();

        OffsetReal offsetReal = new OffsetReal();
        offsetReal.setInputValue("offset", context.getAsRealScalar(generator.getOffset()));
        offsetReal.setInputValue("distribution", exponential);
        offsetReal.initAndValidate();

        MRCAPrior mrcaPrior = new MRCAPrior();
        mrcaPrior.setInputValue("tree", treeValue);
        mrcaPrior.setInputValue("taxonset", taxonSet);
        mrcaPrior.setInputValue("monophyletic", true);
        mrcaPrior.setInputValue("distr", offsetReal);
        mrcaPrior.initAndValidate();

        context.addBEASTObject(mrcaPrior, generator);
        context.addExtraLoggable(mrcaPrior);
        return mrcaPrior;
    }

    @Override
    public Class<OffsetExponentialMRCA> getGeneratorClass() {
        return OffsetExponentialMRCA.class;
    }

    @Override
    public Class<MRCAPrior> getBEASTClass() {
        return MRCAPrior.class;
    }
}
