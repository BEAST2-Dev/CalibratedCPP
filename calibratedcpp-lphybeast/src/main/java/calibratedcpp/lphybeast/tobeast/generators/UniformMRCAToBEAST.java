package calibratedcpp.lphybeast.tobeast.generators;

import beast.base.core.BEASTInterface;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeInterface;
import beast.base.spec.evolution.tree.MRCAPrior;
import beast.base.spec.inference.distribution.Uniform;
import calibratedcpp.lphy.prior.UniformMRCA;
import lphybeast.BEASTContext;
import lphybeast.GeneratorToBEAST;

import static lphybeast.tobeast.TaxaUtils.getTaxonSet;

/**
 * Converts one {@link UniformMRCA} calibration into a BEAST
 * {@code MRCAPrior(monophyletic=true, distr=Uniform(lower,upper))}.
 *
 * <p>Listed in {@code LBcalibratedcppImpl.getExcludedGenerator()}, not {@code getGeneratorToBEASTs()},
 * so the converters call it directly instead of lphybeast's auto-traversal: its {@code taxonset}
 * must be the same {@link TaxonSet} instance the tree model's {@code calibrations} list holds for
 * that clade.
 */
public class UniformMRCAToBEAST implements GeneratorToBEAST<UniformMRCA, MRCAPrior> {

    @Override
    public MRCAPrior generatorToBEAST(UniformMRCA generator, BEASTInterface treeValue, BEASTContext context) {
        TaxonSet taxonSet = getTaxonSet((TreeInterface) treeValue, generator.getTaxa().value());
        return generatorToBEAST(generator, treeValue, taxonSet, context);
    }

    /** Variant that reuses an already-built TaxonSet, so callers can share it with other clades. */
    public MRCAPrior generatorToBEAST(UniformMRCA generator, BEASTInterface treeValue, TaxonSet taxonSet, BEASTContext context) {
        Uniform uniform = new Uniform();
        uniform.setInputValue("lower", context.getAsRealScalar(generator.getLower()));
        uniform.setInputValue("upper", context.getAsRealScalar(generator.getUpper()));
        uniform.initAndValidate();

        MRCAPrior mrcaPrior = new MRCAPrior();
        mrcaPrior.setInputValue("tree", treeValue);
        mrcaPrior.setInputValue("taxonset", taxonSet);
        mrcaPrior.setInputValue("monophyletic", true);
        mrcaPrior.setInputValue("distr", uniform);
        mrcaPrior.initAndValidate();

        context.addBEASTObject(mrcaPrior, generator);
        context.addExtraLoggable(mrcaPrior);
        return mrcaPrior;
    }

    @Override
    public Class<UniformMRCA> getGeneratorClass() {
        return UniformMRCA.class;
    }

    @Override
    public Class<MRCAPrior> getBEASTClass() {
        return MRCAPrior.class;
    }
}
