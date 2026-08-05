package calibratedcpp.lphy.tree;

import lphy.base.distribution.DistributionConstants;
import lphy.base.evolution.birthdeath.BirthDeathConstants;
import lphy.base.evolution.tree.TaxaConditionedTreeGenerator;
import lphy.base.evolution.tree.TimeTree;
import lphy.core.model.GenerativeDistribution;
import lphy.core.model.RandomVariable;
import lphy.core.model.Value;
import lphy.core.model.annotation.Citation;
import lphy.core.model.annotation.GeneratorInfo;
import lphy.core.model.annotation.ParameterInfo;

import java.util.Map;
import java.util.TreeMap;

@Citation(value = "Lambert, A., & Stadler, T. (2013). Birth-death models and coalescent point processes: the shape and probability of reconstructed phylogenies. Theoretical population biology, 90, 113–128. https://doi.org/10.1016/j.tpb.2013.10.002",
        title = "Birth–death models and coalescent point processes: The shape and probability of reconstructed phylogenies",
        DOI = "10.1016/j.tpb.2013.10.002", authors = {"Lambert","Stadler"}, year = 2013)
// TODO: needs validation
public class CPPTree implements GenerativeDistribution<TimeTree>{
    Value<Number> rootAge;
    Value<Number> birthRate;
    Value<Number> deathRate;
    Value<Number> diversification;
    Value<Number> turnover;
    Value<Number> rho;
    Value<Integer> n;
    Value<String[]> taxa;
    Value<Boolean> randomStemAge;
    double conditionAge;
    public final String randomStemAgeName = "randomStemAge";

    public CPPTree(@ParameterInfo(name = BirthDeathConstants.lambdaParamName, description = "per-lineage birth rate.", optional = true) Value<Number> birthRate,
                   @ParameterInfo(name = BirthDeathConstants.muParamName, description = "per-lineage death rate.", optional = true) Value<Number> deathRate,
                   @ParameterInfo(name = BirthDeathConstants.diversificationParamName, description = "diversification rate (lambda - mu), optional alternative to lambda/mu.", optional = true) Value<Number> diversification,
                   @ParameterInfo(name = BirthDeathConstants.turnoverParamName, description = "turnover (mu/lambda), optional alternative to lambda/mu.", optional = true) Value<Number> turnover,
                   @ParameterInfo(name = BirthDeathConstants.rhoParamName, description = "sampling probability") Value<Number> rho,
                   @ParameterInfo(name = TaxaConditionedTreeGenerator.taxaParamName, description = "name for passed in taxa", optional = true) Value<String[]> taxa,
                   @ParameterInfo(name = DistributionConstants.nParamName, description = "the total number of taxa.", optional = true) Value<Integer> n,
                   @ParameterInfo(name = BirthDeathConstants.rootAgeParamName, description = "the root age to be conditioned on optional.", optional = true) Value<Number> rootAge,
                   @ParameterInfo(name = randomStemAgeName, description = "the age of stem of the tree root, default has no stem", optional = true)Value<Boolean> randomStemAge) {
        int count = 0;
        if (birthRate != null) count++;
        if (deathRate != null) count++;
        if (diversification != null) count++;
        if (turnover != null) count++;

        if (count != 2) {
            throw new IllegalArgumentException(
                    "Must specify exactly two of: birthRate, deathRate, diversification, turnover."
            );
        }

        this.birthRate = birthRate;
        this.deathRate = deathRate;
        this.diversification = diversification;
        this.turnover = turnover;
        this.rho = rho;
        this.n = n;
        this.taxa = taxa;
        this.randomStemAge = randomStemAge;
        this.rootAge = rootAge;
    }

    @GeneratorInfo(name="CPP", examples = {"CPPTree.lphy"},
            description = "Generate a tree with a coalescent point process (CPP), node ages drawn i.i.d. and "
                    + "factorised. With a rootAge the tree is conditioned on that root age; with randomStemAge the "
                    + "origin is sampled from the max-of-n node-age law; with no n the number of tips is random.")
    @Override
    public RandomVariable<TimeTree> sample() {
        // Delegate to the shared uncalibrated CPP machinery. randomStemAge overrides rootAge: passing
        // rootAge=null makes the dispatch sample the stem from n·q·Q^(n-1) instead of conditioning on a root.
        Value<Number> effectiveRootAge =
                (getRandomStemAge() != null && getRandomStemAge().value()) ? null : getRootAge();
        CalibratedCPPTree delegate = new CalibratedCPPTree(
                getBirthRate(), getDeathRate(), getDiversificationRate(), getTurnover(),
                getSamplingProbability(), getN(), null, getTaxa(), null, effectiveRootAge);
        RandomVariable<TimeTree> sampled = delegate.sample();
        conditionAge = delegate.getOrigin().value();   // capture the sampled origin/stem for the converter
        return new RandomVariable<>("CPPTree", sampled.value(), this);
    }

    @Override
    public Map<String, Value> getParams() {
        Map<String, Value> map = new TreeMap<>();
        if (birthRate != null) map.put(BirthDeathConstants.lambdaParamName, birthRate);
        if (deathRate != null) map.put(BirthDeathConstants.muParamName, deathRate);
        if (diversification != null) map.put(BirthDeathConstants.diversificationParamName, diversification);
        if (turnover != null) map.put(BirthDeathConstants.turnoverParamName, turnover);
        map.put(BirthDeathConstants.rhoParamName, rho);
        if (rootAge != null) map.put(BirthDeathConstants.rootAgeParamName, rootAge);
        if (n != null) map.put(DistributionConstants.nParamName,n);
        if (taxa != null) map.put(TaxaConditionedTreeGenerator.taxaParamName,taxa);
        if (randomStemAge != null) map.put(randomStemAgeName, randomStemAge);
        return map;
    }

    public void setParam(String paramName, Value value){
        if (paramName.equals(BirthDeathConstants.lambdaParamName)) birthRate = value;
        else if (paramName.equals(BirthDeathConstants.muParamName)) deathRate = value;
        else if (paramName.equals(BirthDeathConstants.diversificationParamName)) diversification = value;
        else if (paramName.equals(BirthDeathConstants.turnoverParamName)) turnover = value;
        else if (paramName.equals(BirthDeathConstants.rhoParamName)) rho = value;
        else if (paramName.equals(BirthDeathConstants.rootAgeParamName)) rootAge = value;
        else if (paramName.equals(DistributionConstants.nParamName)) n = value;
        else if (paramName.equals(TaxaConditionedTreeGenerator.taxaParamName)) taxa = value;
        else if (paramName.equals(randomStemAgeName)) randomStemAge = value;
        else throw new IllegalArgumentException("Unknown parameter name: " + paramName);
    }

    public Value<Integer> getN(){
        return getParams().get(DistributionConstants.nParamName);
    }
    public Value<Number> getBirthRate(){
        return getParams().get(BirthDeathConstants.lambdaParamName);
    }
    public Value<Number> getDeathRate(){
        return getParams().get(BirthDeathConstants.muParamName);
    }
    public Value<Number> getDiversificationRate(){
        return getParams().get(BirthDeathConstants.diversificationParamName);
    }
    public Value<Number> getTurnover(){
        return getParams().get(BirthDeathConstants.turnoverParamName);
    }
    public Value<Number> getSamplingProbability(){
        return getParams().get(BirthDeathConstants.rhoParamName);
    }
    public Value<Number> getRootAge(){
        return getParams().get(BirthDeathConstants.rootAgeParamName);
    }
    public Value<Boolean> getRandomStemAge(){
        return getParams().get(randomStemAgeName);
    }
    public Value<String[]> getTaxa(){
        return getParams().get(TaxaConditionedTreeGenerator.taxaParamName);
    }
    public double getConditionAge(){
        return conditionAge;
    }
}
