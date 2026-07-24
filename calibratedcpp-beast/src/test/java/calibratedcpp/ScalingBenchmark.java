package calibratedcpp;

import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveInt;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.domain.UnitInterval;
import beast.base.spec.inference.distribution.Gamma;
import beast.base.spec.inference.parameter.IntScalarParam;
import beast.base.spec.inference.parameter.RealScalarParam;
import calibratedcpp.distribution.Erlang;
import org.openjdk.jmh.annotations.*;
import org.openjdk.jmh.results.format.ResultFormatType;
import org.openjdk.jmh.runner.Runner;
import org.openjdk.jmh.runner.options.Options;
import org.openjdk.jmh.runner.options.OptionsBuilder;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.TimeUnit;
import java.util.stream.IntStream;

/**
 * Benchmarks likelihood-evaluation time against number of leaves for three models:
 *
 *   bd        — constant-rate birth-death (closed-form, O(n))
 *   erlang    — Erlang lifetime, k=3 (partial-fraction residues, O(n·k))
 *   gammaWarm — Gamma lifetime, VIDE path, cache hit  (O(n·32) GL quadrature)
 *   gammaCold — Gamma lifetime, VIDE path, cache miss (O(N log²N) VIDE + O(n·32))
 *
 * All distributions are parameterised to the same mean lifetime 1/mu so that
 * the models are comparable at equivalent biological timescales.
 */
@State(Scope.Thread)
@BenchmarkMode(Mode.AverageTime)
@OutputTimeUnit(TimeUnit.MICROSECONDS)
@Warmup(iterations = 3, time = 1)
@Measurement(iterations = 5, time = 1)
@Fork(1)
public class ScalingBenchmark {

    @Param({"10", "25", "50", "100", "250", "500", "1000", "10000"})
    public int nLeaves;

    private static final double ORIGIN      = 10.0;
    private static final double BIRTH_RATE  = 1.0;
    private static final double DEATH_RATE  = 0.5;   // mean lifetime = 1/DEATH_RATE = 2.0
    private static final double RHO         = 1.0;
    private static final int    ERLANG_K    = 3;      // Erlang shape
    private static final double GAMMA_ALPHA = 2.5;    // non-integer → VIDE path
    private static final int    GRID_SIZE   = 500;

    private Tree tree;
    private CalibratedBirthDeathModel              bd;
    private CalibratedAgeDependentExtinctionModel erlang;
    private CalibratedAgeDependentExtinctionModel erlang10;
    private CalibratedAgeDependentExtinctionModel gammaModel;

    // Mutable handle on birthRate for the cold VIDE benchmark
    private RealScalarParam<PositiveReal> gammaLambda;
    private boolean toggle = false;

    @Setup(Level.Trial)
    public void setup() {
        // Root at 80% of origin so all node ages are within the VIDE grid
        tree = makeTree(nLeaves, ORIGIN * 0.8);

        // ── Constant-rate birth-death ─────────────────────────────────────────
        bd = new CalibratedBirthDeathModel();
        bd.initByName(
                "tree",      tree,
                "origin",    new RealScalarParam<>(ORIGIN,      PositiveReal.INSTANCE),
                "birthRate", new RealScalarParam<>(BIRTH_RATE,  PositiveReal.INSTANCE),
                "deathRate", new RealScalarParam<>(DEATH_RATE,  NonNegativeReal.INSTANCE),
                "rho",       new RealScalarParam<>(RHO,         UnitInterval.INSTANCE));

        // ── Erlang(k=3) closed-form ───────────────────────────────────────────
        // scale = mean/k = (1/DEATH_RATE)/ERLANG_K
        double erlangScale = 1.0 / (DEATH_RATE * ERLANG_K);
        Erlang erlangDist = new Erlang();
        erlangDist.initByName(
                "shape", new IntScalarParam<>(ERLANG_K,    PositiveInt.INSTANCE),
                "scale", new RealScalarParam<>(erlangScale, PositiveReal.INSTANCE));

        erlang = new CalibratedAgeDependentExtinctionModel();
        erlang.initByName(
                "tree",                 tree,
                "origin",               new RealScalarParam<>(ORIGIN,      PositiveReal.INSTANCE),
                "birthRate",            new RealScalarParam<>(BIRTH_RATE,  PositiveReal.INSTANCE),
                "rho",                  new RealScalarParam<>(RHO,         UnitInterval.INSTANCE),
                "lifetimeDistribution", erlangDist);

        Erlang erlangDist10 = new Erlang();
        erlangDist10.initByName(
                "shape", new IntScalarParam<>(10, PositiveInt.INSTANCE),
                "scale", new RealScalarParam<>(erlangScale, PositiveReal.INSTANCE)
        );
        erlang10 = new CalibratedAgeDependentExtinctionModel();
        erlang10.initByName(
                "tree", tree,
                "origin",               new RealScalarParam<>(ORIGIN,      PositiveReal.INSTANCE),
                "birthRate",            new RealScalarParam<>(BIRTH_RATE,  PositiveReal.INSTANCE),
                "rho",                  new RealScalarParam<>(RHO,         UnitInterval.INSTANCE),
                "lifetimeDistribution", erlangDist10
        );

        // ── Gamma(alpha=2.5) VIDE path ────────────────────────────────────────
        // scale = mean/alpha = (1/DEATH_RATE)/GAMMA_ALPHA  (same mean as Erlang)
        double gammaScale = 1.0 / (DEATH_RATE * GAMMA_ALPHA);
        Gamma gammaDist = new Gamma();
        gammaDist.initByName(
                "alpha", new RealScalarParam<>(GAMMA_ALPHA, PositiveReal.INSTANCE),
                "theta", new RealScalarParam<>(gammaScale,  PositiveReal.INSTANCE));

        gammaLambda = new RealScalarParam<>(BIRTH_RATE, PositiveReal.INSTANCE);
        gammaModel = new CalibratedAgeDependentExtinctionModel();
        gammaModel.initByName(
                "tree",                 tree,
                "origin",               new RealScalarParam<>(ORIGIN, PositiveReal.INSTANCE),
                "birthRate",            gammaLambda,
                "rho",                  new RealScalarParam<>(RHO,    UnitInterval.INSTANCE),
                "lifetimeDistribution", gammaDist,
                "gridSize",             GRID_SIZE);
    }

    // ── Benchmarks ────────────────────────────────────────────────────────────

    /** Constant-rate BD: analytical F(t), O(n) tree traversal. */
    @Benchmark
    public double bd() {
        return bd.calculateTreeLogLikelihood(tree);
    }

    /**
     * Erlang closed-form: Laguerre root-finding + residue computation each call,
     * then O(n·k) complex-exponential evaluations across the tree.
     */
    @Benchmark
    public double erlang() {
        return erlang.calculateTreeLogLikelihood(tree);
    }

    /**
     * Erlang closed-form with shape parameter 10
     */
    @Benchmark
    public double erlang10() {
        return erlang10.calculateTreeLogLikelihood(tree);
    }

    /**
     * Gamma VIDE cold (cache miss): toggles birth rate by 1 ULP each call to
     * force a full O(N log²N) VIDE re-solve followed by the O(n·32) tree
     * traversal. Represents an MCMC move that changes birth rate or lifetime
     * distribution parameters.
     */
    @Benchmark
    public double gammaCold() {
        toggle = !toggle;
        gammaLambda.set(toggle ? BIRTH_RATE : BIRTH_RATE * (1.0 + 1e-7));
        return gammaModel.calculateTreeLogLikelihood(tree);
    }

    // ── Tree generation ───────────────────────────────────────────────────────

    /**
     * Generates a random binary tree with {@code n} leaves via a coalescent
     * simulation. Coalescence times are evenly spaced in (0, maxAge]; random
     * pairs of lineages are merged at each step.
     */
    private static Tree makeTree(int n, double maxAge) {
        Random rng = new Random(12345); // fixed seed → reproducible topology

        List<String> names   = new ArrayList<>();
        List<Double> heights = new ArrayList<>();
        for (int i = 0; i < n; i++) { names.add("leaf_" + i); heights.add(0.0); }

        for (int c = 0; c < n - 1; c++) {
            double t = maxAge * (c + 1.0) / n;   // evenly-spaced coalescence times

            int a = rng.nextInt(names.size());
            int b; do { b = rng.nextInt(names.size()); } while (b == a);
            if (b < a) { int tmp = a; a = b; b = tmp; } // ensure b > a for clean removal

            String merged = "(" +
                    names.get(a) + ":" + (t - heights.get(a)) + "," +
                    names.get(b) + ":" + (t - heights.get(b)) + ")";

            names.remove(b);    heights.remove(b);
            names.set(a, merged); heights.set(a, t);
        }

        Tree tree = new TreeParser();
        tree.initByName("newick", names.get(0) + ";",
                "IsLabelledNewick", true, "adjustTipHeights", false);
        return tree;
    }

    // ── Entry point ───────────────────────────────────────────────────────────

    public static void main(String[] args) throws Exception {
        String[] sizes = IntStream.of(10, 25, 50, 100, 250, 500, 1000, 10000)
                .mapToObj(String::valueOf).toArray(String[]::new);

        Options opt = new OptionsBuilder()
                .include(ScalingBenchmark.class.getSimpleName())
                .param("nLeaves", sizes)
                .resultFormat(ResultFormatType.CSV)
                .result("calibratedcpp-beast/validation/calibratedcpp/heled-and-drummond_validation_and_benchmark/scaling_benchmark.csv")
                .build();

        new Runner(opt).run();
    }
}
