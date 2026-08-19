# Likelihood Benchmark and Validation

We benchmark the computation time for the likelihood calculation of a single fully balanced tree with 128 leaves as we increase the number of cherries for which we want to condition the age. The time is compared with the [Heled and Drummond (2015)](https://doi.org/10.1093/sysbio/syu089) implementation of the 'Calibrated Birth-Death process'.

We also validate the likelihood computation by comparing log-likelihood values for a fixed tree varying the birth-rate with both the Heled and Drummond and Calibrated CPP likelihood implementations.

## Results

**Likelihood Benchmark: CPP v Heled and Drummond**
<!-- start table -->

|   Number of Calibrations |   Time (μs) Calibrated CPP |   Time (μs) Heled & Drummond |
|--------------------------|----------------------------|------------------------------|
|                      1.0 |                       2.50 |                         6.38 |
|                      2.0 |                       3.46 |                       229.27 |
|                      3.0 |                       5.99 |                      9611.06 |
|                      4.0 |                      14.59 |                    338659.46 |
|                      5.0 |                      45.92 |                  11306571.95 |
|                      6.0 |                     154.16 |                 270543394.43 |

<!-- end table -->

**Benchmark and Validation Plot**
![](./combined_plot.png)

## Usage

Run benchmark (takes approximately 50 minutes):
```bash
MAVEN_OPTS="-Xss32m" mvn -pl calibratedcpp-beast test-compile exec:java \
  -Dexec.mainClass=calibratedcpp.LikelihoodBenchmark \
  -Dexec.classpathScope=test
```

Run validation:
```bash
mvn -pl calibratedcpp-beast test \
  -Dtest=CalibratedCoalescentPointProcessTest#heledAndDrummondComparison
```

Generate plot:
```bash
python calibratedcpp-beast/validation/calibratedcpp/heled-and-drummond_validation_and_benchmark/plot_output.py
```

## Output Files

 - [`benchmark_results.csv`](benchmark_results.csv).
 - [`validation_results.csv`](validation_results.csv).
 - [`combined_benchmark_validation.png`](./combined_benchmark_validation.png).
