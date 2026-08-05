# PhyloData Calibration Forest Constraint Trees

Calibration constraint trees extracted from BEAST2 XML files in the PhyloData repository.
These represent real-world calibration structures used in published phylogenetic analyses.

## Data Pipeline

Starting from 222 BEAST2 XML files downloaded via [phylodata](https://pypi.org/project/phylodata/):

| Stage | Count | Reason |
|-------|-------|--------|
| PhyloData BEAST2 XMLs | 222 | Raw input |
| Files with MRCAPrior calibrations | 157 | 65 have no node-age calibrations |
| After removing single-taxon calibrations | 136 | Tip-dating leaf constraints only |
| After removing tip-dated analyses | 128 | Non-contemporaneous taxa not supported |
| Unique calibration forests | 96 | Deduplicated by calibration structure |
| After removing overlapping calibrations | 93 | Non-nested overlaps are methodologically invalid |

**Filtering details:**
- **Tip-dated analyses**: CalibratedCPP assumes contemporaneous tips
- **Single-taxon calibrations**: These are tip-dating constraints, not clade calibrations
- **Overlapping calibrations**: 3 files (hara-2021, harrington-wn6u, serrano-zc3e) have calibrations with non-nested taxon set overlaps - likely researcher errors in the original XML
- **Multi-gene duplicates**: Automatically deduplicated (same calibration applied to multiple gene trees → keep one)

## Format

Annotated Newick format with BEAST-style metadata:
```
((A,B)[&name=Clade1,dist=LogNormal,M=5.2,S=0.5,offset=33.9],(C,D)[&name=Clade2])[&name=Root];
```

## Benchmark Results

All 93 constraint trees benchmark successfully with CalibratedCoalescentPointProcess.

![](benchmark_analysis_plot.svg)

**Regression model** (R² = 0.99):
```
time(μs) = 0.010 × complexity + 0.070 × taxa + 94.3
```

Where complexity = Σ 2^k × k³ for each calibration node with k children.

**Top 5 by execution time:**
<!-- start of table -->
| File | Taxa | Calibrations | Complexity | Time (μs) |
|------|------|--------------|------------|-----------|
| burridge-2020-migration-nbmn_Galaxioid_Diversification | 725 | 17 | 2,726,108 | 27,958 |
| harrington-2024-dispersal-we1f_BEAST_divtimes_set3 | 5,940 | 22 | 1,024,106 | 11,050 |
| trotta-2018-community-a38u_Trotta_PineRockland | 1,080 | 50 | 19,816 | 2,680 |
| v-2016-new-g4qu_Scleractinia | 1,158 | 8 | 1,058 | 992 |
| marburger-2018-whole-8zqj_cory_final_bmodel | 1,105 | 7 | 43,904 | 536 |
<!-- end of table -->
## Usage

Run benchmark:
```bash
java -Xss32m -cp "lib/*:target/classes:target/test-classes" \
    calibratedcpp.SimpleConstraintBenchmark validation/calibratedcpp/phylodata_calibration_forests/
```

Generate plot:
```bash
python validation/calibratedcpp/phylodata_calibration_forests/plot_benchmark.py
```

## Output Files

- `benchmark_results.csv` - Full benchmark results with timing statistics
- `benchmark_analysis_plot.png` - Visualization of time vs complexity
- `*.newick` - 93 constraint tree files

## Source

Generated from [calibrated-time-tree-dist](https://github.com/adru001/calibrated-time-tree-dist) analysis.
