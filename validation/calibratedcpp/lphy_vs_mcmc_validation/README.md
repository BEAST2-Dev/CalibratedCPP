# LPhy vs MCMC validation

Validates `CalibratedCPPTree` by checking that its forward-simulation prior
matches the BEAST sample-from-prior posterior (via Welch's t-test) for four
scenarios.

## Scenarios

| ID                     | n    | λ   | μ | ρ   | Condition    | Calibrations                                                         |
|------------------------|------|-----|---|-----|--------------|----------------------------------------------------------------------|
| `fixStem`              | 100  | 2.0 | 1 | 0.1 | stemAge = 3  | (leaf_1,leaf_2)=1, leaf_45–51=1.2, leaf_45–55=1.5, leaf_87–90=2    |
| `fixRoot`              | 100  | 2.0 | 1 | 0.1 | rootAge = 3  | all 100 taxa = 3, (leaf_1,leaf_2)=1, leaf_45–51=1.2, leaf_45–55=1.5, leaf_87–90=2 |
| `fix4LeafStem`         | 4    | 2.0 | 1 | 0.1 | stemAge = 3  | leaf_1–2=1, leaf_3–4=1.5                                            |
| `fix4LeafRoot`         | 4    | 2.0 | 1 | 0.1 | rootAge = 3  | all 4 taxa = 3, leaf_1–2=1, leaf_3–4=1.5                           |

Each scenario is run both conditioned and unconditioned (`conditionOnCalibrations="false"`)
to confirm that conditioning changes the distribution.

---

## Prerequisites

| Requirement | Notes |
|---|---|
| Java 25 | Must be the JDK used by Maven (see developer guide) |
| CalibratedCPP built | `mvn install -DskipTests` from repo root |
| BEAST + `applauncher` | On `$PATH` |
| Python ≥ 3.9 | `numpy`, `scipy`, `matplotlib` |

---

## Step 1 — Generate BEAST XMLs and LPhy trees

From the project root, run `lphybeast convert` for each scenario to produce the BEAST XML, then `SLPhy` to forward-simulate 10 000 trees. The SLPhy output is named `<base>_tree.trees`; rename it to match the `LPhy` convention.

The `_not_cond` XMLs are hand-authored variants with `conditionOnCalibrations="false"` and are already present in this directory.

```bash
VAL=$(pwd)/calibratedcpp-beast/validation/calibratedcpp/lphy_vs_mcmc_validation

# fixStem
mvn -pl calibratedcpp-lphybeast-launcher exec:exec -P-all \
  -Dlphybeast.args="convert -l 20000000 -le 2000 -o $VAL/fixStemMCMC.xml $VAL/fixStemMCMC.lphy"
mvn -pl calibratedcpp-lphy-studio exec:exec \
  -Dlphy.args="-r 10000 $VAL/fixStemMCMC.lphy"
mv $VAL/fixStemMCMC_tree.trees $VAL/fixStemLPhy.trees
find $VAL -name "fixStemMCMC_r*" -delete

# fixRoot
mvn -pl calibratedcpp-lphybeast-launcher exec:exec -P-all \
  -Dlphybeast.args="convert -l 20000000 -le 2000 -o $VAL/fixRootMCMC.xml $VAL/fixRootMCMC.lphy"
mvn -pl calibratedcpp-lphy-studio exec:exec \
  -Dlphy.args="-r 10000 $VAL/fixRootMCMC.lphy"
mv $VAL/fixRootMCMC_tree.trees $VAL/fixRootLPhy.trees
find $VAL -name "fixRootMCMC_r*" -delete

# fix4LeafStem
mvn -pl calibratedcpp-lphybeast-launcher exec:exec -P-all \
  -Dlphybeast.args="convert -l 20000000 -le 2000 -o $VAL/fix4LeafStemMCMC.xml $VAL/fix4LeafStemMCMC.lphy"
mvn -pl calibratedcpp-lphy-studio exec:exec \
  -Dlphy.args="-r 10000 $VAL/fix4LeafStemMCMC.lphy"
mv $VAL/fix4LeafStemMCMC_tree.trees $VAL/fix4LeafStemLPhy.trees
find $VAL -name "fix4LeafStemMCMC_r*" -delete

# fix4LeafRoot
mvn -pl calibratedcpp-lphybeast-launcher exec:exec -P-all \
  -Dlphybeast.args="convert -l 20000000 -le 2000 -o $VAL/fix4LeafRootMCMC.xml $VAL/fix4LeafRootMCMC.lphy"
mvn -pl calibratedcpp-lphy-studio exec:exec \
  -Dlphy.args="-r 10000 $VAL/fix4LeafRootMCMC.lphy"
mv $VAL/fix4LeafRootMCMC_tree.trees $VAL/fix4LeafRootLPhy.trees
find $VAL -name "fix4LeafRootMCMC_r*" -delete
```

---

## Step 2 — Run BEAST sample-from-prior

Run all eight XMLs (four conditioned + four unconditioned).
Each run samples the prior (no alignment); chain = 20 M, logged every 2 000 steps → 10 000 trees.

```bash
VAL=$(pwd)/calibratedcpp-beast/validation/calibratedcpp/lphy_vs_mcmc_validation

mvn -pl calibratedcpp-beast exec:exec -Dbeast.args="-overwrite -prefix $VAL/ $VAL/fixStemMCMC.xml"
mvn -pl calibratedcpp-beast exec:exec -Dbeast.args="-overwrite -prefix $VAL/ $VAL/fixStemMCMC_not_cond.xml"
mvn -pl calibratedcpp-beast exec:exec -Dbeast.args="-overwrite -prefix $VAL/ $VAL/fixRootMCMC.xml"
mvn -pl calibratedcpp-beast exec:exec -Dbeast.args="-overwrite -prefix $VAL/ $VAL/fixRootMCMC_not_cond.xml"
mvn -pl calibratedcpp-beast exec:exec -Dbeast.args="-overwrite -prefix $VAL/ $VAL/fix4LeafStemMCMC.xml"
mvn -pl calibratedcpp-beast exec:exec -Dbeast.args="-overwrite -prefix $VAL/ $VAL/fix4LeafStemMCMC_not_cond.xml"
mvn -pl calibratedcpp-beast exec:exec -Dbeast.args="-overwrite -prefix $VAL/ $VAL/fix4LeafRootMCMC.xml"
mvn -pl calibratedcpp-beast exec:exec -Dbeast.args="-overwrite -prefix $VAL/ $VAL/fix4LeafRootMCMC_not_cond.xml"
```

Produces one `.trees` file per XML (8 files total).

---

## Step 3 — Compute tree statistics

Uses BEAST's `TreeStat2` (via `applauncher`) to compute tree length, height,
gamma, Colless, B1, and cherry count for every .trees file.

```bash
bash calibratedcpp-beast/validation/calibratedcpp/lphy_vs_mcmc_validation/compute_tree_stats.sh
```

Produces one `*.treestreestats.log` file per trees file (12 files total).

## Step 4 — Visualise distributions

```bash
cd calibratedcpp-beast/validation/calibratedcpp/lphy_vs_mcmc_validation/
python compare_distributions.py
```

Reads the `*.treestreestats.log` files and writes three PDFs:
- `fixStem_distributions.pdf`
- `fixRoot_distributions.pdf`
- `fix4Leaf_distributions.pdf`  (combined 2×3 figure, both 4-leaf scenarios)

Also prints Welch's t-test p-values to stdout for each statistic and scenario.

---

## Quick-reference: expected output files

```
fixStemMCMC.xml                              ← generated by Step 1 (lphybeast convert)
fixRootMCMC.xml
fix4LeafStemMCMC.xml
fix4LeafRootMCMC.xml
fixStemMCMC_not_cond.xml                     ← hand-authored (already present)
fixRootMCMC_not_cond.xml
fix4LeafStemMCMC_not_cond.xml
fix4LeafRootMCMC_not_cond.xml

fixStemLPhy.trees                            ← generated by Step 1 (SLPhy)
fixRootLPhy.trees
fix4LeafStemLPhy.trees
fix4LeafRootLPhy.trees

fixStemMCMC.trees                            ← output of BEAST runs (Step 2)
fixStemMCMC_not_cond.trees
fixRootMCMC.trees
fixRootMCMC_not_cond.trees
fix4LeafStemMCMC.trees
fix4LeafStemMCMC_not_cond.trees
fix4LeafRootMCMC.trees
fix4LeafRootMCMC_not_cond.trees

fixStemLPhy.treestreestats.log               ← output of compute_tree_stats.sh (Step 3)
fixRootLPhy.treestreestats.log
fix4LeafStemLPhy.treestreestats.log
fix4LeafRootLPhy.treestreestats.log
fixStemMCMC.treestreestats.log
fixStemMCMC_not_cond.treestreestats.log
fixRootMCMC.treestreestats.log
fixRootMCMC_not_cond.treestreestats.log
fix4LeafStemMCMC.treestreestats.log
fix4LeafStemMCMC_not_cond.treestreestats.log
fix4LeafRootMCMC.treestreestats.log
fix4LeafRootMCMC_not_cond.treestreestats.log

fixStem_distributions.pdf                    ← output of compare_distributions.py (Step 4)
fixRoot_distributions.pdf
fix4Leaf_distributions.pdf
```
