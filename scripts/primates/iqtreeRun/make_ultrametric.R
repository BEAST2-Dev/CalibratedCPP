#!/usr/bin/env Rscript
# Convert the IQ-TREE ML tree (genetic-distance branch lengths) into an ultrametric
# chronogram (time-based branch lengths, tips aligned at present) using ape::chronos(),
# supplying the same 13 fossil calibration bounds used throughout this project's BEAST
# analyses (from primates.lphy) as calibration constraints, rather than just rescaling
# relative to one arbitrary point.
#
# Usage: Rscript make_ultrametric.R [input.treefile] [output.nwk]

suppressMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
input_path  <- if (length(args) >= 1) args[1] else "primates_iqtree.treefile"
output_path <- if (length(args) >= 2) args[2] else "primates_iqtree_startingTree.nwk"

tree <- read.tree(input_path)
cat("Loaded tree:", length(tree$tip.label), "tips\n")

# IQ-TREE's -o flag only repositions the named outgroup for display; the underlying tree
# is still an unrooted bifurcation with a basal trifurcation (BigClade, Mus_musculus,
# Tupaia_chinensis) rather than the properly nested (Mus_musculus,(Tupaia_chinensis,...))
# structure our calibration scheme assumes. Without this, getMRCA() for "all tips except
# Mus_musculus" resolves to the same node as the full-tree root, corrupting the
# Euarchontoglires/Euarchonta calibrations. Explicitly re-root and force full bifurcation.
tree <- root(tree, outgroup = "Mus_musculus", resolve.root = TRUE)
cat("Re-rooted on Mus_musculus with forced bifurcation\n")

# Same 13 calibration clades (taxa + bounds, in Ma) as calibratedcpp-beast/examples/primates/primates.lphy
calibrations <- list(
  list(taxa = tree$tip.label, lower = 65.79, upper = 125.816),  # Euarchontoglires (root)
  list(taxa = setdiff(tree$tip.label, "Mus_musculus"), lower = 65.79, upper = 125.816),  # Euarchonta
  list(taxa = c("Otolemur_garnettii", "Microcebus_murinus", "Propithecus_coquereli",
                "Carlito_syrichta", "Cebus_capucinus_imitator", "Saimiri_boliviensis",
                "Aotus_nancymaae", "Callithrix_jacchus", "Nomascus_leucogenys", "Pongo_abelii",
                "Gorilla_gorilla", "Homo_sapiens", "Pan_paniscus", "Pan_troglodytes",
                "Colobus_angolensis_palliatus", "Piliocolobus_tephrosceles", "Rhinopithecus_bieti",
                "Rhinopithecus_roxellana", "Chlorocebus_sabaeus", "Macaca_nemestrina",
                "Macaca_fascicularis", "Macaca_mulatta", "Cercocebus_atys", "Mandrillus_leucophaeus",
                "Papio_anubis", "Theropithecus_gelada"), lower = 55.935, upper = 66.095),  # Primates
  list(taxa = c("Otolemur_garnettii", "Microcebus_murinus", "Propithecus_coquereli"),
       lower = 18.5, upper = 55.8),  # Lorisiformes
  list(taxa = c("Colobus_angolensis_palliatus", "Piliocolobus_tephrosceles", "Rhinopithecus_bieti",
                "Rhinopithecus_roxellana", "Chlorocebus_sabaeus", "Macaca_nemestrina",
                "Macaca_fascicularis", "Macaca_mulatta", "Cercocebus_atys", "Mandrillus_leucophaeus",
                "Papio_anubis", "Theropithecus_gelada"), lower = 12.47, upper = 25.235),  # Cercopithecidae
  list(taxa = c("Colobus_angolensis_palliatus", "Piliocolobus_tephrosceles", "Rhinopithecus_bieti",
                "Rhinopithecus_roxellana"), lower = 8.125, upper = 15.0),  # Colobinae
  list(taxa = c("Chlorocebus_sabaeus", "Macaca_nemestrina", "Macaca_fascicularis", "Macaca_mulatta",
                "Cercocebus_atys", "Mandrillus_leucophaeus", "Papio_anubis", "Theropithecus_gelada"),
       lower = 6.5, upper = 15.0),  # Cercopithecinae
  list(taxa = c("Macaca_nemestrina", "Macaca_fascicularis", "Macaca_mulatta", "Cercocebus_atys",
                "Mandrillus_leucophaeus", "Papio_anubis", "Theropithecus_gelada"),
       lower = 5.33, upper = 12.51),  # Papionini
  list(taxa = c("Nomascus_leucogenys", "Pongo_abelii", "Gorilla_gorilla", "Homo_sapiens",
                "Pan_paniscus", "Pan_troglodytes"), lower = 13.4, upper = 25.235),  # Hominoidea
  list(taxa = c("Pongo_abelii", "Gorilla_gorilla", "Homo_sapiens", "Pan_paniscus", "Pan_troglodytes"),
       lower = 12.3, upper = 25.235),  # Hominidae
  list(taxa = c("Homo_sapiens", "Pan_paniscus", "Pan_troglodytes"), lower = 4.631, upper = 15.0),  # Homo+Pan
  list(taxa = c("Cebus_capucinus_imitator", "Saimiri_boliviensis", "Aotus_nancymaae", "Callithrix_jacchus"),
       lower = 13.183, upper = 34.5),  # Callitrichidae+Cebidae
  list(taxa = c("Cebus_capucinus_imitator", "Saimiri_boliviensis"), lower = 13.032, upper = 34.5)  # Cebidae
)

calib_rows <- do.call(rbind, lapply(calibrations, function(cal) {
  node <- if (length(cal$taxa) == length(tree$tip.label)) {
    Ntip(tree) + 1  # root node
  } else {
    getMRCA(tree, cal$taxa)
  }
  data.frame(node = node, age.min = cal$lower, age.max = cal$upper, soft.bounds = FALSE)
}))

cat("Calibration nodes resolved:\n")
print(calib_rows)

chrono_tree <- chronos(tree, calibration = calib_rows, model = "relaxed")

# Post-hoc numerical safety margin -- NOT a change to the fit. chronos() fits ages against
# the exact bounds above, and correctly lands exactly on a boundary when that bound is the
# binding constraint (the mathematically correct constrained optimum). But writing an
# exactly-on-the-boundary value to text and having a *different* program (BEAST/Java)
# re-sum branch lengths in its own floating-point order can push the reconstructed age a
# few ULPs outside the hard Uniform bound, which BEAST then rejects outright (zero prior
# density) when trying to initialise. Nudge only calibrated nodes that land within EPS of
# their own bound, by EPS, toward the interior -- purely so the value survives the R -> text
# -> Java round trip; the calibration-constrained fit itself (chronos() call above) is
# completely untouched.
EPS <- 1e-4  # 100 years; negligible at Ma timescales
depths <- node.depth.edgelength(chrono_tree)  # root-to-node distance, ape node numbering
root_depth <- max(depths[1:Ntip(chrono_tree)])  # tip depth = root-to-tip distance (ultrametric)

n_nudged <- 0
for (i in seq_len(nrow(calib_rows))) {
  node <- calib_rows$node[i]
  lo <- calib_rows$age.min[i]
  hi <- calib_rows$age.max[i]
  age <- root_depth - depths[node]
  # age = root_depth - depths[node], so increasing age requires DECREASING depth, and
  # decreasing age requires INCREASING depth -- opposite sign from age itself.
  if (age - lo < EPS) {
    depths[node] <- depths[node] - EPS  # decrease depth -> increase age, away from lower bound
    n_nudged <- n_nudged + 1
  } else if (hi - age < EPS) {
    depths[node] <- depths[node] + EPS  # increase depth -> decrease age, away from upper bound
    n_nudged <- n_nudged + 1
  }
}
cat("Nudged", n_nudged, "boundary-touching node(s) by", EPS, "Ma for numerical safety\n")

# Rebuild branch lengths from the (locally nudged) depths so the tree stays fully consistent
chrono_tree$edge.length <- depths[chrono_tree$edge[, 2]] - depths[chrono_tree$edge[, 1]]
if (any(chrono_tree$edge.length < 0)) {
  stop("Nudging produced a negative branch length -- two adjacent nudges collided; reduce EPS.")
}

write.tree(chrono_tree, output_path)
cat("Wrote ultrametric tree to", output_path, "\n")
cat("Root age:", max(node.depth.edgelength(chrono_tree)), "\n")
