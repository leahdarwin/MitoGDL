############################################################
# Script: CIcontact_sourceData.R
# Author: Leah Darwin
# Purpose: Standalone (no PyMOL) reproduction of the mito/nuclear
#          contact-interface computation done in CIcontact_plot.py,
#          for exporting the Fig 5 source-data table. PyMOL is
#          expensive to run, so this parses the mmCIF structure
#          directly and computes the same "within 4.0 Angstrom"
#          contact test as a nearest-neighbor distance check.
############################################################

## ---------------------------------------------------------
## Load required packages (install if missing)
## ---------------------------------------------------------
packages <- c("RANN")
installed <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed)) {
    install.packages(p, dependencies = TRUE)
  }
}
lapply(packages, library, character.only = TRUE)
source("scripts/figures/sourceData_helpers.R")

## ---------------------------------------------------------
## Parse ATOM/HETATM records directly out of the mmCIF file
## ---------------------------------------------------------
# _atom_site loop column order in data/8B9Z.cif:
# group_PDB id type_symbol label_atom_id label_alt_id label_comp_id
# label_asym_id label_entity_id label_seq_id ins_code x y z occupancy
# b_iso charge auth_seq_id auth_comp_id chain auth_atom_id model_num
cif_cols <- c("group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id",
              "label_comp_id", "label_asym_id", "label_entity_id", "label_seq_id",
              "ins_code", "x", "y", "z", "occupancy", "b_iso", "charge",
              "auth_seq_id", "auth_comp_id", "chain", "auth_atom_id", "model_num")

cif_lines <- readLines("data/8B9Z.cif")
atom_lines <- cif_lines[startsWith(cif_lines, "ATOM") | startsWith(cif_lines, "HETATM")]

atoms <- read.table(text = atom_lines, stringsAsFactors = FALSE, col.names = cif_cols)
atoms$x <- as.numeric(atoms$x)
atoms$y <- as.numeric(atoms$y)
atoms$z <- as.numeric(atoms$z)
atoms$auth_seq_id <- as.integer(atoms$auth_seq_id)
# PyMOL's "chain" selector on a loaded mmCIF is the auth_asym_id,
# which read.table already placed in the "chain" column above.

## ---------------------------------------------------------
## Read TSV of AA changes (chain_id, AA_POS, MITO)
## ---------------------------------------------------------
df <- read.delim("data/CI_chain_snps.tsv", stringsAsFactors = FALSE)

mt_chains <- unique(df$chain_id)
mt_atoms  <- atoms[atoms$chain %in% mt_chains, ]
nuc_atoms <- atoms[!(atoms$chain %in% mt_chains), ]

## ---------------------------------------------------------
## Compute contact interface: mito-chain atoms within 4.0 A of any
## nuclear-chain atom (equivalent to PyMOL's
## "mt_chain within 4.0 of nuc_chain" selection)
## ---------------------------------------------------------
nn <- RANN::nn2(
  data  = as.matrix(nuc_atoms[, c("x", "y", "z")]),
  query = as.matrix(mt_atoms[, c("x", "y", "z")]),
  k = 1
)
mt_atoms$is_interface <- nn$nn.dists[, 1] <= 4.0

interface_residues <- unique(paste(mt_atoms$chain[mt_atoms$is_interface],
                                    mt_atoms$auth_seq_id[mt_atoms$is_interface]))

## ---------------------------------------------------------
## For each mito-lineage option (BZ/siI/yak), reproduce the
## highlight + overlap logic from CIcontact_plot.py and export
## the per-residue source-data table
## ---------------------------------------------------------
options <- c("BZ", "siI", "yak")

for (option in options) {
  match <- if (option == "BZ") {
    grepl("B", df$MITO) | grepl("Z", df$MITO)
  } else {
    grepl(option, df$MITO, fixed = TRUE)
  }

  residue_key <- paste(df$chain_id, df$AA_POS)
  overlap <- rep(NA, nrow(df))
  overlap[match] <- residue_key[match] %in% interface_residues

  out <- data.frame(
    chain_id    = df$chain_id,
    AA_POS      = df$AA_POS,
    MITO        = df$MITO,
    highlighted = match,
    overlap     = overlap
  )

  save_source_data(out, paste0("CIcontact_fig5_", option), "main_figs")
}
