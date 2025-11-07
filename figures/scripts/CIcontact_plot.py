############################################################
# Script: CIcontact_plot.R
# Author: Leah Darwin
# Date: 2025-11-7
# Purpose: Use pymol to highlight contact sites for mtDNA and
#           nuclear proteins given a csv of AA changes.
# Dependencies: pymol, pandas
############################################################


from pymol import cmd
import pandas as pd

cmd.load("../extra_data/8B9Z.cif")

# ===== 1. Read TSV file =====
# Example TSV structure:
# chain_id   AA_POS   MITO
# A          45       yak_mt
# B          109      mito_hs
tsv_file = "../extra_data/CI_chain_snps.tsv"
df = pd.read_csv(tsv_file, sep="\t")

##Select D.mel="BZ", D.sim="siI", or D.yak="yak" to highlight AA changes
#option="BZ"
#option="yak"
option="siI"

# ===== 2. Get chain lists =====
mt_chains = df["chain_id"].unique().tolist()

# ===== 3. Make selections =====
cmd.select("mt_chain", "chain " + "+".join(mt_chains))
cmd.select("nuc_chain", "not mt_chain")

# Color them
cmd.color("green", "mt_chain")
cmd.color("tv_yellow", "nuc_chain")

##only show ribbon for Nuc proteins 
cmd.show("cartoon", "nuc_chain")
cmd.hide("spheres", "nuc_chain")
cmd.hide("sticks", "nuc_chain")
cmd.hide("sticks", "mt_chain")

##set transparency 
cmd.set("cartoon_transparency", 0.4, "nuc_chain")
cmd.set("cartoon_transparency", 0.4, "mt_chain")

# ===== 4. Compute contact interface =====
cmd.select("mt_if_atoms", "mt_chain within 4.0 of nuc_chain")
#cmd.select("nuc_if_atoms", "nuc_chain within 4.0 of mt_chain")

# Optional: show as spheres
cmd.show("spheres", "mt_if_atoms")
cmd.set("sphere_scale", 0.8)
cmd.color("brightorange", "mt_if_atoms")
cmd.set("sphere_transparency", 0.2, "mt_if_atoms")

# ===== 5. Iterate through TSV and highlight residues =====
for _, row in df.iterrows():
    chain = row["chain_id"]
    pos = int(row["AA_POS"])
    mito_tag = str(row["MITO"])

    match = any(c in mito_tag for c in option) if option == "BZ" else option in mito_tag
    if match:
        # make a selection for the residue
        res_sel = f"chain {chain} and resi {pos}"

        # Check if this residue overlaps with contact sites
        overlap = cmd.count_atoms(f"{res_sel} and mt_if_atoms") > 0

        # Decide color logic
        if overlap:
            color = "red"
            print(f"Residue {res_sel} in chain {chain} overlaps with mt_if_atoms")
        else:
            color = "green"

        # Apply color
        cmd.show("spheres", res_sel)
        cmd.color(color, res_sel)
        cmd.set("sphere_scale", 1.5, res_sel)
        
        with open("residue_log_" + option + ".txt", "a") as f:
                    f.write(f"{pos},{chain},{overlap}\n")

cmd.set_view((
     0.188485950,    0.287537843,   -0.939039052,
     0.443885118,   -0.877878726,   -0.179711923,
    -0.876035929,   -0.382951707,   -0.293101013,
     0.000000000,    0.000000000, -826.158935547,
   231.672821045,  217.831970215,  244.280120850,
   651.349975586, 1000.967895508,  -20.000000000
))

# White background
cmd.bg_color("white")
cmd.set("ray_opaque_background", 0)

# --- Lighting / contrast settings ---
cmd.set("ambient", 0.25)
cmd.set("direct", 0.9)
cmd.set("reflect", 0.25)
cmd.set("shininess", 40)
cmd.set("specular", 0.5)
cmd.set("light_count", 10)
cmd.set("two_sided_lighting", "on")
cmd.set("spec_power", 80)
cmd.set("spec_reflect", 0.4)

# --- Rendering clarity tweaks ---
cmd.set("transparency_mode", 2)
cmd.set("depth_cue", 0)
cmd.set("antialias", 2)
cmd.set("ray_trace_mode", 1)
cmd.set("ray_shadows", 1)
cmd.set("ambient_occlusion_mode", 2)


# Ray trace and save PNG
cmd.png("../main_figs/CIcontact_fig5_" + option + ".png", dpi=500, width=4000, height=4000, ray=1, quiet=0)

