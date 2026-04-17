# Vina Docking Integration

P2Rank integrates with [AutoDock Vina](https://vina.scripps.edu/) (and compatible
forks such as Vina-GPU) through the `vina-dock` command.  The pipeline runs
protein preparation, optional P2Rank binding-site prediction, Vina docking, and
unified result reporting in a single command.

## Quick Start

```bash
# P2Rank-guided docking (default): predict pockets, dock into top-3
prank vina-dock -f protein.pdb -vina_ligand ligand.pdbqt

# Normal Vina docking with an explicit search box
prank vina-dock -f protein.pdb -vina_ligand ligand.pdbqt \
      -vina_use_p2rank_pockets 0 \
      -vina_center_x 12.3 -vina_center_y 45.6 -vina_center_z 7.8 \
      -vina_size_x 20 -vina_size_y 20 -vina_size_z 20

# Run on a dataset of proteins
prank vina-dock proteins.ds -vina_ligand ligand.pdbqt
```

## Requirements

| Tool | Purpose | Installation |
|------|---------|--------------|
| **AutoDock Vina ≥ 1.2** | Docking engine | [vina.scripps.edu](https://vina.scripps.edu/) or `conda install -c conda-forge autodock-vina` |

P2Rank checks Vina availability at startup and prints a clear error if the
binary is not found.  You can point P2Rank to a custom Vina binary with the
`-vina_command` parameter.

## Modes

### P2Rank-guided docking (default)

When `vina_use_p2rank_pockets 1` (the default), the pipeline:

1. **Predicts binding pockets** using the P2Rank ML model.
2. Takes the top `vina_max_pockets` predicted pockets (default 3).
3. Runs **one Vina docking job per pocket**, centred on each pocket centroid.
4. Outputs per-pocket results and a consolidated CSV.

This mode is the recommended approach because it automates box placement and
removes the need to manually inspect the structure to pick a search box.

### Normal Vina docking

Set `-vina_use_p2rank_pockets 0` to supply an explicit search box:

```bash
prank vina-dock -f protein.pdb -vina_ligand ligand.pdbqt \
      -vina_use_p2rank_pockets 0 \
      -vina_center_x <X> -vina_center_y <Y> -vina_center_z <Z> \
      -vina_size_x <SX> -vina_size_y <SY> -vina_size_z <SZ>
```

## Pre-docking Preparation

`vina-dock` automatically:

- **Removes water molecules** (HOH/WAT) — controlled by `-vina_prep_remove_water`.
- **Removes HETATM records** (heteroatoms/ligands) — controlled by
  `-vina_prep_remove_heteroatoms`.  Set to `0` to keep cofactors.
- Writes a cleaned `.pdb` and a minimal `.pdbqt` receptor file suitable for Vina.
- Saves a `prep_metadata.txt` file documenting the preparation settings for
  reproducibility.

> **Note:** The built-in PDBQT converter assigns zero partial charges and does not
> compute full atom-type torsion trees.  For high-accuracy docking, pre-prepare
> your receptor with [MGLTools](https://ccsb.scripps.edu/mgltools/) or
> [OpenBabel](https://openbabel.org/) and pass the `.pdbqt` directly.

## Output

For each protein a subdirectory `<outdir>/<protein_label>/` is created:

```
<outdir>/
  <protein_label>/
    prep/
      <protein_label>_cleaned.pdb         # cleaned receptor PDB
      <protein_label>_receptor.pdbqt      # Vina-ready receptor
      prep_metadata.txt                   # preparation settings
    docking_pocket1/
      docked_poses.pdbqt                  # Vina docked poses
      vina.log                            # Vina stdout/stderr log
    docking_pocket2/
      …
    <protein_label>_docking_results.csv   # consolidated results
```

### Results CSV columns

| Column | Description |
|--------|-------------|
| `protein` | Protein label |
| `pocket` | Pocket/box identifier (e.g. `pocket1`, `box1`) |
| `box_center_x/y/z` | Search box centre (Å) |
| `box_size_x/y/z` | Search box dimensions (Å) |
| `num_poses` | Number of docked poses returned by Vina |
| `best_affinity` | Best binding affinity (kcal/mol, most negative = best) |
| `all_affinities` | Semicolon-separated affinities for all poses |
| `output_pdbqt` | Path to the output PDBQT file with all poses |
| `log_file` | Path to the Vina log file |

## Parameters Reference

Override any parameter on the command line with `-param value`.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `vina_command` | `vina` | Vina binary command or full path |
| `vina_ligand` | *(required)* | Path to ligand `.pdbqt` (or `.sdf`/`.mol2`) |
| `vina_use_p2rank_pockets` | `true` | Use P2Rank pockets as docking centres |
| `vina_max_pockets` | `3` | Max P2Rank pockets to dock into |
| `vina_center_x/y/z` | `0.0` | Explicit box centre (used when guided mode off) |
| `vina_size_x/y/z` | `20.0` | Box size in each dimension (Å) |
| `vina_num_modes` | `9` | Max output poses |
| `vina_energy_range` | `3.0` | Energy range (kcal/mol) for pose filtering |
| `vina_exhaustiveness` | `8` | Search exhaustiveness |
| `vina_cpu` | `0` | CPU cores for Vina (0 = auto) |
| `vina_keep_output` | `true` | Keep prep + docking files after run |
| `vina_prep_remove_water` | `true` | Remove water molecules before docking |
| `vina_prep_remove_heteroatoms` | `true` | Remove HETATM records before docking |

## Examples

```bash
# Dock a ligand into the top predicted pocket
prank vina-dock -f test_data/1fbl.pdb -vina_ligand ligand.pdbqt

# Dock into top-5 pockets with higher exhaustiveness
prank vina-dock -f test_data/1fbl.pdb -vina_ligand ligand.pdbqt \
      -vina_max_pockets 5 -vina_exhaustiveness 16

# Keep cofactors (do not strip HETATM)
prank vina-dock -f test_data/1fbl.pdb -vina_ligand ligand.pdbqt \
      -vina_prep_remove_heteroatoms 0

# Use a custom Vina binary
prank vina-dock -f test_data/1fbl.pdb -vina_ligand ligand.pdbqt \
      -vina_command /opt/vina/bin/vina

# Batch docking on a dataset
prank vina-dock proteins.ds -vina_ligand ligand.pdbqt -o docking_results/
```

## Tips

- For large-scale screens, increase `-threads` to process proteins in parallel
  (each thread handles one protein; Vina itself is single-threaded by default).
- Use `-vina_cpu 4` to allocate 4 CPU cores to each Vina invocation.
- If you see "no pockets found", try lowering the P2Rank scoring threshold via
  the standard `-<param>` mechanism or increase the protein resolution.
