# Imitation — voxel-level lesion-symptom maps

Voxel-based lesion-symptom mapping (VLSM) of two imitation scores (hand and
finger), and of their difference, plus the "Achilles" re-analysis that asks
how much the maps change with the covariates you control for.

## Layout

```
voxel_maps/
  README.md
  config_voxel_maps.m      every setting (columns, cut-offs, FDR, permutations, covariate sets)
  setup_voxel_maps.m       adds the folders below to the MATLAB path
  example_voxel_maps.m     a worked session, section by section
  analysis/
    vlsm_analysis.m              hand, finger and difference maps, with FDR and FWER thresholds
    compare_covariate_control.m  the Achilles comparison: same patients, different covariate sets
    interaction_analysis.m       voxels where the lesion effect depends on a factor
    save_vlsm_maps.m             writes the maps as .nii images
    vlsm_ttest.m  vlsm_ranksum.m  vlsm_bayes.m  vlsm_lpd.m  vlsm_interaction.m   the statistics
  shared/                  data loading, normalisation, Bayes factors, FDR, permutations, image I/O
  tests/VoxelMapsTest.m    runtests('VoxelMapsTest')
```

## Setup

MATLAB R2019b or later with the Statistics and Machine Learning Toolbox. SPM12
on the path to read and write images. Parallel Computing Toolbox is optional
(`cfg.workers` for permutations).

```matlab
run('...\voxel_maps\setup_voxel_maps.m')
runtests('VoxelMapsTest')          % synthetic data; no SPM needed except the image tests
edit example_voxel_maps
```

## Data

`load_imitation_data(cfg)` reads:

- a spreadsheet with one row per patient (header row), with an ID column, the
  two score columns and any covariates. A covariate can be a sum of columns:
  `'TimePost=TimeYears+TimeMonths'`;
- a folder of binary lesion images. Each file is matched to a patient by the
  ID taken from its name (`cfg.data.lesion_id_regex`), exactly and ignoring
  case. Patients without an image are listed in `D.missing_lesions`, images
  without a patient in `D.unmatched_files`, so nothing is dropped silently;
- optionally a brain mask, and the GetLesions preprocessing (fill gaps between
  slices, flip every lesion to one hemisphere) via `cfg.data.preprocess`.

If you still have the old `Data` cell array from `ConnectLesionsToData`,
`legacy_cell_to_data` converts it (see the example).

## What each map is

| method    | statistic per voxel | sign |
|-----------|---------------------|------|
| `ttest`   | Welch t (or pooled, `ttest_vartype`), unlesioned − lesioned | t > 0: lesion lowers the score |
| `ranksum` | Wilcoxon rank-sum z, tie- and continuity-corrected | same |
| `bayes`   | log Bayes factor for a group difference (JZS prior by default) | > 0 favours an effect |
| `lpd`     | lesion–impairment association (Rudrauf et al. 2008), impaired = raw score below `impairment_cutoffs` | |

Scores are normalised first (`cfg.scores.method`: rank-based inverse normal
after regressing out `cfg.scores.covariates`, as `Norm_Resid_MinusThresh`).
The difference map uses hand − finger of the normalised scores (Achilles.m)
or normalises the raw difference (`cfg.vlsm.difference`, Achilles2017.m).

Voxels are tested if lesioned in at least `min_lesion_fraction` of patients
and at least `min_lesion_count`, and spared in at least `min_lesion_count`.
FDR (BH or BY) counts only tested voxels. With `n_perm > 0`, family-wise
thresholds come from permuting whole patients, so hand, finger and the
difference stay paired.

