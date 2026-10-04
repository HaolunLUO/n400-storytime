# Where the data are

The scripts do not read data from this git repository. They use absolute paths on the ORCD pool, account `haolun52`. Moving the code to a new repository does not move the files.

Pool root:

```
/orcd/pool/005/haolun52
```

## Inputs the scripts open

| What | Path under the pool root |
|---|---|
| Pre-ICA word-locked EEG | `extracted_sections_wordlocked_shared` |
| Shared word features | `extracted_sections_wordlocked_shared/_shared_wordlocked_features` |
| Subject list, 63 children | `dyslexia_natualistics_listing/b2b_decoding/subjects.txt` |
| Group labels | `dyslexia_natualistics_listing/b2b_decoding/joint_v4/cohort_groups.csv` |
| Behaviour workbook | `dyslexia_natualistics_listing/behavioral_data.xlsx` |
| First-pipeline word table and LMM input | `n400_storytime_v1/stimulus/` |
| FastText Chinese vectors | `n400_storytime_v1/stimulus/embeddings/wiki.zh.vec` |
| v1-strict word table, trials, models | `n400_storytime_v1/v1_strict/` |
| Revision LMM and behaviour tables | `n400_storytime_v1/revision/` |
| Unfold features | `stimulus_features_fs100/trf_auditory_v1` |

`n400_common.py` defines the first-pipeline paths. `v1_strict/assemble_trials.py` and `v1_strict/extract_channels.py` repeat the subject list, the group file, and the behaviour workbook. The `b2b_decoding` folder name in those paths is only the directory on disk. These scripts do not import the decoder.

## Where the Slurm jobs run the code

`v1_strict/run_extract.sbatch`, `run_lmm.sbatch`, and `run_plot.sbatch` set:

```
ROOT=/orcd/pool/005/haolun52/_work_Dyslexiab2b_n400rev/b2b_decoding/n400_storytime/v1_strict
```

A clone of this repository in another directory is not what those jobs execute until `ROOT` is changed.

## Aligned eye-ICA epochs

The trigger-aligned eye-ICA refit is not what `extract_channels.py` reads. That script reads `extracted_sections_wordlocked_shared`.

The aligned epochs and the aligned model output are:

```
/orcd/pool/005/haolun52/extracted_sections_wordlocked_ica_aligned
/orcd/pool/005/haolun52/n400_storytime_v1/v1_strict_surprisal_ica_aligned
```

Those paths are not wired into the scripts in this repository.
