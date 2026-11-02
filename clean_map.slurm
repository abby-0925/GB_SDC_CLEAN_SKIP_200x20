#!/usr/bin/env bash
#SBATCH --job-name=SDCC_MAP
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=/public/home/tanghui/GB_SDC_TrackVis/logs/sdcc_map_%A_%a.out
#SBATCH --error=/public/home/tanghui/GB_SDC_TrackVis/logs/sdcc_map_%A_%a.err
set -euo pipefail
BASE=/public/home/tanghui/GB_SDC_TrackVis
PY=/public/home/tanghui/software/sdc_fast_env/bin/python
export FSLDIR=/public/apps/neuroscience/fsl/6.0.4
export FSLOUTPUTTYPE=NIFTI_GZ
export PATH="$(dirname "$PY"):$FSLDIR/bin:$PATH"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
unset PYTHONPATH PYTHONHOME || true
: "${GROUP:?}"; : "${WAVE:?}"; : "${SLURM_ARRAY_TASK_ID:?}"
RI=$(( WAVE*20 + SLURM_ARRAY_TASK_ID ))
(( RI<601 )) || { echo 'Invalid reference index'; exit 2; }
"$PY" "$BASE/scripts/clean_map.py" --group "$GROUP" --reference-index "$RI" --patient-chunk 50
