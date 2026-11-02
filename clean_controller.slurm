#!/usr/bin/env bash
#SBATCH --job-name=SDCC_CTL
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=512M
#SBATCH --time=00:05:00
#SBATCH --output=/public/home/tanghui/GB_SDC_TrackVis/logs/sdcc_ctl_%j.out
#SBATCH --error=/public/home/tanghui/GB_SDC_TrackVis/logs/sdcc_ctl_%j.err
set -euo pipefail
BASE=/public/home/tanghui/GB_SDC_TrackVis
ROOT="$BASE/ref_major_clean_v1"
GROUP=${1:?Group required}; WAVE=${2:?Wave required}
[[ "$GROUP" =~ ^g[[:alnum:]]+$ && "$WAVE" =~ ^(0|[1-9][0-9]*)$ ]] || exit 2
(( WAVE<=30 )) || exit 2
MANIFEST="$ROOT/manifests/$GROUP.patients.txt"
[[ -s "$MANIFEST" ]] || exit 2
PATIENTS=$(awk 'NF{n++}END{print n+0}' "$MANIFEST")
(( PATIENTS>=1 && PATIENTS<=200 )) || { echo "Unexpected patient count: $PATIENTS";exit 2; }
[[ $(awk 'NF{n++}END{print n+0}' "$ROOT/manifests/references.txt") == 601 ]] || exit 2
mkdir -p "$ROOT/submissions/$GROUP" "$ROOT/states/$GROUP"
exec 9>"$ROOT/submissions/$GROUP/controller.lock"
flock -n 9 || { echo 'Other controller active'; exit 3; }
HASH=$(sha256sum "$MANIFEST" | cut -d' ' -f1)
HF="$ROOT/submissions/$GROUP/manifest.sha256"
if [[ -e "$HF" ]]; then [[ $(cat "$HF") == "$HASH" ]] || { echo 'Manifest modified';exit 4; }; else echo "$HASH" > "$HF";fi
if (( WAVE>0 )); then
 PREV=$((WAVE-1))
 for ((c=0;c<(PATIENTS+49)/50;c++));do
  printf -v FN 'c%03d_wave%03d.ok' "$c" "$PREV"
  [[ -s "$ROOT/states/$GROUP/$FN" ]] || { echo "Missing prior $FN";exit 5; }
 done
fi
printf -v TAG 'w%03d' "$WAVE"
REC="$ROOT/submissions/$GROUP/$TAG.jobs"
[[ ! -e "$REC" ]] || { echo "Already submitted $REC"; exit 6; }
REMAIN=$((601-20*WAVE)); REFS=$(( REMAIN<20 ? REMAIN : 20 ))
MAP=$(sbatch --parsable --array="0-$((REFS-1))%20" --export="ALL,GROUP=$GROUP,WAVE=$WAVE" "$BASE/scripts/clean_map.slurm")
printf 'MAP_JOB=%s\n' "$MAP" > "$REC"
RED=$(sbatch --parsable --dependency="afterok:$MAP" --export="ALL,GROUP=$GROUP,WAVE=$WAVE" "$BASE/scripts/clean_reduce.slurm")
printf 'REDUCE_JOB=%s\n' "$RED" >> "$REC"
echo "[SUBMITTED] wave=$WAVE refs=$REFS map=$MAP reduce=$RED"
if (( WAVE<30 ));then
 NEXT=$((WAVE+1))
 CTL=$(sbatch --parsable --dependency="afterok:$RED" "$BASE/scripts/clean_controller.slurm" "$GROUP" "$NEXT")
 printf 'NEXT_CONTROLLER=%s\n' "$CTL" >> "$REC"
 echo "[AUTO_NEXT] controller=$CTL wave=$NEXT"
else
 echo '[FINAL_WAVE_SUBMITTED] Check 601/601 and final outputs after Reduce succeeds.'
fi
