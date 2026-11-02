GB SDC CLEAN / 200 patients * 20 references / quarantine on any failure

Uses prior validated fast_batch.py and shared_invwarp untouched.
Make new manifest in ref_major_clean_v1/manifests. Never copy old map/state.
A patient with EMPTY_SOURCE, EMPTY_FA_ROI or TECHNICAL_ERROR is quarantined for the rest of the group (and excluded from final output). Other patients complete 601 references.
The first failed pair is recorded; ongoing waves skip that patient once prior wave Reduce commits.
A failure in FSL executable or reference-wide track loading can still cause a Map job to fail; no blanket suppression of infrastructure failures.
The method still uses legacy two-step interpolation, with onepass fallback for empty ROI. Validate affected regions scientifically.
Run a 20-ref pilot before allowing the auto controller to continue if possible.
