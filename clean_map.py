#!/usr/bin/env python3
"""One healthy reference per process, 50 patients per independently resumable chunk.
Uses the same fiber voxelization and FSL transformation helpers as the existing FAST prototype.
"""
import argparse,os,sys,tempfile,time,json
from pathlib import Path
import numpy as np
import nibabel as nib
p=argparse.ArgumentParser();p.add_argument('--base',default='/public/home/tanghui/GB_SDC_TrackVis');p.add_argument('--group',required=True);p.add_argument('--reference-index',type=int,required=True);p.add_argument('--patient-chunk',type=int,default=50);p.add_argument('--sample-step',type=float,default=.45);p.add_argument('--retry-failed',action='store_true');a=p.parse_args()
b=Path(a.base);sys.path.insert(0,str(b/'scripts'));import fast_batch as fast
root=b/'ref_major_clean_v1'; pats=[x.strip() for x in (root/'manifests'/f'{a.group}.patients.txt').read_text().splitlines() if x.strip()];refs=[x.strip() for x in (root/'manifests/references.txt').read_text().splitlines() if x.strip()]
if not 0<=a.reference_index<len(refs):raise SystemExit('bad reference index')
ref=refs[a.reference_index]; d=fast.HROOT/ref; fa=d/f'{ref}_FA.nii.gz';tr=d/f'{ref}_track_s.trk';t1=d/f'{ref}_ses-001_T1w_biascorr_brain.nii.gz';t1fa=d/f'{ref}_T1_to_FA_mat';fat1=d/f'{ref}_FA_to_T1_mat';fwd=d/f'{ref}_T1_to_MNI2mm_warp.nii.gz';inv=b/'shared_invwarp'/ref/'invwarp_T1.nii.gz';template=fast.HROOT/refs[0]/f'{refs[0]}_T1_to_MNI2mm.nii.gz'
for f in [fa,tr,t1,t1fa,fat1,fwd,inv,template]:
 if not f.is_file() or not f.stat().st_size:raise RuntimeError(f'Missing {f}')
target=root/'maps'/a.group/f'r{a.reference_index:03d}';target.mkdir(parents=True,exist_ok=True)
# immutable patient chunk size per group
meta=target/'config.json'; conf={'group':a.group,'ref_index':a.reference_index,'ref':ref,'chunk':a.patient_chunk,'patients':pats,'sample_step':a.sample_step}
if meta.exists() and json.loads(meta.read_text())!=conf:raise RuntimeError('Changed map configuration')
if not meta.exists():meta.write_text(json.dumps(conf))
faimg=nib.load(str(fa));tmpl=nib.load(str(template));assert tmpl.ndim==3
inc,nfib,diff=fast.build_index(tr,faimg,a.sample_step);print(f'[INDEX] {ref} fibers={nfib} nnz={inc.nnz} difference_mm={diff:.5f}',flush=True)
for ci,first in enumerate(range(0,len(pats),a.patient_chunk)):
 deferred_file=root/'states'/a.group/'deferred_patients.tsv'
 deferred=set()
 if deferred_file.exists():
  for line in deferred_file.read_text().splitlines()[1:]:
   if line.strip(): deferred.add(line.split('\t')[0])
 out=target/f'c{ci:03d}.npz'; ok=target/f'c{ci:03d}.ok';subset=pats[first:first+a.patient_chunk]
 # A committed chunk is immutable unless retry-failed is explicitly selected.
 old_masks=None; old_status=None; old_selected=None
 if ok.exists() and out.exists() and out.stat().st_size:
  with np.load(out,allow_pickle=False) as z:
   if list(z['patients'])!=subset or str(z['ref'])!=ref or z['mask'].shape!=(len(subset),)+tmpl.shape:
    raise RuntimeError('Committed chunk mismatch: '+str(out))
   if 'status' not in z.files:
    raise RuntimeError('Legacy chunk without status; migrate or use a NEW group: '+str(out))
   old_status=np.asarray(z['status']).astype('U32')
   old_masks=z['mask'].copy()
   old_selected=z['selected'].copy()
  if np.all(np.isin(old_status,['SUCCESS_NONZERO','SUCCESS_ZERO'])):
   print('[SKIP_VALID]',ref,ci,flush=True);continue
  if not a.retry_failed:
   print('[SKIP_HAS_FAILURES]',ref,ci,'failed=',int(np.sum(~np.isin(old_status,['SUCCESS_NONZERO','SUCCESS_ZERO']))),flush=True)
   continue
 elif ok.exists() != out.exists():
  raise RuntimeError('Broken committed chunk: '+str(out))
 masks=old_masks if old_masks is not None else np.zeros((len(subset),)+tmpl.shape,dtype=np.uint8)
 selected_arr=old_selected if old_selected is not None else np.zeros(len(subset),dtype=np.int32)
 status=old_status if old_status is not None else np.full(len(subset),'PENDING',dtype='U32')
 issues=[];t0=time.time()
 with tempfile.TemporaryDirectory(prefix='sdc_ref_',dir=os.environ.get('SLURM_TMPDIR') if os.path.isdir(os.environ.get('SLURM_TMPDIR','')) else None) as td:
  tmp=Path(td)
  for pi,pat in enumerate(subset):
   if pat in deferred:
    status[pi]='DEFERRED_PATIENT';masks[pi].fill(0);selected_arr[pi]=0;continue
   if status[pi] in ('SUCCESS_NONZERO','SUCCESS_ZERO'):continue
   # Clear stale counts before retrying a failed combination.
   masks[pi].fill(0);selected_arr[pi]=0;status[pi]='PENDING'
   try:
    lesion=fast.LROOT/pat/f'{pat}_lesion_mask_to_MNI152_2mm.nii.gz'
    le=nib.load(str(lesion))
    if le.shape!=tmpl.shape or not np.allclose(le.affine,tmpl.affine,atol=1e-3):raise RuntimeError('MNI mismatch')
    if not np.any(np.asarray(le.dataobj)>0):
     status[pi]='EMPTY_SOURCE';continue
    lt=tmp/'lesion_T1.nii.gz';lf=tmp/'lesion_FA.nii.gz';fm=tmp/'disconnect_FA.nii.gz';mm=tmp/'disconnect_MNI.nii.gz'
    fast.run(['/public/apps/neuroscience/fsl/6.0.4/bin/applywarp','--in='+str(lesion),'--ref='+str(t1),'--warp='+str(inv),'--interp=nn','--out='+str(lt)])
    fast.run(['/public/apps/neuroscience/fsl/6.0.4/bin/flirt','-in',str(lt),'-ref',str(fa),'-applyxfm','-init',str(t1fa),'-interp','nearestneighbour','-out',str(lf)])
    mask=nib.load(str(lf)).get_fdata(dtype=np.float32)>0.5
    if mask.shape!=faimg.shape[:3]:raise RuntimeError('FA grid mismatch')
    if not np.any(mask):
     # Recover tiny lesions lost during two successive nearest-neighbour resamplings.
     # Warp is MNI->T1, postmat is T1->FA, evaluated once on FA grid.
     onepass=tmp/'lesion_FA_onepass.nii.gz'
     fast.run(['/public/apps/neuroscience/fsl/6.0.4/bin/applywarp',
       '--in='+str(lesion),'--ref='+str(fa),'--warp='+str(inv),
       '--postmat='+str(t1fa),'--interp=nn','--out='+str(onepass)])
     recovered=nib.load(str(onepass))
     if recovered.shape[:3]!=faimg.shape[:3] or not np.allclose(recovered.affine,faimg.affine,atol=1e-4):
      raise RuntimeError('One-pass FA geometry mismatch')
     mask=np.asarray(recovered.dataobj)>0.5
     if not np.any(mask):
      status[pi]='EMPTY_FA_ROI';continue
     print('[ONEPASS_RECOVERED]',pat,ref,'voxels='+str(int(mask.sum())),flush=True)
    vox,selected=fast.select_full_fibers(inc,mask)
    selected_arr[pi]=int(selected)
    if selected:
     ni=nib.Nifti1Image(vox.reshape(faimg.shape[:3]),faimg.affine,header=faimg.header.copy());ni.set_data_dtype(np.uint8);nib.save(ni,str(fm))
     fast.run(['/public/apps/neuroscience/fsl/6.0.4/bin/applywarp','--in='+str(fm),'--ref='+str(template),'--warp='+str(fwd),'--premat='+str(fat1),'--interp=nn','--out='+str(mm)])
     im=nib.load(str(mm))
     if im.shape!=tmpl.shape or not np.allclose(im.affine,tmpl.affine,atol=1e-3):raise RuntimeError('Output geometry mismatch')
     masks[pi]=np.asarray(im.dataobj)>0
     status[pi]='SUCCESS_NONZERO'
    else:status[pi]='SUCCESS_ZERO'
   except Exception as exc:
    status[pi]='TECHNICAL_ERROR'
    issues.append((pat,type(exc).__name__,str(exc).replace('\n',' ')[:500]))
   if status[pi] not in ('SUCCESS_NONZERO','SUCCESS_ZERO'):
    print('[DEFER]',ref,pat,status[pi],flush=True)
 # Publish complete chunk with per-patient validity, even when some fail.
 tmpout=target/f'.c{ci:03d}.{os.getpid()}.npz'
 np.savez_compressed(tmpout,mask=masks,patients=np.array(subset),ref=np.array(ref),selected=selected_arr,status=status)
 os.replace(tmpout,out)
 # Write optional error details, then marker LAST; marker is not full-success evidence.
 issue_path=target/f'c{ci:03d}.issues.tsv'
 issue_path.write_text('patient\ttype\tmessage\n'+''.join(f'{x}\t{y}\t{z}\n' for x,y,z in issues))
 ok.write_text(f'{ref}\t{ci}\t{len(subset)}\tvalid={int(np.sum(np.isin(status,["SUCCESS_NONZERO","SUCCESS_ZERO"]))) }\n')
 print(f'[CHUNK_DONE] ref={ref} chunk={ci} n={len(subset)} valid={int(np.sum(np.isin(status,["SUCCESS_NONZERO","SUCCESS_ZERO"]))) } failed={int(np.sum(~np.isin(status,["SUCCESS_NONZERO","SUCCESS_ZERO"]))) } seconds={time.time()-t0:.1f}',flush=True)
print('[REFERENCE_DONE]',a.group,ref,flush=True)
