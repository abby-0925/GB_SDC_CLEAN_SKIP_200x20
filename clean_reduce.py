#!/usr/bin/env python3
"""Patient-level quarantine reducer: failed patient never contributes to final 601 denominator."""
import argparse, os, tempfile
from pathlib import Path
import numpy as np
import nibabel as nib
p=argparse.ArgumentParser();p.add_argument('--base',default='/public/home/tanghui/GB_SDC_TrackVis');p.add_argument('--group',required=True);p.add_argument('--wave',type=int,required=True);p.add_argument('--chunk',type=int,default=50)
a=p.parse_args();b=Path(a.base);root=b/'ref_major_clean_v1'
pats=[x.strip() for x in (root/'manifests'/f'{a.group}.patients.txt').read_text().splitlines() if x.strip()]
refs=[x.strip() for x in (root/'manifests/references.txt').read_text().splitlines() if x.strip()]
assert len(refs)==601 and len(pats)==len(set(pats)) and 1<=len(pats)<=200
r0=a.wave*20;r1=min(601,r0+20)
assert 0<=r0<r1
hroot=Path('/public/home/tanghui/PRECISE/DIFF601')
template=nib.load(str(hroot/refs[0]/f'{refs[0]}_T1_to_MNI2mm.nii.gz'))
state=root/'states'/a.group;state.mkdir(parents=True,exist_ok=True)
def atomic_npz(path,**kwargs):
 tmp=path.with_name('.'+path.stem+'.'+str(os.getpid())+'.npz')
 try:
  np.savez_compressed(tmp,**kwargs);os.replace(tmp,path)
 finally:
  if tmp.exists():tmp.unlink()
def atomic_text(path,body):
 tmp=path.with_name('.'+path.name+'.'+str(os.getpid())+'.tmp')
 tmp.write_text(body);os.replace(tmp,path)
# Never change deferred list while wave map tasks are running; update only after all maps complete.
prior_deferred=state/'deferred_patients.tsv'
excluded=set()
if prior_deferred.exists():
 for line in prior_deferred.read_text().splitlines()[1:]:
  if line.strip():excluded.add(line.split('\t')[0])
new_fail=[]
for ci,first in enumerate(range(0,len(pats),a.chunk)):
 subset=pats[first:first+a.chunk];path=state/f'c{ci:03d}.npz';marker=state/f'c{ci:03d}_wave{a.wave:03d}.ok'
 if marker.exists():
  with np.load(path,allow_pickle=False) as z:
   assert int(z['next_ref'])==r1 and list(z['patients'])==subset
  print('[REDUCE_SKIP]',ci,flush=True);continue
 if path.exists():
  with np.load(path,allow_pickle=False) as z:
   acc=z['counts'].copy();prior_bad=np.asarray(z['excluded'],dtype=bool).copy();done=int(z['next_ref']);saved=list(z['patients'])
  if done==r1 and saved==subset:
   atomic_text(marker,f'{r0}\t{r1}\n');continue
  assert done==r0 and saved==subset,(ci,done,r0)
 else:
  assert r0==0,f'Missing previous checkpoint c{ci:03d}'
  acc=np.zeros((len(subset),)+template.shape,dtype=np.uint16);prior_bad=np.zeros(len(subset),dtype=bool)
 assert acc.shape==(len(subset),)+template.shape
 bad=prior_bad.copy()
 for i,pat in enumerate(subset):
  if pat in excluded:bad[i]=True
 for ri in range(r0,r1):
  loc=root/'maps'/a.group/f'r{ri:03d}'/f'c{ci:03d}.npz'
  if not loc.is_file() or not loc.with_suffix('.ok').is_file():raise RuntimeError('Missing Map: '+str(loc))
  with np.load(loc,allow_pickle=False) as z:
   if list(z['patients'])!=subset or str(z['ref'])!=refs[ri]:raise RuntimeError('ID mismatch: '+str(loc))
   masks=np.asarray(z['mask'],dtype=np.uint8)
   statuses=[str(x) for x in z['status']]
   assert masks.shape==acc.shape and len(statuses)==len(subset)
   for i,(pat,st) in enumerate(zip(subset,statuses)):
    if st not in ('SUCCESS_NONZERO','SUCCESS_ZERO','DEFERRED_PATIENT'):
     bad[i]=True;new_fail.append((pat,refs[ri],f'wave{a.wave:03d}',st))
   for i,st in enumerate(statuses):
    if not bad[i] and st in ('SUCCESS_NONZERO','SUCCESS_ZERO'):
     acc[i]+=masks[i]
 # A patient that fails later in the same wave loses all prior counts, ensuring no partial SDC release.
 acc[bad]=0
 atomic_npz(path,counts=acc,patients=np.asarray(subset),next_ref=np.asarray(r1),excluded=bad)
 atomic_text(marker,f'{r0}\t{r1}\n')
 for i,flag in enumerate(bad):
  if flag:excluded.add(subset[i])
 print('[REDUCED]',ci,r1,'deferred',int(bad.sum()),flush=True)
# After all chunks committed, update immutable-for-next-wave quarantine list.
# Rebuild list from current state, so restart after partial Reduce is safe.
all_bad=set()
for ci,first in enumerate(range(0,len(pats),a.chunk)):
 with np.load(state/f'c{ci:03d}.npz',allow_pickle=False) as z:
  assert int(z['next_ref'])==r1
  all_bad.update(str(pid) for pid,flag in zip(z['patients'],z['excluded']) if flag)
old_rows={}
if prior_deferred.exists():
 for line in prior_deferred.read_text().splitlines()[1:]:
  if line.strip(): old_rows[line.split('\t')[0]]=line
for pat,ref,wave,st in new_fail:
 old_rows.setdefault(pat,f'{pat}\t{wave}\t{ref}\t{st}')
for pat in all_bad:old_rows.setdefault(pat,f'{pat}\tprevious\t-\tDEFERRED_PATIENT')
atomic_text(prior_deferred,'patient\tfirst_wave\treference\tstatus\n'+''.join(old_rows[k]+'\n' for k in sorted(old_rows)))
if r1==601:
 completed=0
 for ci,first in enumerate(range(0,len(pats),a.chunk)):
  with np.load(state/f'c{ci:03d}.npz',allow_pickle=False) as z:
   assert int(z['next_ref'])==601
   assert list(z['patients'])==pats[first:first+a.chunk]
   for i,pat in enumerate(z['patients']):
    pat=str(pat)
    if bool(z['excluded'][i]):continue
    arr=z['counts'][i]
    dest=root/'final'/pat;dest.mkdir(parents=True,exist_ok=True)
    prob=arr.astype(np.float32)/601.0
    for n,data,dtype in [('SDC_reference_count',arr,np.uint16),('SDC_probability',prob,np.float32),('SDC_binary_p50',(prob>=.5).astype(np.uint8),np.uint8)]:
     im=nib.Nifti1Image(np.asarray(data),template.affine,header=template.header.copy());im.set_data_dtype(dtype)
     tmp=dest/f'.{n}.part.nii.gz';nib.save(im,str(tmp));os.replace(tmp,dest/f'{n}.nii.gz')
    atomic_text(dest/'MERGE_SUCCESS.txt',f'MERGE_SUCCESS {pat} 601\n');completed+=1
 atomic_text(state/'GROUP_DONE.txt',f'{a.group} 601 completed={completed} deferred={len(all_bad)}\n')
 print('[FINAL]',completed,'patients; deferred',len(all_bad),flush=True)
print('[WAVE_DONE]',a.group,a.wave,'deferred',len(all_bad),flush=True)
