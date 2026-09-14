"""Pinned native Bowtie 2, reference and alignment caches; no environment installs.

All new bytes are written under caller-selected tool/cache directories. Downloads
use TLS and a release-asset SHA-256. Existing executables are version-checked and
hashed; their digest is provenance, not a claim of a verified publisher binary.
"""
from __future__ import annotations
import gzip
import hashlib
import json
import mmap
import os
from pathlib import Path,PurePosixPath
import platform
import re
import shutil
import stat
import subprocess
import tempfile
import time
import urllib.request
import zipfile
from .common import DataError,contains,sha256

SUFFIXES=['.1.bt2','.2.bt2','.3.bt2','.4.bt2','.rev.1.bt2','.rev.2.bt2']

def digest_json(obj):return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()

def free_check(path,gib):
    if shutil.disk_usage(path).free<gib*2**30:raise DataError('Insufficient free space on cache/output volume: '+str(path))

def version_of(path,required):
    if not path.is_file():raise DataError('Native executable missing: '+str(path))
    try:r=subprocess.run([str(path),'--version'],capture_output=True,text=True,timeout=30,check=True)
    except (OSError,subprocess.SubprocessError) as exc:raise DataError('Native executable does not run: '+str(path)+'; '+str(exc)) from exc
    output=r.stdout+'\n'+r.stderr
    if not re.search(r'\bversion\s+'+re.escape(required)+r'(?:\s|$)',output,re.I):raise DataError('Bowtie 2 version is not the pinned '+required+': '+str(path))
    return output.strip()


def check_tool_dir(folder,version):
    ext='.exe' if os.name=='nt' else ''
    a=folder/('bowtie2-align-s'+ext);b=folder/('bowtie2-build-s'+ext)
    versions={x.name:version_of(x,version) for x in (a,b)}
    return {'aligner':a,'builder':b,'versions':versions,'sha256':{x.name:sha256(x) for x in (a,b)}}


def safe_tool_extract(archive,dest):
    with zipfile.ZipFile(archive) as z:
        total=0;seen=set()
        for m in z.infolist():
            p=PurePosixPath(m.filename);parts=m.filename.rstrip('/').split('/')
            if p.is_absolute() or '\\' in m.filename or any(x in ('','..','.') or ':' in x for x in parts):raise DataError('Unsafe native-tool ZIP member.')
            if stat.S_IFMT(m.external_attr>>16) not in (0,stat.S_IFREG,stat.S_IFDIR):raise DataError('Links/special entries in native-tool ZIP.')
            if m.filename.casefold() in seen:raise DataError('Duplicate native-tool ZIP member.')
            seen.add(m.filename.casefold());total+=m.file_size
            if total>1024**3:raise DataError('Native-tool archive exceeds extraction bound.')
            target=dest.joinpath(*p.parts)
            if not contains(dest,target):raise DataError('Tool extraction path escapes destination.')
            if m.is_dir():target.mkdir(parents=True,exist_ok=True)
            else:
                target.parent.mkdir(parents=True,exist_ok=True)
                with z.open(m) as src,target.open('xb') as out:shutil.copyfileobj(src,out,1024*1024)
                if os.name!='nt' and (m.external_attr>>16)&0o111:target.chmod(0o755)


def resolve_tools(toolroot,cfg,explicit,offline,log):
    ver=cfg['alignment']['version'];spec=cfg['tools'];ext='.exe' if os.name=='nt' else ''
    if explicit is not None:
        obj=check_tool_dir(explicit.resolve(),ver);obj['origin']='explicit_existing_directory';return obj
    probes=[]
    installed=toolroot/('bowtie2-'+ver)
    if installed.is_dir():probes.extend(p.parent for p in installed.rglob('bowtie2-align-s'+ext))
    x=shutil.which('bowtie2-align-s'+ext)
    if x:probes.append(Path(x).parent)
    seen=set()
    for p in probes:
        if p in seen:continue
        seen.add(p)
        if not (p/('bowtie2-align-s'+ext)).is_file():continue
        try:
            obj=check_tool_dir(p.resolve(),ver);obj['origin']='existing_version_checked_directory';log.info('Reusing native Bowtie 2: %s',p);return obj
        except DataError as e:log.warning('Not using executable candidate: %s',e)
    if offline:raise DataError('Pinned native Bowtie 2 not found. Supply --bowtie2-dir or run without --offline to download to configured tools.')
    machine=platform.machine().lower();machine='x86_64' if machine in ('amd64','x86_64') else machine
    label=platform.system()+'-'+machine
    if label not in spec['platforms']:raise DataError('No automatic native binary for '+label+'; supply --bowtie2-dir.')
    asset=spec['platforms'][label];toolroot.mkdir(parents=True,exist_ok=True);free_check(toolroot,1.0)
    target=toolroot/asset['asset'];url=spec['base_url']+asset['asset']
    if target.exists() and sha256(target)!=asset['sha256']:raise DataError('Existing native-tool archive hash differs; it was not overwritten: '+str(target))
    if not target.exists():
        tmp=toolroot/(asset['asset']+'.partial_'+str(os.getpid()))
        log.info('Downloading native Bowtie 2 %.1f MiB to %s',asset['bytes']/2**20,target)
        try:
            req=urllib.request.Request(url,headers={'User-Agent':'CacaoMarkerHarmonization/3.0.0'})
            with urllib.request.urlopen(req,timeout=90) as response,tmp.open('xb') as out:
                if not response.geturl().startswith('https://'):raise DataError('Non-HTTPS tool redirect.')
                n=0;last=0
                while True:
                    chunk=response.read(1024*1024)
                    if not chunk:break
                    n+=len(chunk)
                    if n>asset['bytes']:raise DataError('Tool download larger than pinned release asset.')
                    out.write(chunk)
                    if n-last>=10*2**20:log.info('Tool download: %.1f MiB',n/2**20);last=n
            if tmp.stat().st_size!=asset['bytes'] or sha256(tmp)!=asset['sha256']:raise DataError('Native-tool download failed pinned SHA-256/size verification.')
            if target.exists():raise DataError('Native-tool archive appeared concurrently.')
            tmp.rename(target)
        except Exception as exc:
            raise DataError('Native Bowtie 2 download failed; no environment was changed. '+str(exc)+' URL: '+url) from exc
        finally:
            if tmp.exists():tmp.unlink()
    if installed.exists():raise DataError('Tool cache exists but no valid pinned executable was found. Preserve it and supply --bowtie2-dir.')
    with tempfile.TemporaryDirectory(prefix='.tool_',dir=toolroot) as temp:
        stage=Path(temp);safe_tool_extract(target,stage)
        aligns=list(stage.rglob('bowtie2-align-s'+ext))
        if len(aligns)!=1:raise DataError('Native-tool ZIP has an unexpected executable layout.')
        check_tool_dir(aligns[0].parent,ver)
        stage.rename(installed)
    found=next(installed.rglob('bowtie2-align-s'+ext));obj=check_tool_dir(found.parent,ver)
    obj['origin']='downloaded_pinned_release_sha256';obj['asset_sha256']=asset['sha256'];return obj


def execute(cmd,logfile,env,cwd,timeout,log):
    logfile.parent.mkdir(parents=True,exist_ok=True)
    log.info('External command: %s',subprocess.list2cmdline([str(c) for c in cmd]))
    start=time.monotonic();proc=None
    try:
        with logfile.open('xb') as f:
            proc=subprocess.Popen([str(c) for c in cmd],cwd=str(cwd),env=env,stdout=f,stderr=subprocess.STDOUT)
            last=start
            while proc.poll() is None:
                now=time.monotonic()
                if now-start>timeout:raise DataError('External command timed out; see '+str(logfile))
                if now-last>=15:
                    log.info('Native process running: %.0f seconds; log %.1f KiB',now-start,logfile.stat().st_size/1024);last=now
                    free_check(cwd,0.75)
                time.sleep(0.3)
            if proc.returncode:
                tail=logfile.read_bytes()[-6000:].decode('utf-8',errors='replace')
                raise DataError('Native command failed (exit '+str(proc.returncode)+').\n'+tail)
    finally:
        if proc is not None and proc.poll() is None:proc.kill();proc.wait()


def validate_cache(folder,expected):
    meta=folder/'CACHE.json'
    if not meta.is_file():raise DataError('Incomplete cache directory (not reused): '+str(folder))
    obj=json.loads(meta.read_text(encoding='utf-8'))
    if obj['identity']!=expected:raise DataError('Cache identity mismatch.')
    for name,h in obj['files'].items():
        if Path(name).name!=name or not (folder/name).is_file() or (folder/name).is_symlink() or sha256(folder/name)!=h:raise DataError('Cache file is missing/changed: '+str(folder/name))
    return obj


def finish_cache(stage,dest,identity):
    files={p.name:sha256(p) for p in sorted(stage.iterdir()) if p.is_file()}
    (stage/'CACHE.json').write_text(json.dumps({'identity':identity,'files':files},indent=2)+'\n', encoding='utf-8')
    if dest.exists():raise DataError('Cache destination appeared concurrently.')
    stage.rename(dest)




def reference_cache(source,asset,inventory,cache,log):
    identity={'kind':'reference_fasta','source_sha256':asset['sha256'],'inventory':digest_json(inventory),'layout':'unwrapped-v1'}
    dest=cache/('reference_'+digest_json(identity)[:20]);cache.mkdir(parents=True,exist_ok=True)
    if dest.exists():validate_cache(dest,identity);return dest/'genome.fa',dest/'sequences.json'
    expected={r['sequence_id']:r for r in inventory}
    if len(expected)!=len(inventory):raise DataError('Duplicate reference sequence ID.')
    with tempfile.TemporaryDirectory(prefix='.reference_',dir=cache) as temp:
        stage=Path(temp);entries={};name=None;seq=bytearray()
        with (stage/'genome.fa').open('xb') as out:
            def emit():
                if name is None:return
                r=expected.get(name)
                if r is None or name in entries:raise DataError('Unexpected/duplicate FASTA sequence.')
                if len(seq)!=r['length_bp'] or hashlib.sha256(seq).hexdigest()!=r['sequence_sha256']:raise DataError('FASTA sequence differs from the checked exact-tag inventory: '+name)
                if not re.fullmatch(r'[A-Za-z0-9_.-]+',name):raise DataError('Unsupported FASTA identifier for SAM.')
                out.write(('>'+name+'\n').encode());offset=out.tell();out.write(seq);out.write(b'\n')
                entries[name]={'offset':offset,'length':len(seq),'chromosome':r.get('chromosome'),'sequence_sha256':r['sequence_sha256']}
            op=gzip.open if source.suffix=='.gz' else open
            with op(source,'rt',encoding='utf-8-sig') as f:
                for line in f:
                    line=line.strip()
                    if not line:continue
                    if line.startswith('>'):emit();name=line[1:].split()[0];seq=bytearray()
                    else:
                        if name is None or not re.fullmatch('[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+',line):raise DataError('Malformed FASTA sequence.')
                        seq.extend(line.upper().encode('ascii'))
                        if len(seq)>200000000:raise DataError('Reference contig exceeds 200 Mb bounded read.')
                emit()
        if set(entries)!=set(expected):raise DataError('Incomplete reference FASTA.')
        (stage/'sequences.json').write_text(json.dumps(entries,indent=2)+'\n', encoding='utf-8');finish_cache(stage,dest,identity)
    log.info('Validated reference prepared in configured cache: %s',dest)
    return dest/'genome.fa',dest/'sequences.json'


class Reference:
    def __init__(self,fasta,index):
        self.entries=json.loads(index.read_text(encoding='utf-8'));self.f=fasta.open('rb');self.m=mmap.mmap(self.f.fileno(),0,access=mmap.ACCESS_READ)
    def slice(self,seqid,start,end):
        r=self.entries.get(seqid)
        if r is None or start<0 or end<start or end>r['length']:raise DataError('Reference slice outside sequence bounds.')
        return self.m[r['offset']+start:r['offset']+end]
    def close(self):self.m.close();self.f.close()
    def __enter__(self):return self
    def __exit__(self,*args):self.close()


def build_index(fasta,cache,tools,threads,env,cfg,log):
    if ',' in str(fasta):raise DataError('Bowtie 2 FASTA paths cannot contain commas.')
    ident={'kind':'bowtie2_index','fasta_sha256':sha256(fasta),'builder_sha256':tools['sha256'][tools['builder'].name],
           'version':cfg['alignment']['version'],'threads':threads,'offrate':5}
    dest=cache/('index_'+digest_json(ident)[:20])
    if dest.exists():validate_cache(dest,ident);log.info('Reusing Bowtie 2 index: %s',dest);return dest/'genome'
    with tempfile.TemporaryDirectory(prefix='.index_',dir=cache) as temp:
        stage=Path(temp);prefix=stage/'genome'
        cmd=[tools['builder'],'--threads',str(threads),'--offrate','5',str(fasta),str(prefix)]
        execute(cmd,stage/'build.log',env,stage,cfg['alignment']['timeout_seconds'],log)
        if not all((stage/('genome'+s)).is_file() and (stage/('genome'+s)).stat().st_size>0 for s in SUFFIXES):raise DataError('Bowtie 2 did not produce the complete small-index file set.')
        finish_cache(stage,dest,ident)
    return dest/'genome'


def align_cached(query_fasta,prefix,cache,tools,threads,env,cfg,log):
    a=cfg['alignment']
    params=['--wrapper','basic-0','--end-to-end','--very-sensitive','-f','--ignore-quals',
            '-D',str(a['effort_D']),'-R',str(a['effort_R']),'-N',str(a['seed_mismatches']),'-L',str(a['seed_length']),
            '-i',a['seed_interval'],'--score-min',a['score_min'],'--mp',f"{a['mismatch_penalty']},{a['mismatch_penalty']}",
            '--rdg',f"{a['gap_open']},{a['gap_extend']}",'--rfg',f"{a['gap_open']},{a['gap_extend']}",
            '--np','1','--gbar','4','-k',str(a['max_reported_hits']),'--seed',str(a['seed']),'-p',str(threads),'--reorder']
    ident={'kind':'bowtie2_sam','query_sha256':sha256(query_fasta),'index_manifest_sha256':sha256(prefix.parent/'CACHE.json'),
           'aligner_sha256':tools['sha256'][tools['aligner'].name],'parameters':params}
    dest=cache/('alignment_'+digest_json(ident)[:20])
    if dest.exists():validate_cache(dest,ident);log.info('Reusing verified compressed SAM: %s',dest);return dest/'tags.sam.gz',ident
    if ',' in str(query_fasta):raise DataError('Bowtie 2 query paths cannot contain commas.')
    with tempfile.TemporaryDirectory(prefix='.alignment_',dir=cache) as temp:
        stage=Path(temp);sam=stage/'tags.sam'
        cmd=[tools['aligner'],*params,'-x',str(prefix),'-U',str(query_fasta),'-S',str(sam)]
        execute(cmd,stage/'alignment.log',env,stage,a['timeout_seconds'],log)
        if not sam.is_file() or sam.stat().st_size==0:raise DataError('No SAM output from native aligner.')
        with sam.open('rb') as src,(stage/'tags.sam.gz').open('xb') as raw,gzip.GzipFile(fileobj=raw,mode='wb',mtime=0,filename='') as out:shutil.copyfileobj(src,out,1024*1024)
        sam.unlink() # This stage owns the uncompressed SAM; originals are untouched.
        (stage/'COMMAND.json').write_text(json.dumps([str(x) for x in cmd],indent=2)+'\n', encoding='utf-8')
        finish_cache(stage,dest,ident)
    return dest/'tags.sam.gz',ident
