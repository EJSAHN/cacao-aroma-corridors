"""A tiny native-tool software test, isolated from every biological input."""
from __future__ import annotations
import hashlib,json,random,tempfile
from pathlib import Path
from .common import DataError,sha256
from .tags import Tag,revcomp,write_fasta
from .sam import parse_sam
from .toolchain import reference_cache,Reference,build_index,align_cached,digest_json


def native_smoke(tools,cache,threads,env,cfg,log):
    key=digest_json({'binary_hashes':tools['sha256'],'policy':cfg['alignment'],'smoke_fixture':2,'threads':threads})
    base=cache/'native_smoke_test';base.mkdir(exist_ok=True)
    report=base/('PASS_'+key+'.json')
    if report.is_file():
        obj=json.loads(report.read_text(encoding='utf-8'))
        if obj.get('test_key')==key and obj.get('status')=='PASS':log.info('Native smoke test already passed for these executable hashes and settings.');return obj
    rng=random.Random(612924);dna=lambda n:''.join(rng.choice('ACGT') for _ in range(n))
    tags={};genome=dna(100)
    expected={};expected_positions={}
    for i,kind in enumerate(['exact','reverse','distal_mismatch','distal_insertion','distal_deletion','duplicated','unmapped']):
        seq=dna(69) if kind!='unmapped' else 'A'*69
        altbase=next(b for b in 'ACGT' if b!=seq[33]);tag=Tag('smoke_'+kind,i+2,kind,kind,seq,seq[:33]+altbase+seq[34:],33);tags[tag.query_id]=tag
        refseq=revcomp(seq) if kind=='reverse' else seq
        if kind=='distal_mismatch':refseq=seq[:8]+next(b for b in 'ACGT' if b!=seq[8])+seq[9:]
        if kind=='distal_insertion':refseq=seq[:8]+seq[9:]
        if kind=='distal_deletion':refseq=seq[:8]+'G'+seq[8:]
        if kind!='unmapped':
            focal_index=(len(seq)-1-tag.snp0) if kind=='reverse' else tag.snp0
            if kind=='distal_insertion':focal_index-=1
            if kind=='distal_deletion':focal_index+=1
            if kind!='duplicated':expected_positions[tag.query_id]=(len(genome)+focal_index+1, '-' if kind=='reverse' else '+')
            genome+=refseq+dna(150)
            if kind=='duplicated':genome+=refseq+dna(150)
        expected[tag.query_id]=kind not in ('duplicated','unmapped')
    log.info('Native smoke test: seven known synthetic software cases; not cacao data.')
    with tempfile.TemporaryDirectory(prefix='inputs_',dir=base) as d:
        src=Path(d)/'smoke.fa';src.write_text('>smoke_contig\n'+genome+'\n',encoding='ascii')
        inv=[{'sequence_id':'smoke_contig','chromosome':'test','length_bp':len(genome),'sequence_sha256':hashlib.sha256(genome.encode()).hexdigest()}]
        fa,ix=reference_cache(src,{'sha256':sha256(src)},inv,base,log)
        q=Path(d)/'tags.fa';write_fasta(q,tags)
        idx=build_index(fa,base,tools,threads,env,cfg,log);sam,_=align_cached(q,idx,base,tools,threads,env,cfg,log)
        with Reference(fa,ix) as ref:rows,_,_,_=parse_sam(sam,tags,ref,cfg['alignment'],{'smoke_contig':'test'},log)
        actual={r['query_id']:bool(r['accepted_for_mapping']) for r in rows}
        if actual!=expected:raise DataError('Native tool/software smoke test did not match known expectations. No cacao alignment was started. '+json.dumps(rows))
        for row in rows:
            if row['query_id'] in expected_positions and (row.get('mapped_snp_bp'),row.get('mapped_strand'))!=expected_positions[row['query_id']]:
                raise DataError('Native smoke-test SNP coordinate or strand is incorrect. No cacao alignment was started.')
    obj={'test_key':key,'status':'PASS','purpose':'Synthetic software verification only; excluded from cacao scientific outputs.','cases':actual,'expected_retained_snp_positions_and_strands':expected_positions}
    report.write_text(json.dumps(obj,indent=2)+'\n', encoding='utf-8');log.info('Native smoke test PASS (seven cases).');return obj
