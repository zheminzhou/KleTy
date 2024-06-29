import tempfile, numpy as np, pandas as pd
import click, os, logging
from _collections import OrderedDict


try :
    from .configure import dbs, scheme_info
    from .uberBlast import uberBlast, readFastq, rc
except :
    from configure import dbs, scheme_info
    from uberBlast import uberBlast, readFastq, rc


logging.basicConfig(format='%(asctime)s - %(message)s', level=20)


refseqs = dbs['refseqs']
core_genes = dbs['core_genes']
repr = dbs['repr']
hcc = dbs['hcc']
species = dbs['species']



import hashlib
def get_md5(value, dtype=str):
    m = hashlib.md5(str(value).encode()).hexdigest()
    if dtype == str:
        return m
    else:
        return int(m, 16)





@click.command()
@click.option('-q', '--query', help='input assembly. fasta or fastq format, can be gziped.')
@click.option('-o', '--outfile', help='alleles output file')
@click.option('-n', '--n_thread', help='n_threads [default:8]', default=8, type=int)
def main(query, outfile, n_thread) :
    alleles, hiercc = cgmlst(query, n_thread)
    with open(outfile, 'wt') as fout :
        for gene, allele in sorted(alleles.items()) :
            fout.write('>{gene_name} value_md5={value_md5} CIGAR={CIGAR} accepted={flag} reference={reference} identity={identity:.3f} coordinates={coordinates}\n{sequence}\n'.format(**allele))
    return outfile


def cgmlst(query, n_thread=8) :
    alleles = nomenclature(query, n_thread)

    repr_profile = pd.read_parquet(repr)
    repr_hcc = pd.read_csv(hcc, sep='\t').set_index('#ST_id')
    repr_species = pd.read_csv(species, sep='\t').set_index('ID')
    
    profile = np.array([('-' if v.startswith('-') else v) for v in [v.get('value_md5', '-') for k, v in sorted(alleles.items())]])
    relshare = np.sum((repr_profile.values[:, 1:] == profile), 1)/np.max((np.sum((repr_profile.values[:, 1:] != '') & (profile != '-'), 1), np.sum((repr_profile.values[:, 1:] != ''), 1)*0.97), 0)*profile.size
    max_idx = np.argmax(relshare)
    min_dist = int(profile.size-relshare[max_idx]+0.5)
    ref_repr = repr_profile.values[max_idx, 0][:15]
    hc = repr_hcc.loc[ref_repr].values.astype(str).tolist()
    sp = repr_species.loc[ref_repr].values[0].replace(' ', '_') if min_dist < 2900 else 'Novel_species'
    hc[:min_dist] = ['ND'] * min_dist
    hiercc = {'reference':ref_repr, 'distance':min_dist, 'species':sp}
    hiercc['HC1360.500.200.100.50.20.10.5.2'] = '.'.join([hc[d] for d in (1360,500,200,100,50,20,10,5,2)])
    return alleles, hiercc
    
def nomenclature(query, n_thread=8) :
    core = {}
    with open(core_genes, 'rt') as fin :
        for line in fin :
            core[line.strip().split()[0]] = 1
            
    with tempfile.TemporaryDirectory(prefix='NM_', dir='.') as dirname :
        blastab = uberBlast(
            '-r {0} -q {1} -f --blastn --min_id {2} --min_ratio {3} -t {4} -p -s 2 -e 9,12 -m --merge_gap 300'.format(
                query, refseqs, scheme_info['min_iden']-0.05, 0.05, n_thread).split())
    merges = {}
    for b in blastab.T[16] :
        if len(b) > 4 :
            key = tuple(b[3:])
            merges[key] = b[:3]
    bsn = { b[15]:b[:15] for b in blastab }
    for ids, (score, iden, size) in merges.items() :
        bs = np.array([ bsn[i] for i in ids ], dtype=object)
        if np.unique(bs.T[1]).size == 1 :
            bs[0][2], bs[0][3], bs[0][7], bs[0][9], bs[0][11], bs[0][14] = iden, size, bs[-1][7], bs[-1][9], score, 'COMPLEX'
        for i in ids[1:] :
            bsn.pop(i)
        bsn[ids[0]] = bs[0]
    blastab = np.array(list(bsn.values()))

    blastab = blastab[(blastab.T[2] >= scheme_info['min_iden']) & ((blastab.T[7]-blastab.T[6]+1) >= scheme_info['min_frag'] * blastab.T[12])]
    blastab.T[11] = blastab.T[11]*(blastab.T[7] - blastab.T[6]+1)/blastab.T[12]
    blastab = blastab[np.lexsort((blastab[:, 8], blastab[:, 1]))]
    for i, b0 in enumerate(blastab[:-1]) :
        if b0[0] == '' :
            continue
        s0, e0 = sorted(b0[8:10])
        todel = []
        for b1 in blastab[i+1:] :
            s1, e1 = sorted(b1[8:10])
            if b0[1] != b1[1] or e0 < s1 :
                break
            ovl = min(e0, e1) - max(s0, s1) + 1
            if ovl >= 0.5 * (e0-s0+1) or ovl >= 0.5 * (e1-s1+1) :
                sc0, sc1 = abs(b0[11]), abs(b1[11])
                g0, g1 = b0[0].rsplit('_', 1)[0], b1[0].rsplit('_', 1)[0]
                if b0[2] < b1[2]*scheme_info['max_iden'] or (b1[2] >= b0[2]*scheme_info['max_iden']
                    and (sc0 < sc1 or (sc0 == sc1 and b0[0] > b1[0]))) :
                    b0[11] = -sc0
                    if g0 == g1 or sc0 < sc1 * scheme_info['max_iden'] or b0[2] < b1[2]*scheme_info['max_iden'] :
                        b0[0] = ''
                        break
                else :
                    b1[11] = -sc1
                    if g0 == g1 or sc1 < sc0 * scheme_info['max_iden'] or b1[2] < b0[2]*scheme_info['max_iden'] :
                        todel.append(b1)
        if b0[0] and len(todel) :
            for b1 in todel :
                b1[0] = ''
    blastab = blastab[blastab.T[0] != '']
    blastab = blastab[np.lexsort([-blastab.T[11], [b.rsplit('_', 1)[0] for b in blastab.T[0]]])]
    alleles = OrderedDict()
    for bsn in blastab :
        gene = bsn[0].rsplit('_', 1)[0]
        if gene in alleles :
            if alleles[gene]['score']*scheme_info['max_iden'] > bsn[11] :
                continue
            alleles[gene]['coordinates'].append((bsn[1], bsn[8], bsn[9]))
            alleles[gene]['flag'] |= 32
            continue
        flag = 0
        if bsn[6] > 1 or bsn[7] < bsn[12] :
            flag = 64
            if bsn[14] != 'COMPLEX' :
                if bsn[6] > 1 :
                    bsn[14] = '{0}D{1}'.format(bsn[6]-1, bsn[14])
                if bsn[7] < bsn[12] :
                    bsn[14] = '{0}{1}D'.format(bsn[14], bsn[12]-bsn[7])
        alleles[gene] = {'gene_name': gene, 'CIGAR':bsn[0]+':'+bsn[14], 'reference':os.path.basename(query),
                         'identity': bsn[2], 'coordinates':[(bsn[1], bsn[8], bsn[9])], 'flag':flag, 'score':bsn[11]}

    seq, qual = readFastq(query)
    for gene, allele in sorted(alleles.items()) :
        if allele['flag'] & 96 == 96 :
            alleles.pop(gene)
            continue
        if allele['flag'] & 32 > 0 :
            allele['sequence'] = 'DUPLICATED'
        else :
            c, s, e = allele['coordinates'][0]
            ss = seq[c][s - 1:e] if s < e else rc(seq[c][e - 1:s])
            qs = (min(qual[c][s - 1:e] if s < e else qual[c][e-1:s])) if len(qual) else 0
            
            allele['sequence'] = ss
            if qs < 10 :
                allele['flag'] |= 2
            if allele['flag'] == 0 :
                allele['flag'] = 1
        allele['coordinates'] = ','.join(['{0}:{1}..{2}'.format(*c) for c in  allele['coordinates']])
        allele['value_md5'] = ('' if allele['flag'] < 16 else '-') + get_md5(allele['sequence'])
    return {g:alleles.get(g, {"gene_name":g, "value_md5":"-"}) for g in core }


if __name__ == '__main__':
    main()
