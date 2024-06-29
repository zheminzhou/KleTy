import os, click, pandas as pd, numpy as np, re
from Bio import Seq
try :
    from uberBlast import uberBlast, readFastq
    from configure import dbs
except :
    from .uberBlast import uberBlast, readFastq
    from .configure import dbs


dbs = dbs['other_dbs']

def get_tops(bsn) :
    bsn = bsn[np.argsort(-bsn.T[11]*bsn.T[3]/bsn.T[12])]
    covs = { c:[] for c in np.unique(bsn.T[1]) }
    for p in bsn :
        s, e = p[8:10] if p[8] < p[9] else (p[9], p[8])
        x = np.zeros(e-s+1, dtype=np.uint8)
        for c in covs[p[1]] :
            if c[0] > e or c[1] < s :
                continue
            ovl = [max(c[0], s), min(c[1], e)]
            if ovl[1] - ovl[0] + 1 >= 0.3 * (c[1] - c[0] + 1) :
                p[1] = ''
                break
            x[ovl[0]-s:ovl[1]-s+1] = 1
            if np.sum(x) >= 0.3 * (e - s + 1) :
                p[1] = ''
                break
        else :
            covs[p[1]].append([s, e])
    return bsn[bsn.T[1] != '']
    

def rc(s) :
    c = {'A':'T', 'T':'A', 'G':'C', 'C':'G', '-':'-'}
    return ''.join([ c.get(x, 'N') for x in s[::-1].upper() ])


def gene_search(query, g_table=11, n_proc=8) :
    qry, _ = readFastq(query)
    res = []

    for title, (db, _, min_iden, min_cov) in dbs.items() : 
        bsn = uberBlast('-r {0} -q {1} -f --blastn --min_id {2} --min_ratio {3} -t {4} -p -s 2 -e 0,3 -m --merge_gap 300'.format(
            query, db, min_iden/100., min_cov/100., n_proc).split())            
        bsn = get_tops(bsn)

        for p in bsn:
            if title == 'nuc_stress' :
                ref = p[0].split(' ')[0].split('|')
                res.append([p[1], p[8], p[9], p[2]*100, np.round(100.*(p[7]-p[6]+1)/p[12], 2), ref[1], 'BLASTX', ref[4], ref[5], ref[8], ref[9], ref[10]])
            else :
                if title == 'nuc_vir' :
                    if p[0].find('delete') >= 0 :
                        continue
                    if p[0].find('__rmp') >= 0 :
                        p[0] = re.split('__', p[0])[2]
                    x = p[0].rsplit('_', 1)[0]
                    ref = ['', p[0], '', '', '', p[0], x, x, x]
                    category = p[0][:3]
                    if category in ('fyu', 'irp') :
                        category = 'ybt'
                    elif category in ('iut') :
                        category = 'iuc'
                elif title.endswith('plasmids') :
                    x = p[0].rsplit('|', 1)[-1]
                    ref = ['', p[0], '', '', '', x, x, x, x]
                    category = 'MOB_TYPE' if x.startswith('MOB') else ('MPF_TYPE' if x.startswith('MPF') else 'INC_TYPE')
                else :
                    x = p[0].rsplit('_', 1)[0]
                    ref = ['', p[0], '', '', '', p[0], x, x, x]
                    category = p[0][:3]

                qry_na = Seq.Seq(qry[p[1]][p[8]-1:p[9]]) if p[8] < p[9] else Seq.Seq(qry[p[1]][p[9]-1:p[8]]).reverse_complement()
                if title == 'nuc_plasmids' :
                    res.append([p[1], p[8], p[9], p[2]*100, np.round(100.*(p[7]-p[6]+1)/p[12], 2), ref[1], 'BLASTN', ref[5], ref[6], category, category, ref[7]])

                elif (abs(p[9] - p[8]) + 1) % 3 > 0 :
                    if (max([len(x) for x in qry_na[:int((len(qry_na))/3)*3].translate(g_table).split('*')])*3) >= 0.85 * len(qry_na) or \
                        (max([len(x) for x in qry_na[len(qry_na)%3:].translate(g_table).split('*')])*3) >= 0.85 * len(qry_na) :
                        res.append([p[1], p[8], p[9], p[2]*100, np.round(100.*(p[7]-p[6]+1)/p[12], 2), ref[1], 'BLASTN', ref[5], ref[6], category, category, ref[7]])
                    else :
                        res.append([p[1], p[8], p[9], p[2]*100, np.round(100.*(p[7]-p[6]+1)/p[12], 2), ref[1], 'BLASTN', ref[5] + '(*Frameshift)', ref[6], category, category, ref[7]])
                else :
                    qry_aa = qry_na.translate(g_table)
                    if '*' in qry_aa[:-1] and max([len(x) for x in qry_aa.split('*')])*3 < 0.85 * len(qry_na) :
                        res.append([p[1], p[8], p[9], p[2]*100, np.round(100.*(p[7]-p[6]+1)/p[12], 2), ref[1], 'BLASTN', ref[5] + '(*Premature)', ref[6], category, category, ref[7]])
                    else :
                        res.append([p[1], p[8], p[9], p[2]*100, np.round(100.*(p[7]-p[6]+1)/p[12], 2), ref[1], 'BLASTN', ref[5], ref[6], category, category, ref[7]])

    res = sorted([ r for r in res if r[9] != '' ])
    for j, r2 in list(enumerate(res))[1:] :
        todel = []
        for r1 in res[j-1::-1] :
            if r1[0] == '' :
                continue
            if r1[0] != r2[0] :
                break
            s1, e1 = sorted(r1[1:3])
            s2, e2 = sorted(r2[1:3])
            ovl = min(e2, e1) - max(s1, s2) + 1
            if ovl >= 0.5 * (e1-s1+1) or ovl >= 0.5 * (e2-s2+1) :
                sc1 = (e1-s1+1) * r1[3] * r1[4]/10000.
                sc2 = (e2-s2+1) * r2[3] * r2[4]/10000.
                if sc1 >= sc2 :
                    r2[0] = ''
                    break
                else :
                    todel.append(r1) #r1[0] = ''
        if r2[0] != '' :
            for r1 in todel :
                r1[0] = ''
    return [r for r in res if r[0] != '']


@click.command()
@click.option('-q', '--query', help='query file in FASTA format')
def main(query) :
    gene_search(query)



if __name__ == '__main__' :
    main()