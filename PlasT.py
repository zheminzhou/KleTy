#! /usr/bin/python3
import numpy as np, os, subprocess, tempfile, re
import click, pandas as pd, gzip
try :
    from .configure import exe, dbs
except :
    from configure import exe, dbs



group = dbs['plast_clu']
db = dbs['plast_db']
blastn = exe['blastn']



@click.command()
@click.option('-i', '--ani', help='minimum ANI for clustering. default: 85', type=float, default=85)
@click.option('-c', '--afr', help='minimum proportion of reference in alignments. default: 50', type=float, default=50.)
@click.option('-C', '--afq', help='minimum proportion of query contig in alignments. default: 30', type=float, default=30.)
@click.option('-f', '--frag_ref', help='minimum proportion of reference in alignments as a fragment of plasmid. default: 30', type=float, default=30.)
@click.option('-F', '--frag_qry', help='minimum proportion of query contig in alignments as a fragment of plasmid. default: 50', type=float, default=50.)
@click.option('-b' ,'--bsn', help='blastn results', default=None)
@click.argument('query')
def main(query, bsn, ani, afr, afq, frag_ref, frag_qry) :
    res = plast(query, bsn, ani, afr, afq, frag_ref, frag_qry)
    for cls, cover, iden, ctg, info in res :
        print('{0}\tcover={1}\tiden={2}\t{3}\t{4}'.format(cls, cover, iden, ctg, info))


def plast(query, bsn=None, ani=85., afr=50., afq=30., frag_ref=30., frag_qry=50., n_proc=8) :
    with tempfile.TemporaryDirectory(dir='.') as dname :
        groups, cluster_info, weights, partial_ref = parseGroup(group)
        if not bsn :
            if query.lower().endswith('.gz') :
                if query.lower().find('fastq') > 0 :
                    with open(os.path.join(dname, 'query'), 'wt') as fout, gzip.open(query, 'rt') as fin :
                        for i, line in enumerate(fin) :
                            if i % 4 == 0 :
                                fout.write('>' + line[1:])
                            elif i % 4 == 1 :
                                fout.write(line)

                else :
                    with open(os.path.join(dname, 'query'), 'wt') as fout, gzip.open(query, 'rt') as fin :
                        fout.write(fin.read())
                query = os.path.join(dname, 'query')
            else :
                if query.lower().find('fastq') > 0 :
                    with open(os.path.join(dname, 'query'), 'wt') as fout, open(query, 'rt') as fin :
                        for i, line in enumerate(fin) :
                            if i % 4 == 0 :
                                fout.write('>' + line[1:])
                            elif i % 4 == 1 :
                                fout.write(line)
                    query = os.path.join(dname, 'query')

            bsn = os.path.join(dname, 'bsn')
            subprocess.Popen('{3} -query {0} -num_threads {4} -max_target_seqs 5000 -out {2} -db {1} -outfmt'.format(
                query, db, bsn, blastn, n_proc).split() +
                    ['6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'],
                                ).communicate()
        res = parseBsn(bsn, groups, cluster_info, weights, ani, afr, afq, frag_ref, frag_qry, partial_ref)
        return res


def parseGroup(group) :
    groups = {}
    cluster_info = {}
    weights = {}
    partial_ref = {}
    with open(group, 'rt') as fin :
        for line in fin :
            p = line.strip().split('\t')
            if not p[0].startswith('#') :
                grp = '\t'.join([p[1], p[2]])
                groups[p[0]] = grp
                partial_ref[p[0]] = 1 if len(p) > 5 and p[5] == 'partial' else 0
                if p[1].startswith('P') and p[4].startswith('V') and grp in cluster_info : 
                    pass
                else :
                    cluster_info[grp] = p[4]
                if grp not in weights :
                    weights[grp] = 1
                else :
                    weights[grp] += 1
    return groups, cluster_info, {g:(np.log(w)/np.log(2)+1.) for g, w in weights.items() }, partial_ref


def eachBsn(bsn, contigs, af_r, af_q) :
    if len(bsn) == 0 :
        return 0
    bsn = bsn[np.argsort(bsn.T[6])]
    bsn = bsn[np.argsort(bsn.T[0], kind='mergesort')]

    cont_info = {bsn[0][0]:[bsn[0][12], bsn[0][7]-bsn[0][6]+1]}
    region = [[bsn[0][0], bsn[0][6], bsn[0][7]]]
    for p in bsn[1:] :
        if p[0] != region[-1][0] or p[6] > region[-1][2]+101 :
            region.append([p[0], p[6], p[7]])
            if p[0] not in cont_info :
                cont_info[p[0]] = [p[12], p[7]-p[6]+1]
            else :
                cont_info[p[0]][1] += p[7] - p[6] + 1
        elif p[7] > region[-1][2] :
            cont_info[p[0]][1] += p[7] - region[-1][2]
            region[-1][2] = p[7]
    cont_info = {c:r[1]/r[0] for c, r in cont_info.items() if r[0]*af_q/100. <= r[1]}
    bsn = bsn[ [c in cont_info for c in bsn.T[0]] ]
    if len(bsn) == 0 :
        return 0

    bsn = bsn[np.argsort(np.min(bsn[:, 8:10], 1))]
    region = [[bsn[0][8], bsn[0][9]]] if bsn[0][8] < bsn[0][9] else [[bsn[0][9], bsn[0][8]]]
    for p in bsn[1:] :
        s, e = [p[8], p[9]] if p[8] < p[9] else [p[9], p[8]]
        if s > region[-1][1]+101 :
            region.append([s, e])
        elif region[-1][1] < e :
            region[-1][1] = e
    tot_region = sum([ r[1] - r[0] + 1 for r in region ])
    if tot_region < af_r/100. * bsn[0][-1] :
        return tot_region/bsn[0][-1]

    for p in bsn :
        if p[0] not in contigs :
            contigs[p[0]] = [p]
        else :
            contigs[p[0]].append(p)
    return tot_region/bsn[0][-1]


def parseBsn(bsn, groups, cluster_info, weights, ani, af_r, af_q, frag_ref, frag_qry, partial_ref) :
    data = pd.read_csv(bsn, sep='\t', header=None)
    data = data.loc[data[2] >= ani]
    data = data.loc[data[3] >= 400]
    data = data.loc[[r in groups for r in data[1]]]
    data = data.sort_values(by=1)
    contLen = { c:l for c, l in data[[0, 12]].values }
    data = data.values

    data = data[[groups.get(r, '') != '' for r in data.T[1]]]
    res = []
    while data.shape[0] :
        contigs = {}
        ref_cover = {}
        starts = [0] + (np.where(data[:-1, 1] != data[1:, 1])[0] + 1).tolist() + [len(data)]
        for s, e in zip(starts[:-1], starts[1:]) :
            r = eachBsn(data[s:e], contigs, af_r, af_q)
            c = groups[data[s, 1]]
            if r > ref_cover.get(c, 0) :
                ref_cover[c] = r

        cont_info = {}
        for cont, regions in contigs.items() :
            cont_seqs = { c:np.zeros(contLen[cont]) for c in set([groups[r[1]] for r in regions]) }
            for r in regions :
                p = r[2]
                s = cont_seqs[groups[r[1]]]
                s[r[6]-1:r[7]][s[r[6]-1:r[7]] < p] = p
            cont_info[cont] = {}
            for c, s in cont_seqs.items() :
                cont_info[cont][c] = [np.sum(s), np.sum((s>0)), np.sum((s>0))/s.size*100.]
        clusters = {}
        for cont, info in cont_info.items() :
            for cls, mat in info.items() :
                if mat[2] >= af_q :
                    if cls not in clusters :
                        clusters[cls] = [mat[0], [[mat[0], mat[1], cont]]]
                    else :
                        clusters[cls][0] += mat[0]
                        clusters[cls][1].append([mat[0], mat[1], cont])
        if len(clusters) :
            for cls, (score, mats) in clusters.items() :
                acc = [0, 0]
                for m in mats :
                    acc[0] += m[0]
                    acc[1] += m[1]
                    if acc[0] > score * 0.5 :
                        break
                iden = acc[0]/acc[1]
                mats = [ m for m in mats if m[1] >= iden - 0.1 ]
                clusters[cls][0] = sum([m[0] for m in mats])*weights[cls]*(ref_cover[cls]-af_r/100.+0.05)*(ref_cover[cls]-af_r/100.+0.05)
                clusters[cls][1] = sum([m[0] for m in mats])/sum([m[1] for m in mats])
                clusters[cls].append([m[2] for m in mats])

            cls, (score, iden, conts) = max(clusters.items(), key=lambda c:[c[1][0], c[0]])
            data = data[[groups[r] != cls for r in data.T[1]]]
            data = data[[ c not in conts for c in data.T[0] ]]
            res.append([cls, '{0:.1f}'.format(ref_cover[cls]*100.), '{0:.1f}'.format(iden), ','.join(sorted(conts)), cluster_info[cls]])
        else :
            break

    while data.shape[0] :
        contigs = {}
        ref_cover = {}
        starts = [0] + (np.where(data[:-1, 1] != data[1:, 1])[0] + 1).tolist() + [len(data)]
        for s, e in zip(starts[:-1], starts[1:]) :
            ref = data[s][1]
            if partial_ref.get(ref, 0) == 1 :
                continue
            r = eachBsn(data[s:e], contigs, frag_ref, frag_qry)
            c = groups[data[s, 1]]
            if r > ref_cover.get(c, 0) :
                ref_cover[c] = r


        cont_info = {}
        for cont, regions in contigs.items() :
            cont_seqs = { c:np.zeros(contLen[cont]) for c in set([groups[r[1]] for r in regions]) }
            for r in regions :
                p = r[2]
                s = cont_seqs[groups[r[1]]]
                s[r[6]-1:r[7]][s[r[6]-1:r[7]] < p] = p
            cont_info[cont] = {}
            for c, s in cont_seqs.items() :
                cont_info[cont][c] = [np.sum(s), np.sum((s>0)), np.sum((s>0))/s.size*100.]
        clusters = {}
        for cont, info in cont_info.items() :
            for cls, mat in info.items() :
                if mat[2] >= frag_qry :
                    if cls not in clusters :
                        clusters[cls] = [mat[0], [[mat[0], mat[1], cont]]]
                    else :
                        clusters[cls][0] += mat[0]
                        clusters[cls][1].append([mat[0], mat[1], cont])
        if len(clusters) :
            for cls, (score, mats) in clusters.items() :
                acc = [0, 0]
                for m in mats :
                    acc[0] += m[0]
                    acc[1] += m[1]
                    if acc[0] > score * 0.5 :
                        break
                iden = acc[0]/acc[1]
                mats = [ m for m in mats if m[1] >= iden - 0.1 ]
                clusters[cls][0] = sum([m[0] for m in mats])*weights[cls]*(ref_cover[cls]-frag_ref/100.+0.05)*(ref_cover[cls]-frag_ref/100.+0.05)
                clusters[cls][1] = sum([m[0] for m in mats])/sum([m[1] for m in mats])
                clusters[cls].append([m[2] for m in mats])

            cls, (score, iden, conts) = max(clusters.items(), key=lambda c:[c[1][0], c[0]])
            data = data[[groups[r] != cls for r in data.T[1]]]
            data = data[[ c not in conts for c in data.T[0] ]]
            res.append([cls+' (Fragment)', '{0:.1f}*'.format(ref_cover[cls]*100.), '{0:.1f}'.format(iden), ','.join(sorted(conts)), cluster_info[cls]])
        else :
            break

    return res



if __name__ == '__main__' :
    main()
