#! /usr/bin/python3

import click, tempfile, gzip, os, sys

from configure import logging
from amr_search import amr_search
from gene_search import gene_search
from cgMLST import cgmlst
from PlasT import plast


def readFile(fname) :
    fin = gzip.open(fname, 'rt') if fname.upper().endswith('.GZ') else open(fname, 'rt')
    text = fin.readlines()
    fin.close()
    if text[0].startswith('@') :
        text = [ (t if i % 4 == 1 else '>' + t[1:]) for i, t in enumerate(text) if i % 4 < 2 ]
    return ''.join(text)


def write_cgMLST(prefix, profiles) :
    genes = sorted({gene for genome, prof in profiles for gene in prof.keys()})
    fout = '{0}.cgMLST.profile.gz'.format(prefix)
    with gzip.open(fout, 'at') as fout :
        fout.write('#Query\t{0}\n'.format('\t'.join(genes)))
        for genome, alleles in profiles :
            fout.write('{0}\t{1}\n'.format(genome, 
                                        '\t'.join([alleles.get(g, {}).get('value_md5', '-') for g in genes])))
    


@click.command()
@click.option('-q', '--query', help='query genome in fasta or fastq format. May be gzipped.')
@click.option('--ql', help='a list of query files. One query per line.')
@click.option('-o', '--prefix', help='prefix for output. Only work when there is only one query. default: query filename', default=None)
@click.option('-n', '--n_proc', help='number of process to use. default: 8', default=8, type=int)
@click.option('-m', '--skip_gene', help='flag to skip AMR/VF searching. default: False', default=False, is_flag=True)
@click.option('-g', '--skip_cgmlst', help='flag to skip cgMLST. default: False', default=False, is_flag=True)
@click.option('-p', '--skip_plasmid', help='flag to skip plasmid typing. default: False', default=False, is_flag=True)
def klety(query, ql, prefix, skip_gene, skip_cgmlst, skip_plasmid, n_proc) :
    queries = [] if not query else [query]
    if ql :
        with open(ql, 'rt') as fin :
            for line in fin :
                queries.append(line.strip().split()[0])
    
    if not prefix :
        if ql :
            prefix = os.path.basename(ql).rsplit('.', 1)[0]
        else :
            prefix = os.path.basename(query).rsplit('.', 1)[0]

    with open(prefix + '.KleTy', 'wt') as klety_out :
        resistances = 'AMR:AMINOGLYCOSIDE|AMR:BETA-LACTAM|AMR:CARBAPENEM|AMR:ESBL|AMR:INHIBITOR-RESISTANT|AMR:COLISTIN|AMR:FOSFOMYCIN|AMR:MACROLIDE|AMR:PHENICOL|AMR:QUINOLONE|AMR:RIFAMYCIN|AMR:GLYCOPEPTIDES|AMR:SULFONAMIDE|AMR:TETRACYCLINE|AMR:TIGECYCLINE|AMR:TRIMETHOPRIM|AMR:BLA_INTRINSIC|STRESS:COPPER|STRESS:MERCURY|STRESS:NICKEL|STRESS:SILVER|STRESS:TELLURIUM|STRESS:ARSENIC|STRESS:FLUORIDE|STRESS:QUATERNARY_AMMONIUM|VIRULENCE:clb|VIRULENCE:iro|VIRULENCE:iuc|VIRULENCE:rmp|VIRULENCE:ybt|Others|REPLICON:INC_TYPE|REPLICON:MOB_TYPE|REPLICON:MPF_TYPE'.split('|')
        fields = ['INPUT', 'REPLICON', 'SPECIES', 'HC1360.500.200.100.50.20.10.5.2', 'REFERENCE', 'PLASTYPE', 'COVERAGE'] + \
            resistances + ['ANNOTATION', 'CONTIGS']
        klety_out.write('\t'.join(fields)+'\n')

        if os.path.isfile(prefix + '.cgMLST.profile.gz') :
            os.unlink(prefix + '.cgMLST.profile.gz')

        profiles = []
        for query in queries :
            logging.info('Running query: {0}'.format(query))
            alleles = klebtyper(query, prefix, klety_out, skip_gene, skip_cgmlst, skip_plasmid, n_proc)
            profiles.append([query, alleles])
        write_cgMLST(prefix, profiles)


def klebtyper(query, prefix, klety_out, skip_gene, skip_cgmlst, skip_plasmid, n_proc) :
    # if prefix == None :
    #     prefix = os.path.basename(query).rsplit('.', 1)[0]
    
    with tempfile.TemporaryDirectory(dir='.', prefix='kt_') as tmpdir :
        tmpfile = os.path.join(tmpdir, 'qry.fna')
        with open(tmpfile, 'wt') as fout :
            fout.write(readFile(query))

        if not skip_gene :
            logging.info('\tSearching VF/STRESS genes...')
            genes = gene_search(tmpfile, n_proc=n_proc)
            logging.info('\tDone.')
            logging.info('\tSearching AMR genes...')
            amrs = amr_search(tmpfile, n_proc=n_proc)
            logging.info('\tDone.')
        else :
            genes, amrs = [], {}
        
        if not skip_plasmid :
            logging.info('\tSearching plasmids...')
            plasmids = plast(tmpfile, n_proc=n_proc)
            plasmids = [p for p in plasmids if p[0].startswith('PT')]
            logging.info('\tDone.')
        else :
            plasmids = []
        
        cont_replicon = {cont:i for i, plasmid in enumerate(plasmids) for cont in plasmid[3].split(',') }

        replicons = {ri:{} for ri in range(-1, len(plasmids)) }
        for gene in genes :
            replicon = cont_replicon.get(gene[0], len(plasmids))
            for category in set(gene[9].split('/') + gene[10].split('/')) :
                if category not in {'', 'MULTIDRUG', 'ORGANOMERCURY', 'PHENYLMERCURY', 'POLYKETIDE', 'ARSENITE', 'ARSENATE'} :
                    if category in {'ARSENIC', 'COPPER', 'SILVER', 'GOLD', 'FLUORIDE', 'MERCURY', 'TELLURIUM', 'NICKEL', 'QUATERNARY_AMMONIUM'} :
                        category = 'STRESS:{0}'.format(category)
                    elif category in ('INC_TYPE', 'MOB_TYPE', 'MPF_TYPE') :
                        category = 'REPLICON:{0}'.format(category)
                    elif category in ('iuc', 'rmp', 'ybt', 'clb', 'iro') :
                        category = 'VIRULENCE:{0}'.format(category)
                    else :
                        category = 'AMR:{0}'.format(category)
                    
                    if replicon not in replicons :
                        replicons[replicon] = {}
                    if category not in replicons[replicon] :
                        replicons[replicon][category] = {gene[7]}
                    else :
                        replicons[replicon][category].update({gene[7]})
                    if category not in replicons[-1] :
                        replicons[-1][category] = {gene[7]}
                    else :
                        replicons[-1][category].update({gene[7]})

        for category, genes in amrs.items() :
            for gene, contig in genes :
                replicon = cont_replicon.get(contig, len(plasmids))
                if replicon not in replicons :
                    replicons[replicon] = {}
                if category not in replicons[replicon] :
                    replicons[replicon][category] = {gene}
                else :
                    replicons[replicon][category].update({gene})
                if category not in replicons[-1] :
                    replicons[-1][category] = {gene}
                else :
                    replicons[-1][category].update({gene})

        if skip_plasmid :
            replicons.pop(0, None)

        resistances = 'AMR:AMINOGLYCOSIDE|AMR:BETA-LACTAM|AMR:CARBAPENEM|AMR:ESBL|AMR:INHIBITOR-RESISTANT|AMR:COLISTIN|AMR:FOSFOMYCIN|AMR:MACROLIDE|AMR:PHENICOL|AMR:QUINOLONE|AMR:RIFAMYCIN|AMR:GLYCOPEPTIDES|AMR:SULFONAMIDE|AMR:TETRACYCLINE|AMR:TIGECYCLINE|AMR:TRIMETHOPRIM|AMR:BLA_INTRINSIC|STRESS:COPPER|STRESS:MERCURY|STRESS:NICKEL|STRESS:SILVER|STRESS:TELLURIUM|STRESS:ARSENIC|STRESS:FLUORIDE|STRESS:QUATERNARY_AMMONIUM|VIRULENCE:clb|VIRULENCE:iro|VIRULENCE:iuc|VIRULENCE:rmp|VIRULENCE:ybt|Others|REPLICON:INC_TYPE|REPLICON:MOB_TYPE|REPLICON:MPF_TYPE'.split('|')
        for r, drugs in replicons.items() :
            for drug in resistances :
                replicons[r][drug] = ','.join(sorted(replicons[r][drug])) if drug in replicons[r] else '-'
            for drug in sorted(drugs.keys()) :
                if drug not in resistances :
                    x = '{0}({1})'.format(','.join(sorted(replicons[r][drug])), drug)
                    if replicons[r]['Others'] == '-' :
                        replicons[r]['Others'] = x
                    else :
                        replicons[r]['Others'] += '|' + x
        
        hiercc = {}
        alleles = {}
        if not skip_cgmlst :
            logging.info('\tRunning cgMLST...')
            alleles, hiercc = cgmlst(tmpfile, n_thread=n_proc)
            # write_cgMLST(query, alleles)
            logging.info('\tDone.')

        fields = ['INPUT', 'REPLICON', 'SPECIES', 'HC1360.500.200.100.50.20.10.5.2', 'REFERENCE', 'PLASTYPE', 'COVERAGE'] + \
            resistances + ['ANNOTATION', 'CONTIGS']
        for rep_id in sorted(replicons.keys()) :
            if rep_id == -1 :
                res = [query, 'ALL', hiercc.get('species', 'ND'), hiercc.get('HC1360.500.200.100.50.20.10.5.2', 'ND'), hiercc.get('reference','ND'), '-', '-'] + \
                    [replicons[-1][r] for r in resistances] + ['-', '-']
            elif rep_id >= len(plasmids) :
                res = [query, 'Others', '-', '-', '-', '-', '-'] + [replicons[rep_id][r] for r in resistances] + \
                    ['-', '-']
            else :
                plasmid = plasmids[rep_id]
                ref, ann = plasmid[-1].split(' ', 1)
                ref = ref.split('_', 2)[2]
                res = [query, 'P{0}'.format(rep_id+1), '-', '-', ref, plasmid[0].replace('\t', ','), plasmid[1]] + \
                    [replicons[rep_id][r] for r in resistances] + [ann.replace(' ', '_'), plasmid[3]]
            klety_out.write('\t'.join(res)+'\n')
    return alleles

if __name__ == '__main__' :
    klety()
