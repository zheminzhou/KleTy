import os, sys, numpy as np, shutil
from subprocess import Popen, PIPE
import logging
logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)


dirname = os.path.dirname(os.path.abspath(__file__))

exe = dict(makeblastdb=shutil.which('makeblastdb'), 
           blastn=shutil.which('blastn'), 
           diamond=shutil.which('diamond'))


db_folder = os.path.join(dirname, 'db')

dbs = dict(refseqs = os.path.join(db_folder, 'klebsiella.refsets.fas'), 
          core_genes = os.path.join(db_folder, 'klebsiella.cgmlst'), 
          repr = os.path.join(db_folder, 'profile.parq'), 
          hcc = os.path.join(db_folder, 'HierCC.tsv.gz'), 
          species = os.path.join(db_folder, 'klebsiella.species'), 
          plast_clu = os.path.join(db_folder, 'plasmids.repr.clu'), 
          plast_db = os.path.join(db_folder, 'plasmids.repr.fas'), 
          amr_dbs = dict(
              amr_class = os.path.join(db_folder, 'kleborate', 'CARD_AMR_clustered.csv'), 
              CARD = os.path.join(db_folder, 'kleborate', 'CARD_v3.1.13.fasta'), 
              qrdr = os.path.join(db_folder, 'kleborate', 'QRDR_120.fasta'), 
              omp = os.path.join(db_folder, 'kleborate', 'OmpK.fasta'), 
              trunc = os.path.join(db_folder, 'kleborate', 'MgrB_and_PmrB.fasta'), 
            ),
          other_dbs = dict(
                nuc_stress=[os.path.join(db_folder, 'stress_CDS'),  None, 90, 60], 
                nuc_vir=[os.path.join(db_folder, 'kleborate', 'virulence.fasta'),  None, 90, 60], 
                nuc_plasmids=[os.path.join(db_folder, 'traditional_plasmid_type.fas'), None, 90, 80],
            )
)

scheme_info = dict(min_iden = 0.65, min_frag=0.6, max_iden=0.95)



blosum62 = np.array(
    [4., -2., 0., -2., -1., -2., 0., -2., -1., 0., -1., -1., -1., -2., 0., -1., -1., -1., 1., 0., -4., 0.,
     -3., 0., -2., -1., 0., 0., 0., 0., 0., 0., -2., 4., -3., 4., 1., -3., -1., 0., -3., 0., 0., -4.,
     -3., 3., 0., -2., 0., -1., 0., -1., -4., -3., -4., -1., -3., 1., 0., 0., 0., 0., 0., 0., 0., -3.,
     9., -3., -4., -2., -3., -3., -1., 0., -3., -1., -1., -3., 0., -3., -3., -3., -1., -1., -4., -1., -2., -2.,
     -2., -3., 0., 0., 0., 0., 0., 0., -2., 4., -3., 6., 2., -3., -1., -1., -3., 0., -1., -4., -3., 1.,
     0., -1., 0., -2., 0., -1., -4., -3., -4., -1., -3., 1., 0., 0., 0., 0., 0., 0., -1., 1., -4., 2.,
     5., -3., -2., 0., -3., 0., 1., -3., -2., 0., 0., -1., 2., 0., 0., -1., -4., -2., -3., -1., -2., 4.,
     0., 0., 0., 0., 0., 0., -2., -3., -2., -3., -3., 6., -3., -1., 0., 0., -3., 0., 0., -3., 0., -4.,
     -3., -3., -2., -2., -4., -1., 1., -1., 3., -3., 0., 0., 0., 0., 0., 0., 0., -1., -3., -1., -2., -3.,
     6., -2., -4., 0., -2., -4., -3., 0., 0., -2., -2., -2., 0., -2., -4., -3., -2., -1., -3., -2., 0., 0.,
     0., 0., 0., 0., -2., 0., -3., -1., 0., -1., -2., 8., -3., 0., -1., -3., -2., 1., 0., -2., 0., 0.,
     -1., -2., -4., -3., -2., -1., 2., 0., 0., 0., 0., 0., 0., 0., -1., -3., -1., -3., -3., 0., -4., -3.,
     4., 0., -3., 2., 1., -3., 0., -3., -3., -3., -2., -1., -4., 3., -3., -1., -1., -3., 0., 0., 0., 0.,
     0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
     0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., 0., -3., -1., 1., -3., -2., -1., -3., 0.,
     5., -2., -1., 0., 0., -1., 1., 2., 0., -1., -4., -2., -3., -1., -2., 1., 0., 0., 0., 0., 0., 0.,
     -1., -4., -1., -4., -3., 0., -4., -3., 2., 0., -2., 4., 2., -3., 0., -3., -2., -2., -2., -1., -4., 1.,
     -2., -1., -1., -3., 0., 0., 0., 0., 0., 0., -1., -3., -1., -3., -2., 0., -3., -2., 1., 0., -1., 2.,
     5., -2., 0., -2., 0., -1., -1., -1., -4., 1., -1., -1., -1., -1., 0., 0., 0., 0., 0., 0., -2., 3.,
     -3., 1., 0., -3., 0., 1., -3., 0., 0., -3., -2., 6., 0., -2., 0., 0., 1., 0., -4., -3., -4., -1.,
     -2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
     0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -1., -2., -3., -1.,
     -1., -4., -2., -2., -3., 0., -1., -3., -2., -2., 0., 7., -1., -2., -1., -1., -4., -2., -4., -2., -3., -1.,
     0., 0., 0., 0., 0., 0., -1., 0., -3., 0., 2., -3., -2., 0., -3., 0., 1., -2., 0., 0., 0., -1.,
     5., 1., 0., -1., -4., -2., -2., -1., -1., 3., 0., 0., 0., 0., 0., 0., -1., -1., -3., -2., 0., -3.,
     -2., 0., -3., 0., 2., -2., -1., 0., 0., -2., 1., 5., -1., -1., -4., -3., -3., -1., -2., 0., 0., 0.,
     0., 0., 0., 0., 1., 0., -1., 0., 0., -2., 0., -1., -2., 0., 0., -2., -1., 1., 0., -1., 0., -1.,
     4., 1., -4., -2., -3., 0., -2., 0., 0., 0., 0., 0., 0., 0., 0., -1., -1., -1., -1., -2., -2., -2.,
     -1., 0., -1., -1., -1., 0., 0., -1., -1., -1., 1., 5., -4., 0., -2., 0., -2., -1., 0., 0., 0., 0.,
     0., 0., -4., -4., -4., -4., -4., -4., -4., -4., -4., 0., -4., -4., -4., -4., 0., -4., -4., -4., -4., -4.,
     1., -4., -4., -4., -4., -4., 0., 0., 0., 0., 0., 0., 0., -3., -1., -3., -2., -1., -3., -3., 3., 0.,
     -2., 1., 1., -3., 0., -2., -2., -3., -2., 0., -4., 4., -3., -1., -1., -2., 0., 0., 0., 0., 0., 0.,
     -3., -4., -2., -4., -3., 1., -2., -2., -3., 0., -3., -2., -1., -4., 0., -4., -2., -3., -3., -2., -4., -3.,
     11., -2., 2., -3., 0., 0., 0., 0., 0., 0., 0., -1., -2., -1., -1., -1., -1., -1., -1., 0., -1., -1.,
     -1., -1., 0., -2., -1., -1., 0., 0., -4., -1., -2., -1., -1., -1., 0., 0., 0., 0., 0., 0., -2., -3.,
     -2., -3., -2., 3., -3., 2., -1., 0., -2., -1., -1., -2., 0., -3., -1., -2., -2., -2., -4., -1., 2., -1.,
     7., -2., 0., 0., 0., 0., 0., 0., -1., 1., -3., 1., 4., -3., -2., 0., -3., 0., 1., -3., -1., 0.,
     0., -1., 3., 0., 0., -1., -4., -2., -3., -1., -2., 4., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
     0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.])


def check_executables() :
  noFound = []
  if not os.path.isfile(exe['blastn']) :
    noFound.append('blastn')
  if not os.path.isfile(exe['diamond']) :
    noFound.append('diamond')
  if noFound :
    logging.error('{0} is not found.'.format(' and '.format(noFound)))
    sys.exit(-1)


def check_db_files() :
  logging.info('Checking database files.')

  if not os.path.isfile(dbs['refseqs']) :
    if os.path.isfile(dbs['refseqs']+'.gz') :
      Popen('gzip -d {0}'.format(dbs['refseqs']+'.gz'), shell=True).communicate()
    else :
      logging.error('{0} is not found.'.format(dbs['refseqs']))
      sys.exit(-1)

  if not os.path.isfile(dbs['other_dbs']['nuc_vir'][0]) :
    if os.path.isfile(dbs['other_dbs']['nuc_vir'][0]+'.gz') :
      Popen('gzip -d {0}'.format(dbs['other_dbs']['nuc_vir'][0]+'.gz'), shell=True).communicate()
    else :
      logging.error('{0} is not found.'.format(dbs['other_dbs']['nuc_vir'][0]))
      sys.exit(-1)

  if not os.path.isfile(dbs['plast_clu']) :
    if os.path.isfile(dbs['plast_clu']+'.gz') :
      Popen('gzip -d {0}'.format(dbs['plast_clu']+'.gz'), shell=True).communicate()
    else :
      logging.error('{0} is not found.'.format(dbs['plast_clu']))
      sys.exit(-1)

  if not os.path.isfile(dbs['other_dbs']['nuc_stress'][0]) :
    tmp_gz = dbs['other_dbs']['nuc_stress'][0]+'.gz'
    if os.path.isfile(tmp_gz) :
      Popen('gzip -d {0}'.format(tmp_gz), shell=True).communicate()
      Popen('{0} -in {1} -dbtype nucl'.format(exe['makeblastdb'], dbs['other_dbs']['nuc_stress'][0]), shell=True).communicate()
    else :
      logging.error('{0} is not found.'.format(dbs['other_dbs']['nuc_stress'][0]))
      sys.exit(-1)

  if not os.path.isfile(dbs['other_dbs']['nuc_plasmids'][0]) :
    tmp_gz = dbs['other_dbs']['nuc_plasmids'][0]+'.gz'
    if os.path.isfile(tmp_gz) :
      Popen('gzip -d {0}'.format(tmp_gz), shell=True).communicate()
      Popen('{0} -in {1} -dbtype nucl'.format(exe['makeblastdb'], dbs['other_dbs']['nuc_plasmids'][0]), shell=True).communicate()
    else :
      logging.error('{0} is not found.'.format(dbs['other_dbs']['nuc_plasmids'][0]))
      sys.exit(-1)

  if not os.path.isfile(dbs['plast_db']+'.nhr') :
    logging.info('Plasmid reference file is not found. Will download it from https://zenodo.org/records/12590507/files/plasmids.repr.fas.gz.')
    logging.info('This may take a while but will run only once.')
    tmp_gz = dbs['plast_db']+'.gz'

    Popen('wget -O {0} https://zenodo.org/records/12590507/files/plasmids.repr.fas.gz'.format(tmp_gz), shell=True).communicate()
    Popen('gzip -cd {0} > {1}'.format(tmp_gz, dbs['plast_db']), shell=True).communicate()
    Popen('{0} -in {1} -dbtype nucl'.format(exe['makeblastdb'], dbs['plast_db']), shell=True, stdout=PIPE, stderr=PIPE).communicate()
    os.unlink(dbs['plast_db'])


check_executables()
check_db_files()