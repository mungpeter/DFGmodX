#!/usr/bin/env python3

import sys,os,re
from shutil import copyfile
from Bio import SeqIO

msg = '''
  > {0}
\t[ kinase fasta ]\n\t[ model type: cidi|cido|codi|codo|wcd ]\n\t[ output folder ]
\t[ -c: Check completion | -p: Copy Over ]\n'''.format(sys.argv[0])
if len(sys.argv) != 5 : sys.exit(msg)

##########################################################################
def main( fasta_database, model, out_fold, option ):

  unip_list = '/home/pmung/Dropbox/9_scripts/3_program/structures/3_DFGmodx/x_dataset/kinase_unprot_id.txt'
  Uniprot = CollectUniProt(unip_list)
  Fasta   = [ fasta for fasta in SeqIO.parse(fasta_database, 'fasta') ]

  for fasta in Fasta:
    # Format of fasta name from KinBase: >TK:SRC/ABL1 xxxxxx xxxx xxxx
    family, name = fasta.id.split()[0].split('/')

    if name in Uniprot:
      unip = Uniprot[name]
    else:
      print('Different Name: '+name)
      continue

    if re.search(r'-p', option):
      CopyFiles( name, model, out_fold, unip )
    else:
      CheckFiles( name, model, out_fold, unip )

##########################################################################
## If the result folder has no *top_pdb.tar.bz2, i.e. the run failed, then:
## if no *.y1.corr.fasta, then the PIR file is not viable, check *y1.fasta
## if *.y1.corr.fasta presents, then PIR file didnt run correctly

def CheckFiles( name, model, out_fold, unip ):

  if not os.path.isfile('./{0}/1_result/{1}.{0}.top_pdb.tar.bz2'.format(name, model)):
    if os.path.isfile('./{0}/0_tempfile/_TEMP.{1}.{0}.y1.corr.fasta'.format(name, model)):
      print('pir didn\'t run: '+name)
    else:
      print('No viable pir: '+name)


##########################################################################

def CopyFiles( name, model, out_fold, unip ):
  if not os.path.exists('./{0}/{1}'.format(out_fold, name)):
    os.makedirs('./{0}/{1}'.format(out_fold, name))

  if not os.path.isfile('./{1}/{0}/{3}.{2}.{0}.top_pdb.tar.bz2'.format(name, out_fold, unip, model)):

    if os.path.isfile('./{0}/1_result/{1}.{0}.top_pdb.tar.bz2'.format(name, model)):
      copyfile('./{0}/1_result/{1}.{0}.top_pdb.tar.bz2'.format(name, model),
               './{1}/{0}/{3}.{2}.{0}.top_pdb.tar.bz2'.format(name, out_fold, unip, model))
    else:
      print('{0}   ./{0}/1_result/{1}.{0}.top_pdb.tar.bz2 not found'.format(name, model))
      
  if not os.path.isfile('./{1}/{0}/{3}.{2}.{0}.sample.pdb'.format(name, out_fold, unip, model)):
    if os.path.isfile('./{0}/1_result/{1}.{0}.sample.pdb'.format(name, model)):
      copyfile('./{0}/1_result/{1}.{0}.sample.pdb'.format(name, model),
               './{1}/{0}/{3}.{2}.{0}.sample.pdb'.format(name, out_fold, unip, model))
    

##########################################################################
def CollectUniProt( unip_list ):

  Uniprot = {}
  with open(unip_list, 'r') as fi:
    for l in fi:
      if re.search(r'^#', l):
        continue
      x = l.rstrip().split()
      if len(x) != 2:
        x.append('NaN')
      Uniprot[ x[0] ] = x[1]
  
  return Uniprot

  
##########################################################################
if __name__ == '__main__':
  main( sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4] ) 

##########################################################################
#
#   Peter M.U. Ung @ MSSM
#
#   v1.0    18.08.13
#
#   convert human protein kinase name to UniProt ID using existing data
#   also an option to append UniProt ID to filename and save the designated
#   files to certain folder
#
