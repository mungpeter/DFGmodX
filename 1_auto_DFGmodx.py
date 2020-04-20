#!/usr/bin/env python3

import sys

##########################################################################
#
#   Peter Man-Un Ung    @ MSSM
#
#   Purpose: automate running of multiple kinases (all in same conformation)
#            within a fasta database, instead of submitting them manually
#
##########################################################################

msg = '''\n\t> {0}
\t  [ FASTA file of kinase(s) to be modeled ]          * no gap in seq
\t  [ Conformation: cidi | cido | codi | codo | wcd ]  * capital or lowercase
\t  [ Number of models to generate ]
\t  [ Number of top-models to take ]
\t  [ CPU per Modeller run ]
\t  [ Step: 0 - Setup homology modeling (default)                   ]
\t  [       1 - Edit paths when change from local to HPC directory  ]
#\t  [       1 - Setup homology modeling (_TEMP.y.fasta corrected)   ]
#\t  [           #( rename *.y1.fasta to corrected *.y1.corr.fasta ) ]
\t  [       2 - Run homology modeling                               ]
\t  [       3 - Run volume calculation and selection                ]\n
\t  [ Optional Conformation: III\t- CODI MEK1-like for type-III inh ]
\t  [                        MET\t- CODI cMET-like altered A-loop   ]
\t  [                        EGFR\t-CODI EGFR-like altered A-loop    ]

e.g.> ~/Dropbox/9_scripts/3_program/structures/3_DFGmodx/1_auto_DFGmodx.py ~/Dropbox/9_scripts/3_program/structures/3_DFGmodx/z_dataset/MD_human_kinome_alignment.2019.200331.be_modeled.nogap.fasta cidi 50 10 10 0
'''.format(sys.argv[0])
if len(sys.argv) != 7: sys.exit(msg)

import re,os
import pickle

from Bio import SeqIO
from tqdm import tqdm
from pathos import multiprocessing

from x_variables import DefaultVariables
from x_variables_run import VariableSetup

##########################################################################
def main( fasta_file, conformation, 
          num_model, num_top, cpu_num, mod_step ):

  # current directory
  home_dir = str(os.getcwd())

  # work on each kinase sequence in supplied fasta file
  Data = list(SeqIO.parse(fasta_file, 'fasta'))[2:]
  print('\033[31m## Reading FASTA file of kinases to be modeled:\033[0m {0}'.format(len(Data)))

  ## Make .pir for each fasta sequence in their own directory
  if mod_step <= 1:
    if mod_step == 0:
      corr = False
    else:
      corr = True  # only for editing input when changing from local to HPC directory

    PIR  = RunPIRSetup( conf=conformation, cwd=home_dir, cpu=cpu_num, 
                        model=num_model, top=num_top, pref=None, corr=corr )

#    tmp = [ PIR(fasta) for fasta in tqdm(Data) ]
    mpi = multiprocessing.Pool()
    tmp = [x for x in tqdm(mpi.imap(PIR, Data), total=len(Data))]
    mpi.close()
    mpi.join()

  else:
  ## run homology modeling, skip if kinase directory is not already there
  ## OR run volume calculation and model sorting
    for idx, fasta in enumerate(Data):
      # Format of fasta: ><name>|<family>|<alt-name>/<resnum>
      name, family = fasta.id.split()[0].split('/')[0].split('|')[:2]

      if not os.path.isdir(name):
        print('  > #1# Warning: Directory not found. Skip: '+name)
        continue
      os.chdir(name)
      SetupVars = pickle.load( open(name+'.setup.pkl', 'rb') )

      if mod_step == 2:
        print('\033[34m## Running model generation:\033[0m Remaining \033[31m{0:3d} - \033[36m{1}\033[0m'.format(
                      len(Data)-idx, name))

#        ## for use with HPC only
#        os.system("sed 's/MODE/{0}/g' {1}/x_template.dfgmodx.lsf | sed 's/JBNAME/{2}/g' > {2}.dfgmodx.lsf; bsub < {2}.dfgmodx.lsf".format(
#                      'mod', SetupVars['ScriptDirectory'], name))

        command = '{0}/1_run_single_DFGmodx.py {1} -mod -pass > {2} '.format(
                      SetupVars['ScriptDirectory'], name+'.setup.pkl', name+'.mod_log')
        os.system(command)    # local machine only

        os.chdir(home_dir)    # return to home directory

      elif mod_step == 3:
        print('\033[34m## Running Volume Calculation and Model Sorting: \033[0mRemaining \033[31m{0:3d} - \033[36m{1}\033[0m'.format(
                      len(Data)-idx, name))

#        ## for use with HPC only
#        os.system("sed 's/MODE/{0}/g' {1}/x_template.dfgmodx.lsf | sed 's/JBNAME/{2}/g' > {2}.dfgmodx.lsf; bsub < {2}.dfgmodx.lsf".format(
#                      'vol', SetupVars['ScriptDirectory'], name))

        command = '{0}/1_run_single_DFGmodx.py {1} -vol -pass > {2} '.format(
                      SetupVars['ScriptDirectory'], name+'.setup.pkl', name+'.vol_log')
        os.system(command)    # local machine only

        os.chdir(home_dir)    # return to home directory)
      else:
        sys.exit('  > #2# ERROR: Wrong Step input (only accept 1-3): '+str(mod_step))


##########################################################################
# Generate PIR for each input fasta sequence
# 'prot' is biopython object of input fasta data
class RunPIRSetup( object ):
  def __init__( self, cwd='', conf='', cpu='', model='', 
                top='', pref='', corr='' ):
    self.cwd  = cwd
    self.conf = conf
    self.cpu  = cpu
    self.mod  = model
    self.top  = top
    self.pref = pref
    self.corr = corr

  def __call__( self, fasta ):
    return pir_setup( self, fasta )

###############
def pir_setup( self, fasta ):

  os.chdir(self.cwd)
  # Format of fasta: ><name>|<family>|<alt-name>/<resnum>
  name, family = fasta.id.split()[0].split('/')[0].split('|')[:2]

  # create and move into the working protein directory
  if not os.path.isdir(name):
    os.makedirs(name)
  os.chdir(name)

#### This is only used for editing the directory paths ####
  if self.corr:
    Vars = DefaultVariables()
    print('\033[34m## Editing Paths to Directories in Setup File: \033[31m{0}.setup\033[0m'.format(name))
    command = '{0}/1_run_single_DFGmodx.py {1} -none -paths '.format(
                  Vars['ScriptDirectory'], name+'.setup')
    print(command)
    os.system(command)
    print('  -- Finished Editing: '+name+'.setup --\n')
    return None

###########################################################


  ## Tell the script to use the sequence found in kinome database but
  ## specifying only the kinase name (in kinome database)
  if os.path.isfile(name+'.setup.pkl'):
    Vars = VariableSetup( pickle.load( open(name+'.setup.pkl', 'rb') ),
                            family, name, self.conf, self.cpu, 
                            self.mod, self.top, self.pref )
  else:
    Vars = VariableSetup( DefaultVariables(), 
                            family, name, self.conf, self.cpu, 
                            self.mod, self.top, self.pref )

  ## write out temporary setup files
  pickle.dump( Vars, open( name+'.setup.pkl', 'wb' ) )
  with open(name+'.setup', 'w') as fo:
    for key in Vars:
      fo.write('{0:20s} {1}\n'.format(key, Vars[key]))

  with open(name+'.fasta', 'w') as fo:
    SeqIO.write(fasta, fo, 'fasta')

  print('\033[34m## Running PIR generation: \033[31m{0}\033[0m'.format(name))
  command = '{0}/1_run_single_DFGmodx.py {1} -pir -restart > {2} '.format(
                  Vars['ScriptDirectory'], name+'.setup.pkl', name+'.pir_log')
  os.system(command)
  print('  -- Finished PIR generation: '+name+' --\n')


##########################################################################
if __name__ == '__main__':
  main( sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]),
        int(sys.argv[5]), int(sys.argv[6]) )

##########################################################################
#
# Peter M.U. Ung @ MSSM
#
# v1.0	18.03.07
#
# Run multiple DFGmod jobs
#
# Required packages:        use Mini/Anaconda 3
# python3         # 3.7.2+
#   numpy         # 1.17.3
#   pandas        # 0.25.3
#   biopython     # 1.74
#   zlib          # 1.2.11
#   bzip2         # 1.0.6
#   pathos        # 0.2.3
#   tqdm          # 
#   tarfile       #
#
# PyMOL           # 2.0.1+  standalone or conda
# blast+          # 2.6.0   standalone or conda
# clustalo        # 1.2+    standalone or conda
# muscle          # 3.8.15  standalone or conda
# t_coffee        # 11+     standalone or conda
# modeller        # 9.20+   standalone or conda
# POVME           # 2.1     source codes already integrated into DFGmodx
