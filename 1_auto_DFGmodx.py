#!/usr/bin/env python3

import re,sys,os
import pickle
import os.path
from Bio import SeqIO
from pathos import multiprocessing

from x_variables     import DefaultVariables
from x_variables_run import VariableSetup

##########################################################################
#
#   Peter Man-Un Ung    @ MSSM
#
#   Purpose: automate running of multiple kinases (all in same conformation)
#            within a fasta database, instead of submitting them manually
#
##########################################################################

msg = '''\n\t> {0}
\t  [ FASTA file of kinase(s) to be modeled ]
\t  [ Conformation: cidi | cido | codi | codo | wcd ]  *capital or lowercase
\t  [ Number of models to generate ]
\t  [ Number of top-models to take ]
\t  [ CPU per Modeller run ]
\t  [ Step: 0 - Setup homology modeling (default)                   ]
\t  [       1 - Setup homology modeling (_TEMP.y.fasta corrected)   ]
\t  [           #( rename *.y1.fasta to corrected *.y1.corr.fasta ) ]
\t  [       2 - Run homology modeling                               ]
\t  [       3 - Run volume calculation and selection                ]\n
\t  [ Optional Conformation: III\t- CODI MEK1-like for type-III inh ]
\t  [                        MET\t- CODI cMET-like altered A-loop   ]
\t  [                        EGFR\t-CODI EGFR-like altered A-loop    ]
'''.format(sys.argv[0])
if len(sys.argv) != 7: sys.exit(msg)


##########################################################################
def main( fasta_database, conformation, num_model, num_top, cpu_num, mod_step ):

  # current directory
  home_dir = str(os.getcwd())

  # work on each kinase sequence in supplied fasta file
  Data = [ fasta for fasta in SeqIO.parse(fasta_database, 'fasta') ]

  ## Make .pir for each fasta sequence in their own directory
  if mod_step <= 1:
    if mod_step == 0:
      corr = False
    else:
      corr = True

    PIR  = RunPIRSetup( conf=conformation, cwd=home_dir, cpu=cpu_num, 
                    model=num_model, top=num_top, corr=corr )

#    for fasta in Data: PIR(fasta)
    mpi = multiprocessing.Pool(processes = cpu_num)
    mpi.map( PIR, Data )
    mpi.close()
    mpi.join()

  else:
  ## run homology modeling, skip if kinase directory is not already there
  ## OR run volume calculation and model sorting
    for idx, fasta in enumerate(Data):
      # Format of fasta name from KinBase: >TK:SRC/ABL1 xxxxxx xxxx xxxx
      family, name = fasta.id.split()[0].split('/')
      print(name)
      if not os.path.isdir(name):
        print('  > #1# Warning: Directory not found. Skip: '+name)
        continue
      os.chdir(name)
      SetupVars = pickle.load( open('SetupVars.pkl', 'rb') )

      if mod_step == 2:
        print(' ## Running model generation: Remaining {0:3d} - {1}'.format(
                      len(Data)-idx, name))
        os.system('{0}/1_run_single_DFGmodx.py {1} -mod -pass > {2} '.format(
                      SetupVars['ScriptDirectory'],
                      'SetupVars.pkl', name+'.mod_log'))
        os.chdir(home_dir)    # return to home directory
        
      elif mod_step == 3:
        print(' ## Running Volume Calculation and Model Sorting: Remaining {0:3d} - {1}'.format(
                      len(Data)-idx, name))
        os.system('{0}/1_run_single_DFGmodx.py {1} -vol -pass > {2} '.format(
                      SetupVars['ScriptDirectory'],
                      'SetupVars.pkl', name+'.vol_log'))
        os.chdir(home_dir)    # return to home directory)
      else:
        sys.exit('  > #2# ERROR: Wrong Step input (only accept 1-3): '+str(mod_step))


##########################################################################
# Generate PIR for each input fasta sequence
# 'prot' is biopython object of input fasta data
class RunPIRSetup( object ):
  def __init__( self, cwd='', conf='', cpu='', model='', top='', corr=False ):
    self.cwd  = cwd
    self.conf = conf
    self.cpu  = cpu
    self.mod  = model
    self.top  = top
    self.corr = corr

  def __call__( self, fasta ):
    return pir_setup( self, fasta )

###############
def pir_setup( self, fasta ):

  os.chdir(self.cwd)
  # Format of fasta downloaded from 'KinBase.com': >TK:SRC/ABL1 xxxx xx xx
  family, name = fasta.id.split()[0].split('/')

  # create and move into the working protein directory
  if not os.path.isdir(name):
    os.makedirs(name)
  os.chdir(name)

  if os.path.isfile('SetupVars.pkl'):
    Vars = VariableSetup( pickle.load( open('SetupVars.pkl', 'rb') ),
                            family, name, self.conf, self.cpu, 
                            self.mod, self.top, self.corr, None )
  else:
    Vars = VariableSetup( DefaultVariables(), 
                            family, name, self.conf, self.cpu, 
                            self.mod, self.top, self.corr, None )
    
  pickle.dump( Vars, open( 'SetupVars.pkl', 'wb' ) )

  with open(name+'.setup', 'w') as fo:
    for key in Vars:
      fo.write('{0}\t\t{1}\n'.format(key, Vars[key]))

  with open(name+'.fasta', 'w') as fo:
    SeqIO.write(fasta, fo, 'fasta')

  print('## Running PIR generation: '+name+' ##')
  os.system('{0}/1_run_single_DFGmodx.py {1} -pir -pass > {2} '.format(
                  Vars['ScriptDirectory'], 
                  'SetupVars.pkl', name+'.pir_log'))
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
# Required packages:  * use Mini/Anaconda 3
#  python3       # 3.7.2
#  numpy         # 1.17.3
#  pandas        # 0.25.3
#  modeller      # 9.20+, conda or standalone
#  biopython     # 1.74
#  zlib          # 1.2.11
#  bzip2         # 1.0.6
#  pathos        # 0.2.3
#  PyMOL         # 2.0.1+, standalone
#  POVME         # 2 or 3, standalone
