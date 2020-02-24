import os,sys
from os import path
from distutils import spawn

#
#   v1.0  - 15.07       created for regular DFG-out modeling
#   v2.0  - 16.06.07    added C-out parameters -- depreciated
#   v3.0  - 16.08.11    added primary model generation

def CheckEssentFiles(script_directory, database_dir, Settings):

  Vars = [  '1atp.pdb', 
            'stdy_kinase.clean_3.180815.tcoffee_1d.fasta',
            'stdy_kinase.clean_3.180815.tcoffee_1d.nogap.fasta',
            'u_cido.stkinase_templ_pdb.txt', 
            'u_cido.ykinase_templ_pdb.txt',
            'u_codi.templ_pdb.1.txt',
            'u_codi.templ_pdb.2.txt',
            'u_codi.templ_pdb.3.txt',
            'u_codi.templ_pdb.4.txt',
            'u_codi.templ_pdb.III.txt',
            'u_codi.templ_pdb.met.txt',
            'u_codi.templ_pdb.egfr.txt',
            'u_codo.templ_pdb.txt',
            'CommonUtility.py', 
            'x_variables.py',           'x_variables_run.py',
            'x_homolog_templ_check.py', 'x_build_multi_pdb.py',
            'x_pir_gen.py',  'x_pir_multi_align.py',   'x_pir_edit.py', 
            'x_modl_gen.py', 'x_modeller_parallel.py', 'x_mod_class.py', 
            'x_vol_gen.py',  'x_pymol_alignment.py',
          ]

  print( '\n -- Checking Components --')

  ## Check if script path is correct
  if not os.path.exists(script_directory):
    if not os.path.exists(database_dir):
      sys.exit('\n  Error: [ScriptDirectory] is not found')

  ## Check if required files are present
  for var in Vars:
    if os.path.isfile('{0}/{1}'.format(script_directory, var)):
      print('Found\t{0}'.format(var))
    elif os.path.isfile('{0}/{1}'.format(database_dir, var)):
      print( 'Found\t{0}'.format(var))
    else:
      sys.exit('\n  Error: {0} does not exist.'.format(var))

  ## Check if 3rd party software are present
  Soft = ['pymol', 'POVMELocation']
  if not spawn.find_executable(Settings[Soft[1]]):
    sys.exit('\n  Error: Software cannot be found: {0}'.format(Soft[1]))

  ## Check if Python modules are present
  Modules = ['modeller', 'Bio']

  for mod in Modules:
    try:
      __import__(mod)
    except ImportError:
      sys.exit('\n  Error: Python module cannot be found: {0}'.format(mod))
