import sys,os
import re,glob
import shutil
import tarfile
import subprocess
from x_pymol_alignment import PyMOLSuperpose

#####################################################################  
# Run modeller with the modified .pir file  
def GenerateAndAlignModeller(
        script_directory, home_directory, work_directory, result_directory,
        chimera_tmpl_list, mdl_pir_file, number_of_model, 
        number_of_cpu, mdl_output_pref,
        reference_pdb, best_match_struc, superpose_resi,
        pymol_exec, Settings ):

  Vars = [ 'ScriptDirectory', 'ChimeraTemplList', 'ModifiedPIRFile' ]

  print('\n  \033[31m-- DFGmodel Structure Generation --\033[0m')
  for var in Vars:
    if var == 'ChimeraTemplList':
      if os.path.exists(var):
        sys.exit('\n  Error: {0} does not exists.'.format(var))
    if Settings[var] is None:
      sys.exit('\n  Error: \'{0}\' is not specified: {0}'.format(var))

  #### Hard-coded variable
  aligned_mdl_list = mdl_output_pref+'.align_pdb.list'

#################################################

  # Generate DFGmodel
  RunModelGeneration(
        script_directory, work_directory, result_directory,
        chimera_tmpl_list, mdl_pir_file, number_of_model, 
        number_of_cpu, mdl_output_pref )

  # Align all generated models to a reference frame
  RunStructureAlignment(
        script_directory, work_directory, result_directory,
        Settings['DatasetDirectory']+'/'+reference_pdb, best_match_struc, 
        superpose_resi, mdl_output_pref, aligned_mdl_list, number_of_model, 
        pymol_exec )

  os.chdir(home_directory)


##########################################################################
## Run Modeller to generate the models
def RunModelGeneration(
        script_directory, work_directory, result_directory,
        chimera_tmpl_list, mdl_pir_file, number_of_model, 
        number_of_cpu, mdl_output_pref ):

  os.chdir(work_directory)
  mod = 'python3'

  print('\033[34m## Current directory:\033[0m\n{0}\n'.format(os.getcwd()))
  print( [chimera_tmpl_list, mdl_pir_file, number_of_model, 
          number_of_cpu, mdl_output_pref] )

  Log = []
  try:
    print(' >> Modeller: \033[34m{0}\033[0m'.format(mod))
    run_m = mod+' "{0}/x_modeller_parallel.py" "{1}/{2}" "{1}/{3}" {4} {5} {6}'.format(
          script_directory, work_directory, 
          chimera_tmpl_list, mdl_pir_file,
          number_of_model, number_of_cpu, mdl_output_pref )
    print('\n\033[31m{0}\033[0m\n'.format(run_m))
#    os.system(run_m)
    Log = subprocess.check_output(run_m, shell=True, universal_newlines=True)

    with open(mdl_output_pref+'.modeller.log', 'w') as fo:
      for line in Log: fo.write(line)
  except Exception as e:
    print(e.__doc__)
#    print(e.message)
    print('  > #1# WARNING: Modelling failed for: '+mdl_output_pref)
    return None

  # Print Z-DOPE score in the Modeller log file  
  ExtractZDOPEScore(result_directory, mdl_output_pref)

  os.system('rm -r {0}.D000* {0}.V999* *.rsr *.sch *.ini x_modeller_parallel.slave*'.format(mdl_output_pref))

  Model_List = sorted(glob.glob('{0}.B999*.pdb'.format(mdl_output_pref)))

  # Tar orig_pdb results to "result"; Replace existing tar.bz2 if it is in "result_directory"
  if os.path.isfile('{0}/{1}.orig_pdb.tar.bz2'.format( result_directory, mdl_output_pref)):
    os.remove('{0}/{1}.orig_pdb.tar.bz2'.format( result_directory, mdl_output_pref))
  tar = tarfile.open('{0}/{1}.orig_pdb.tar.bz2'.format( 
                      result_directory, mdl_output_pref), 'w:bz2')

  # writing into the tar_file in "result_directory"  
  with open(mdl_output_pref+'.orig_pdb.list', 'w') as fo:
    for model in Model_List:
      fo.write(model.strip()+'\n')
      tar.add(model)
  tar.add(mdl_output_pref+'.orig_pdb.list')
  tar.close()


##########################################################################
## Extract Z-DOPE scores from the Modeller output and sort them
def ExtractZDOPEScore( result_directory, mdl_output_pref ):
  Record = []
  Title  = []
  record = False
  with open(mdl_output_pref+'.modeller.log','r') as fi:
    for line in filter(None, (l.rstrip() for l in fi)):  # skip blank lines
      if record:
        if re.match(r'---', line) or re.match(r'Filename', line):
          Title.append(line+'\n')
          continue
        Items = line.split()
        Record.append(Items)
      if re.search(r'>> Summary of successfully produced models:', line):
        record = True

  Sorted = sorted(Record, key=lambda tup: (tup[5]), reverse=True)

  with open(mdl_output_pref+'.zDOPE.txt', 'w') as fo:
    for line in Title: fo.write(line)
    for Items in Sorted:
      if float(Items[5]) > 0.0:
        print('  > #4# BAD ALIGNMENT WARNING: {0} Model {1} zDOPE: {2:.2f}###'.format(
                  mdl_output_pref, Items[0], float(Items[5])))

      fo.write('{0:28s}\t{1:7.2f}\t{2:7.2f}\t{3:7.2f}\t{4:7.5f}\t{5:10.5f}\n'.format(
              Items[0], float(Items[1]), float(Items[2]),
              float(Items[3]), float(Items[4]), float(Items[5])))

  if os.path.isfile(mdl_output_pref+'.modeller.log.bz2'):
    os.remove(mdl_output_pref+'.modeller.log.bz2')

  os.system('bzip2 {0}'.format(mdl_output_pref+'.modeller.log'))
  os.system('cp *modeller.log.bz2 *.zDOPE.txt "{0}"'.format(result_directory))


##########################################################################
## Align modeller-generated structures to reference kinase (1ATP or 2BDF)
def RunStructureAlignment(
        script_directory, work_directory, result_directory, 
        reference_pdb, best_match_struc, superpose_resi, 
        mdl_output_pref, aligned_mdl_list, number_of_model,
        pymol_exec ):

  os.chdir(work_directory)
  print('\n  \033[31m** Running Structure Superposition **\033[0m\n')
  print('\033[34m## Current directory:\033[0m\n{0}\n'.format(os.getcwd()))

# Align DFGmodel structure to a reference structure
  PyMOLSuperpose( pymol_exec, reference_pdb, best_match_struc, superpose_resi, 
                  mdl_output_pref+'.orig_pdb.list', mdl_output_pref+'.align', 
                  aligned_mdl_list, number_of_model, 'mod.pdb' )

  # Make a copy of sample kinase model for checking and viewing
  Models = sorted(glob.glob('{0}*.mod.pdb'.format(mdl_output_pref)))

  print('\n ** {0} will be added to Result as {1}.sample.pdb'.format(
                    Models[0], Models[0].split('.mod')[0] ))
  os.system('cp {0}/{1} {2}/{3}.sample.pdb'.format(
                work_directory, Models[0], 
                result_directory, Models[0].split('.mod')[0] ))

  # Replace existing tar.bz2
  if os.path.isfile('{0}/{1}.mdl_pdb.tar.bz2'.format( result_directory, mdl_output_pref)):
    os.remove('{0}/{1}.mdl_pdb.tar.bz2'.format( result_directory, mdl_output_pref))
  tar = tarfile.open('{0}/{1}.mdl_pdb.tar.bz2'.format( 
                      result_directory, mdl_output_pref), 'w:bz2')

  with open(mdl_output_pref+'.aligned_pdb.list', 'w') as fo:
    for mdl in Models:
      fo.write(mdl+'\n')
      tar.add(mdl)
  tar.add(mdl_output_pref+'.aligned_pdb.list')
  tar.close()


  #   Copy work-directory items to result-directory
  ### mdl_output_pref.align_pdb.list is hard-coded
  os.system('cp {0} {1}.*pse.bz2 {1}.mod-*.log "{2}"'.format(
            mdl_output_pref+'.align_pdb.list', 
            mdl_output_pref+'.align',
            result_directory))


##########################################################################
##
##  Peter M.U. Ung @ MSSM
##
##  v0.0    ???
##  v1.0    15.12.15
##  v2.0    16.06.12    changed the alignment sequences. For both S/T- and Y-
##                      kinases, use the same template and sequence: catalytic
##                      beta-sheets of C-lobe of 1ATP
##                          1ATP 'resi 122-138+162-183'
##  v2.1    16.08.10    use superpose sequence from global variables
##  v3.0    17.07.11    x_pymol_alignment.py is now Object PyMOLSuperpose
##  v4.0    18.03.09    update for use with CODI multiple runs
##  v5.0    19.12.28    PyMOLSuperose use output PDB extension
##
