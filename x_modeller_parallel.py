#!/usr/bin/env python3

# Homology modeling with multiple templates
from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the automodel class
from modeller.parallel import *
from CommonUtility import *
import os,sys

##########################################################################
msg = """\n  > python[or mod9.x] {0}
        [template PDB prefix list]
        [modeller .pir file]
        [number of model to generate]
        [number of CPU]
        [output prefix: same as 'target' in .pir]\n""".format(sys.argv[0])
#if len(sys.argv) != 6: sys.exit(msg)

print('\033[34m** x_modeller_parallel.py **\033[0m')
print('\033[34m## Current directory:\033[0m\n{0}\n'.format(os.getcwd()))

Templates = remove_remark(file_handle(sys.argv[1]))
pir_file  = sys.argv[2]         # .pir file
model_nm  = int(sys.argv[3])    # No. of model to generate
cpu_num   = int(sys.argv[4])    # No. of CPU used for parallel job
tget_prot = sys.argv[5]         # name of 'target' in .pir file

print(tget_prot)


##########################################################################
## Define variables (cont'd)
env = environ()  # create a new MODELLER environment to build this model in

DIR = ['.','..']                # Directories for input PDB files
env.io.hetatm = False		# Read in non-water HETATM records from PDB
env.io.water  = False		# Read in water records from PDB

mdl_start = 1
mdl_end   = model_nm

min_refine = True	# Turn on/off optimization
refine_cyc = 2		# No. of refinement to do


##########################################################################
## Import AddChainID class from modeller to add chain ID.
## Direct imbeddment of the 'class' is not compatible with with parallel job
## Require a separate .py file to store the class
## Shown here is a sample of the code. A more advanced code is in the .py

from x_mod_class import AddChainID
#class AddChainID(automodel):
#  def special_patches(self, aln):
#    self.rename_segments(segment_ids=['A'])
#    self.chains[0].name = 'A'


##########################################################################
## Setup of homology modeling
# request verbose output
log.verbose()

# directories for input atom files
env.io.atom_files_directory = DIR

# Use class to add chain ID to the models (imported earlier)
a = AddChainID( env, alnfile=pir_file, 
                knowns=(Templates), sequence=tget_prot )

# Model refinement and optimization
if min_refine is True:
	# Very thorough VTFM optimization:
	a.library_schedule = autosched.slow
	a.max_var_iterations = 300

	# Thorough MD optimization:
	a.md_level = refine.slow

	# Repeat the whole cycle 2 times and do not stop unless
        # object function > 1E6
	a.repeat_optimization = refine_cyc
	a.max_molpdf = 1e6

# Determine the number of model to generate
a.starting_model= mdl_start
a.ending_model  = mdl_end  

# Use multiple scoring functions to assess models
a.assess_methods = (assess.GA341, assess.DOPE, assess.DOPEHR,
                    assess.normalized_dope)


##########################################################################
## Run Modeller in parallel mode if CPU > 1
if cpu_num > 1:
  j = job()
  for x in range(cpu_num): 
    j.append(local_slave())     # processor no. 1
  a.use_parallel_job(j)


##########################################################################
## do the actual homology modeling
a.make()


##########################################################################
#
#   Peter M.U. Ung @ MSSM
#
#   v1.0    - 14.02.15
#   v2.0    - 15.10.15
#   v3.0    - 16.07.17  - clean up the code further
#
