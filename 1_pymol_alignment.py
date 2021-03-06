#!/usr/bin/env python3

##########################################################################
#
#   Peter M.U. Ung @ MSSM
#	
#   v1.0    -- 14.03.xx
#   v1.1    -- 15.12.15 - change the way models are numbered to avoid
#                         unix's problem with sort order
#   v1.2    -- 16.07.17 - add pymol log output
#   v2.0    -- 16.07.18 - read pymol log and check alignment result
#   v3.0    -- 16.07.19 - rewrite the pymol alignment in OOP format
#                         add superimpose result checking, use 'align'
#                         if 'super' failed
#   v4.0    -- 16.10.27 - bugfix to the B000xxxx search - single digit fail
#   v4.1    -- 16.11.26 - bugfix if file is not model (B9999000xxx)
#   v5.0    -- 17.01.11 - bugfix if Error arose from no atom after refinement
#   v5.1    -- 17.02.03 - bugfix for input filename with 'B99' in it
#   v6.0    -- 17.02.03 - different alignment sequence for better coverage
#                         retire "1atp and resi 122-138+162-183"
#                         newer  "1atp and resi 124-138+160-183"
#                         actually either is similar, not much changed
#   v6.1    -- 17.06.13 - if no badly aligned structure, skip realignment
#
#   Purpose: Generate a pymol session file that contains the Template 
#   PDB and the model PDBs. The models are superposed onto the Template
#   using a supplied string of template residues (pymol format).
#   If needed, the superposed models can be output with '.mod.pdb' 
#   extension.
#
##########################################################################

import sys,glob,os,re
from CommonUtility import *

msg = """
    > {0}
        [Template PDB] [Template Residues to align]
        [List of Model PDBs] [Output PyMOL session prefix]
	[Name of Aligned PDB List]
        [Number of Model] (for Modeller PDB, B999000xx)\n
  e.g.> x.py 1atp.pdb 'resi 122-138+162-183' model.list align_output _tmp.list 1\n""".format(sys.argv[0])
if len(sys.argv) != 7: sys.exit(msg)


##########################################################################
def main( reference_PDB, reference_resid, org_pdb_list, 
          pymol_align_pref, aligned_mdl_list, number_of_model):
  
  print('  ** Running Structure Superposition **')

  # Digit (X) of model number in the Modeller-generated models
  # (50)=2,(500)=3      B9999000XX (2 digits), B999900XXX (3 digits)
  if int(number_of_model) == 0:
    max_digt = 0
  else:
    max_digt = len(str(number_of_model))

  ## Align models to reference structure using 'superimpose'
  AlignStructures( aligned_mdl_list, pymol_align_pref,
                   reference_PDB, reference_resid, org_pdb_list,
                   'super', max_digt)

  ## Check pymol alignment Log
  Mdls, Atoms = [], []
  # Extract information from pymol-log file
  with open('{0}.{1}.pymol-log'.format(pymol_align_pref,'super'), 'r') as fi:
    for line in fi:
      if re.search(r'PyMOL>load', line):
        curr_pdb = line.split('load ')[1].split(', ')
        Mdls.append(curr_pdb)
      if re.search(r'Executive: RMS =', line):
        atom = int(line.split('=')[1].split('(')[1].split('to')[0])
        Atoms.append([line.split(':')[1], atom])
      if re.search(r'Executive: Error', line):
        Atoms.append([line, 0])
        print('{0} | {1}'.format(curr_pdb, line))

  del Mdls[0]   # remove the first item, the reference structure

  # Write out extracted information and identify bad alignment
  Aligns = []
  with open('{0}.mod-super.log'.format(pymol_align_pref), 'w') as fo:
    # write to alignment, save models that fail to have 'enough atoms'
    for idx, Names in enumerate(Mdls):
      print(Atoms[idx])
      if Atoms[idx][1] < 60:
        print(' ** Insufficient Number of Atom for Structure Superposition **')
        print('    PDB: '+Names[0])
        print('    Atom: {0} -- fewer than threshold 60'.format(Atoms[idx][1]))
        fo.write(' ** Check {0}: atom < 60\n'.format(Names[0]))
        Aligns.append(Names[0])
      fo.write('{0}: {1}\n'.format(Names[1], Atoms[idx][0]))

  ## Redo the alignment of models that failed the 'superimpose' with 'align'
  if len(Aligns) > 0:
    AlignStructures( aligned_mdl_list, pymol_align_pref,
                     reference_PDB, reference_resid, Aligns,
                     'align', max_digt)


##########################################################################
## Write out .pml file for PyMOL-based structure superposition
## Rename the aligned models with a shortened model numbering
def AlignStructures( aligned_mdl_list, pymol_align_pref,
                     reference_PDB, reference_resid, org_pdb_list,
                     align_mode, max_digt ):

  pymol_pref = '{0}.{1}'.format(pymol_align_pref, align_mode)
  l = open(aligned_mdl_list, 'w')
  r = open('{0}.pml'.format(pymol_pref), 'w')
  r.write('load {0}, template\n'.format(reference_PDB))
  r.write('sele ref_resid, template and {0}\n\n'.format(reference_resid))

  # Modeller-generated models are name B9999000xxx, ensuring human-sorting
  # Convert model names in a list into array, if not in array already
  if type(org_pdb_list) is str:
    Models = sorted(remove_remark(file_handle(org_pdb_list)))
  else:
    Models = org_pdb_list

  for mdl in Models:
    orig = mdl.split('/')[-1].split('.pdb')[0]

    # If the input structure is not Modeller-generated, pass the name/prefix.
    # If the input structure is Modeller-generated, it will have B999900xxxx,
    # rename the number according to the digit of model to ensure good sorting
    # add '0' to single-digit model naming to match the total number of digits
    # for sorting purpose
    # Some PDB has 'B99': 4B99. Force recognitoin of 'B999' to avoid mix-up
    # e.g. B001 instead of B1 when there are 100 models (max B100)
    if re.search(r'B99[9]+[0]*', orig):     # Modeller PDB has 'B999' prefix
      print(max_digt)
      print(orig)
      if   max_digt == 0:
        name = orig
      elif max_digt == 1:
        num = re.search(r'(.+)B[9]+[0]*([1-9])', orig).group(2)
        name = re.sub(r'(.+)B[9]+[0]*([1-9])', r'\1B'+num, orig)
      else:
        num = ''.join(re.search(r'(.+)B[9]+[0]*([1-9])([0-9]*)', orig).group(2,3))
        if len(str(num)) != max_digt:
          for x in list(range(len(str(num)), max_digt)):
            num = '0'+num
        name = re.sub(r'(.+)B[9]+[0]*([1-9])([0-9]*)', r'\1B'+num, orig)
    else:
      name = orig
 
    # load model and align or superimpose it to reference residues
    r.write('load {0}, {1}\n'.format(mdl, name))
    if align_mode == 'align':
      r.write('align {0}, ref_resid\n'.format(name))
    else:
      r.write('super {0}, ref_resid\n'.format(name))

    # save the aligned model under the new name
    l.write(name+'.mod.pdb\n')
    r.write('save {0}, {1}\n'.format(name+'.mod.pdb', name))

  # wrap up the pymol session
  r.write('hide everything\nshow ribbon\n')
  r.write('color white, all\n color red, template\n')
  r.write('set pse_export_version, 1.70\n')
  r.write('save {0}.pse'.format(pymol_pref))
  r.close()
  l.close()

  os.system('pymol -c {0}.pml > {0}.pymol-log'.format(pymol_pref))
  os.system('gzip -f {0}.pse'.format(pymol_pref))
  print('  ** Finished Structure Superposition **')


##########################################################################
if __name__ == "__main__":
  
  main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],
       sys.argv[6])
