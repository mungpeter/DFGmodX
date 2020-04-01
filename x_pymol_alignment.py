#!/usr/bin/env python3
import sys,glob,os,re
from CommonUtility import *

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
#   v7.0    -- 17.07.11 - made into non-callable OOP
#   v8.0    -- 18.06.10 - added best_match structure
#   v9.0    -- 19.12.28   this version came from 3_DFGmodx, new output extension#   v9.1    -- 20.02.27   search "Executive: RMS" for both PyMOL 1.x/2.x
#
#   Purpose: Generate a pymol session file that contains the Template 
#   PDB and the model PDBs. The models are superposed onto the Template
#   using a supplied string of template residues (pymol format).
#   If needed, the superposed models can be output with '.mod.pdb' 
#   extension.
#
##########################################################################

msg = """
    > {0}
        [pymol executable]
        [Template PDB] [Template Residues to align]
        [List of Model PDBs] [Output PyMOL session prefix]
	      [Name of Aligned PDB List]
        [Number of Model] (for Modeller PDB, B999000xx)
        [output PDB extension]\n
  e.g.> x.py /usr/bin/pymol 
          1atp.pdb 'resi 122-138+162-183' model.list align_output _tmp.list 1\n""".format(sys.argv[0])
#if len(sys.argv) != 7: sys.exit(msg)


##########################################################################
def PyMOLSuperpose( pymol_exec, reference_pdb, best_match_struc, 
                    superpose_resi, org_pdb_list, pymol_align_pref, 
                    aligned_mdl_list, number_of_model, out_pdb_ext ):

  print('  \033[34m** PyMOL Structure Superposition **\033[0m')

  # Digit (X) of model number in the Modeller-generated models
  # (50)=2,(500)=3      B9999000XX (2 digits), B999900XXX (3 digits)
  if int(number_of_model) == 0:
    max_digt = 0
  else:
    max_digt = len(str(number_of_model))

  ## Align models to reference structure using 'superimpose'
  AlignStructures(  pymol_exec, aligned_mdl_list, pymol_align_pref, 
                    reference_pdb, best_match_struc, superpose_resi, 
                    org_pdb_list, 'super', max_digt, out_pdb_ext )

  ## Check pymol alignment Log
  Mdls, Atoms = [], []
  # Extract information from pymol-log file
  with open('{0}.{1}.pymol-log'.format(pymol_align_pref,'super'), 'r') as fi:
    for line in fi:
      if re.search(r'PyMOL>load', line):
        curr_pdb = line.split('load ')[1].split(', ')
        Mdls.append(curr_pdb)
      if re.search(r'Executive: RMS', line):
        atom = int(line.split('=')[1].split('(')[1].split('to')[0])
        Atoms.append([line.split(':')[1], atom])
      if re.search(r'Executive: Error', line):
        Atoms.append([line, 0])
        print('{0} | {1}'.format(curr_pdb, line))

  del Mdls[0:2]   # remove the first 2 items, reference and best_match PDBs

  # Write out extracted information and identify bad alignment
  Aligns = []
  with open('{0}.mod-super.log'.format(pymol_align_pref), 'w') as fo:
    # write to alignment, save models that fail to have 'enough atoms'
    for idx, Names in enumerate(Mdls):
      try:
        print(Atoms[idx])
      except IndexError:
        continue
      if Atoms[idx][1] < 60:
        print(' ** Insufficient Number of Atom for Structure Superposition **')
        print('    PDB: '+Names[0])
        print('    Atom: {0} -- fewer than threshold 60'.format(Atoms[idx][1]))
        fo.write(' ** Check {0}: atom < 60\n'.format(Names[0]))
        Aligns.append(Names[0])
      fo.write('{0}: {1}\n'.format(Names[1], Atoms[idx][0]))

  ## Redo the alignment of models that failed the 'superimpose' with 'align'
  if len(Aligns) > 0:
    AlignStructures(  pymol_exec, aligned_mdl_list, pymol_align_pref, 
                      reference_pdb, best_match_struc, superpose_resi, 
                      Aligns, 'align', max_digt, out_pdb_ext )


##########################################################################
## Write out .pml file for PyMOL-based structure superposition
## Rename the aligned models with a shortened model numbering
def AlignStructures(  pymol_exec, aligned_mdl_list, pymol_align_pref, 
                      reference_pdb, best_match_struc, superpose_resi, 
                      org_pdb_list, align_mode, max_digt, out_pdb_ext ):

  pymol_pref = '{0}.{1}'.format(pymol_align_pref, align_mode)
  match_pref = best_match_struc.split('/')[-1].split('.pdb')[0]

  l = open(aligned_mdl_list, 'w')
  r = open('{0}.pml'.format(pymol_pref), 'w')
  r.write('load {0}, reference\n'.format(reference_pdb))
  r.write('sele ref_resid, reference and {0}\n\n'.format(superpose_resi))
  r.write('load {0}, best_match_{1}\n'.format(best_match_struc, match_pref))

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
#      print(max_digt)
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
    l.write(name+'.'+out_pdb_ext+'\n')
    r.write('save {0}, {1}\n'.format(name+'.'+out_pdb_ext, name))

  # wrap up the pymol session
  r.write('hide everything\nshow ribbon\n')
  r.write('color white, all\n')
  r.write('color red, template\n')
  r.write('color cyan, best_match_{0}\n'.format(match_pref))
  r.write('util.cbaw\n')
  r.write('set pse_export_version, 1.70\n')
  r.write('save {0}.pse'.format(pymol_pref))
  r.close()
  l.close()

  os.system('{0} -c {1}.pml > {1}.pymol-log'.format(pymol_exec, pymol_pref))
  os.system('bzip2 -f {0}.pse'.format(pymol_pref))
  print('  ** Finished Superposition **')


##########################################################################
#if __name__ == "__main__":

#  main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],
#       sys.argv[6])
