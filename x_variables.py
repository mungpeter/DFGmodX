import sys, os
import re

##########################################################################
## Global constants

####  !!Also check x_check_scripts.py to make sure all links are updated!! ####

script_dir    = '/home/pmung/Dropbox/9_scripts/3_program/structures/3_DFGmodx/'
dataset_dir   = script_dir+'z_dataset/'

povme_exec    = '/home/software/POVME-2.1/POVME2.py'  # lagacy; integrated as python function
pymol_exec    = '/usr/bin/pymol'    # PyMOL executable

stkinase_dir  = '/home/pmung/Dropbox/1_kinase/1_family/1_stdkinases/170109/2_align/'
#ykinase_dir   = '/home/pmung/zzz_data/1_kinase/1_family/2_ykinase/170109/2_align/'
#kstruct_dirs  = stkinase_dir+','+ykinase_dir

ident_thres   = 40.0  # seq iden threshold to switch to flag the sequence for check
align_thres   = 5.0  # seq iden threshold to switch from MUSCLE to EXPRESSO alignment

struct_database=dataset_dir+'stdy_kinase_xtal.all.200329.clean.fasta'
struct_nogap  = dataset_dir+'stdy_kinase_xtal.all.200329.clean.nogap.fasta'
kinome_database=dataset_dir+'MD_human_kinome_alignment.2019.200331.be_modeled.fasta'
kinome_nogap  = dataset_dir+'MD_human_kinome_alignment.2019.200331.be_modeled.nogap.fasta'
# T-coffee has a limiting issue with too long of path+filename, result in coredump

reference_pdb = '1atp.pdb'
superpose_resi= 'resi 122-139+162-183'  # C-lobe catalytic beta-sheet: 180709
#superpose_resi= 'resi 126-181+205-305' # entire C-lobe minus activation-loop

#superpose_resi= 'resi 122-138+162-183'  # C-lobe catalytic beta-sheet: 16xxxx
#superpose_resi= 'resi 124-138+160-183'  # C-lobe catalytic beta-sheet: 170203
## 170203 - either C-lobe select is similar, just different set of PDB failed
#superpose_resi= 'resi 122-127+154-175+180-183'	#C-lobe catalytic: 170810
## 170810 - ditch the short helix after hinge region and add helix after cata
##	    loop


##########################################################################
## Reformat FASTA file to have X-number of AA residue per line
def per_line():
  return 72


##########################################################################
## POVME 2 setup for various conformations
def POVMESetup( conformation ):

  Pts = {
    'gen'  : '''
## Used with POVME_2.1, not for POVME_3

PointsInclusionSphere\t\t  6.70   8.30  -1.80  11  # center-general
PointsExclusionSphere\t\t 16.30   9.70  -9.30  7   # back-activation
PointsExclusionSphere\t\t 20.60   9.50   3.50  7   # front-to-side
PointsExclusionSphere\t\t -2.30   7.00  -3.10  6   # side channel
PointsExclusionSphere\t\t -1.50   6.50  -8.90  6   # back-side channel
ContiguousPocketSeedSphere\t  6.70   8.30  -1.80  8   # center-general\n''',

    'cido' : '''
PointsExclusionSphere\t\t 24.50   8.60  -0.20  13  # CIDO-specfic\n''',

    'type-III' : '''
PointsInclusionSphere\t\t  6.70   8.50  -5.50  13  # type-III center
PointsExclusionSphere\t\t -2.30   7.00  -3.10  5   # III side channel
PointsExclusionSphere\t\t -4.15  12.50  -9.45  5   # III back-side channel
PointsExclusionSphere\t\t  8.60   9.50   5.90  8   # III ATP volume
ContiguousPocketSeedSphere\t  6.70   8.50  -5.50  8  # type-III center\n''',

    'other' : '''
ContiguousPointsCriteria\t10\n\nGridSpacing\t\t0.75\nDistanceCutoff\t\t1.09\n
ConvexHullExclusion\t\ttrue\nSaveIndividualPocketVolumes\ttrue
SavePocketVolumesTrajectory\ttrue\nSaveTabbedVolumeFile\t\ttrue
SaveVolumetricDensityMap\ttrue\n\nCompressOutput\t\t\tfalse
UseDiskNotMemory\t\tfalse\nOutputEqualNumPointsPerFrame\tfalse''',

  }

  if not re.search(r'III', conformation, re.IGNORECASE):
    Choice = Pts['gen']
  else:
    Choice = Pts['type-III']

  ## This exclusion sphere also overlaps with General 'back-act' & 'front-side'
  if re.search(r'cido', conformation, re.IGNORECASE):
    Choice += Pts['cido']

  Choice += Pts['other']

  return Choice


##########################################################################
## Choose which conformation template list for homology modeling
def TemplatePDBList ( conf ):
  Template = {
    'cidi_s_templ' : 'u_cidi.stkinase_templ_pdb.txt',
    'cidi_y_templ' : 'u_cidi.ykinase_templ_pdb.txt',
    'cido_s_templ' : 'u_cido.stkinase_templ_pdb.txt',
    'cido_y_templ' : 'u_cido.ykinase_templ_pdb.txt',
    'codi_templ'   : 'u_codi.templ_pdb.1.txt,u_codi.templ_pdb.2.txt,u_codi.templ_pdb.3.txt,u_codi.templ_pdb.4.txt',
    'III_templ'    : 'u_codi.templ_pdb.III.txt',
    'met_templ'    : 'u_codi.templ_pdb.met.txt',
    'egfr_templ'   : 'u_codi.templ_pdb.egfr.txt',
    'codo_templ'   : 'u_codo.templ_pdb.txt',
    'wcd_templ'    : 'u_wcd.templ_pdb.1.txt,u_wcd.templ_pdb.2.txt',
  }
  return Template[conf]


##########################################################################
## CutSite templates for DFGmodel and N-/C-termini. This limits the size
## of modeled structure to specifically the core catalytic domain only
def NCTermTempls():
  NCCutsites = ['1ATP:QLD|QF/NHK|WF',   # (Refernce columns) 1ATP P05132 (bovine)
                '3KQ7:VPE|RY/AQY|HD',   # (c-in,  d-out, st) p38a Q16539
                '3NPC:VLK|RY/RHP|YI',   # (c-in,  d-out, st) JNK2 P45984
                '3UGC:EER|HL/QIR|DN',   # (c-in,  d-out, y)  JAK2 O60674
                '3ETA:SRE|KI/EVS|FF',   # (c-in,  d-out, y)  IGFR P06213
                '3OCS:DPK|DL/VMD|EN',   # (c-out, d-in,  1)  BTK  Q06187 (1t)
                '3VS6:PRE|SL/TAT|ES',   # (c-out, d-in,  1)  HCK  P08631 (1t)
                '2R4B:KET|EL/RDP|QR',   # (c-out, d-in,  2)  ErbB4 Q15303 (mbum)
                '2H8H:PRE|SL/PQY|QP',   # (c-out, d-in,  2)  SRC  P12931 (mbum)
                '4EHG:PDG|QI/RSL|PK',   # (c-out, d-in,  3)  BRaf P15056 (hbum)
                '4RFM:DPS|EL/AIA|AS',   # (c-out, d-in,  3)  ITK  Q08881 (hbum)
                '1AD5:PRE|SL/TAT|ES',   # (c-out, d-in,  4)  HCK  P08631 (norm)
                '3TV6:PDG|QI/RSL|PK',   # (c-out, d-in,  4)  BRaf P15056 (norm)                
                '4M15:DPS|EL/AEI|AE',   # (c-out, d-in,  s)  ITK  Q08881 (mek)
                '3ZLW:KDD|DF/FAG|WL',   # (c-out, d-in,  s)  MEK1 Q02750 (mek)                
                '2WGJ:SSL|IV/FST|FI',   # (c-out, d-in,  s)  cMET P08581 (met)
                '4XMO:SSL|IV/FST|FI',   # (c-out, d-in,  s)  cMET P08581 (met)
                '4I22:KET|EF/RDP|QR',   # (c-out, d-in,  s)  EGFR P00533 (egfr)
                '2C0O:PRE|SL/TAT|ES',   # (c-out, d-in,  s)  HCK  P08631 (egfr)

                '3EFJ:SSL|IV/GEH|YV',   # (c-out, d-out, y)  cMET P08581
                '3NAX:RPE|DF/AHP|FF' ]  # (c-out, d-out, st) PDK1 O15530
  return NCCutsites

def DFGTempls():
  DFGCutsites = [ '1ATP:AGG|EM/GYI|QV/GTP|EY', # (Refernce columns) 1ATP P05132 (bovine)
                  '3KQ7:MGA|DL/CEL|KI/ATR|WY', # (c-in,  d-out, st) p38a Q16539
                  '3NPC:MDA|NL/CTL|KI/VTR|YY', # (c-in,  d-out, st) JNK2 P45984
                  '3UGC:PYG|SL/NRV|KI/SPI|FW', # (c-in,  d-out, y)  JAK2 O60674
                  '3ETA:AHG|DL/FTV|KI/LPV|RM', # (c-in,  d-out, y)  IGFR P06213
                  '3OCS:ANG|CL/GVV|KV/FPV|RW', # (c-out, d-in,  1)  BTK  Q06187 (1t)
                  '3VS6:AKG|SL/LVC|KI/FPI|KW', # (c-out, d-in,  1)  HCK  P08631 (1t)
                  '2R4B:PHG|CL/NHV|KI/MPI|KW', # (c-out, d-in,  2)  ErbB4 Q15303 (mbum)
                  '2H8H:SKG|SL/LVC|KV/FPI|KW', # (c-out, d-in,  2)  SRC  P12931 (mbum)
                  '4EHG:EGS|SL/LTV|KI/GSI|LW', # (c-out, d-in,  3)  BRaf P15056 (hbum)
                  '4RFM:EHG|CL/QVI|KV/FPV|KW', # (c-out, d-in,  3)  ITK  Q08881 (hbum)
                  '1AD5:AKG|SL/LVC|KI/FPI|KW', # (c-out, d-in,  4)  HCK  P08631 (norm)
                  '3TV6:EGS|SL/LTV|KI/GSI|LW', # (c-out, d-in,  4)  BRaf P15056 (norm)                  
                  '4M15:EHG|CL/QVI|KV/FPV|KW', # (c-out, d-in,  s)  ITK  Q08881 (mek)
                  '3ZLW:DGG|SL/GEI|KL/GTR|SY', # (c-out, d-in,  s)  MEK1 Q02750 (mek)                  
                  '2WGJ:KHG|DL/FTV|KV/LPV|KW', # (c-out, d-in,  s)  cMET P08581 (met)
                  '4XMO:KHG|DL/FTV|KV/LPV|KW', # (c-out, d-in,  s)  cMET P08581 (met)
                  '4I22:PFG|CL/QHV|KI/VPI|KW', # (c-out, d-in,  s)  EGFR P00533 (egfr)
                  '2C0O:AKG|SL/LVC|KI/FPI|KW', # (c-out, d-in,  s)  HCK  P08631 (egfr)

                  '3EFJ:KHG|DL/FTV|KV/LPV|KW', # (c-out, d-out, y)  cMET P08581
                  '3NAX:KNG|EL/MHI|QI/GTA|QY'] # (c-out, d-out, st) PDK1 O15530
  return DFGCutsites


##########################################################################
## Default settings
def DefaultVariables():

  home_dir  = str(os.getcwd())
  Variables = {
    'ScriptDirectory'   : script_dir,
    'DatasetDirectory'  : dataset_dir,
    'PymolExecutable'   : pymol_exec,
    'POVMEExecutable'   : povme_exec,
    'HomeDirectory'     : home_dir,
    'ResultDirectory'   : home_dir+'/1_result',
    'WorkingDirectory'  : home_dir+'/0_tempfile',
    'POVMEDirectory'    : home_dir+'/0_tempfile/1_povme',

    'StructDatabase'    : struct_database,
    'StructNoGap'       : struct_nogap,
    'KinomeDatabase'    : kinome_database,
    'KinomeNoGap'       : kinome_nogap,
    'PDBDirectory'      : stkinase_dir,
    'SuperposeRefResi'  : superpose_resi,
    'ReferencePDB'      : reference_pdb,
    'SeqIdentThres'     : ident_thres,
    'AlignSwitchThres'  : align_thres,
    'BestMatchStruc'    : 'None',
    'SeqIdentity'       : 'None',
    'CorrectFASTAFile'  : 'None',

    'KinaseStructInput' : '???/4L43_A.1atp.pdb',
    'ModelKinaseFasta'  : 'e.g. 4L43_A',

    'OutputPrefix'      : 'cido.4L43_L',                    # variable
    'CHelixDFGModel'    : 'cido',                           # variable
    'KinaseFamily'      : 'serine',                         # variable
    'NumberOfModel'     : '50',                             # variable
    'NumberOfTopModel'  : '10',                             # variable
    'NumberOfCPU'       : '10',                             # variable

    'TemplateList'      : 'None',                           # variable
    'ChimeraTemplList'  : '<conf>.<name>.chimera.list',     # variable
    'ModifiedPIRFile'   : '<conf>.<name>.chimera.pir',      # variable
    'PymolAlignPrefix'  : '<conf>.<name>.align',            # variable
    'AlignedModelList'  : '<conf>.<name>.aligned_pdb.list', # variable
    'POVMEStructure'    : '<conf>.<name>.multi_frame.pdb',  # variable
  }
  return Variables


#####################################################################
#
#   v1.0    15.11.02
#   v2.0    16.06.07    add C-out parameters
#   v3.0    16.11.27    add auto structure search parameters
#   v4.0    16.12.05    remove homolog check
#   v5.0    17.02.03    use new sequence alignment 'resi 124-138+160-183'
#   v6.0    17.03.30    internalize PDBDirectory by type of stkinase/ykinase
#   v7.0    17.06.28    add pre-aligned seq profile
#   v8.0    18.03.06    allow modeling for multiple types of conformations,
#                       and multiple subconformations for CODI and wCD
#   v9.0    18.03.27    simplify printout and readin options
#   v10.0   18.03.28    new function to handle variable setup
#   v11.0   18.04.06    include fasta cutsites for kinase conformations
#   v12.0   18.04.08    separate running scripts to "x_variables_run.py"
#   v13.    18.06.15    update fasta library and APE cut-sites
#   v14     20.03.12    add MD_human_kinome_alignment
#   v15     20.04.05    update CutSite residues based on MD kinome
#