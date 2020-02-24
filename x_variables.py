import sys, os
import re

##########################################################################
## Global constants

####  !!Also check x_check_scripts.py to make sure all links are updated!! ####

script_dir    = '/home/pmung/Dropbox/9_scripts/3_program/structures/3_DFGmodx'
dataset_dir   = script_dir+'/x_dataset'
pymol_dir     = '/usr/bin/pymol'    # PyMOL executable

stkinase_dir  = '/home/pmung/xxx_data/1_kinase/1_family/1_stkinase/170109/2_align'
ykinase_dir   = '/home/pmung/xxx_data/1_kinase/1_family/2_ykinase/170109/2_align'
kstruct_dirs  = stkinase_dir+','+ykinase_dir

ident_thres   = 40.0  # seq iden threshold to switch to flag the sequence for check
align_thres   = 50.0  # seq iden threshold to switch from MUSCLE to EXPRESSO alignment

fasta_data    = dataset_dir+'/stdy_kinase.clean_3.180815.tcoffee_1d.fasta'
kinase_profile= dataset_dir+'/stdy_kinase.clean_3.180815.tcoffee_1d.nogap.fasta'
# T-coffee has a limiting issue with too long of path+filename, result in coredump

template_pdb  = '1atp.pdb'
superpose_resi= 'resi 122-139+162-183'  # C-lobe catalytic beta-sheet: 180709
#superpose_resi= 'resi 126-181+205-305' # entire C-lobe minus activation-loop

#superpose_resi= 'resi 122-138+162-183'  # C-lobe catalytic beta-sheet: 16xxxx
#superpose_resi= 'resi 124-138+160-183'  # C-lobe catalytic beta-sheet: 170203
## 170203 - either C-lobe select is similar, just different set of PDB failed
#superpose_resi= 'resi 122-127+154-175+180-183'	#C-lobe catalytic: 170810
## 170810 - ditch the short helix after hinge region and add helix after cata
##	    loop

povme_dir     = '/home/software/POVME-2.1/POVME2.py'  # replaced as python function


##########################################################################
## Reformat FASTA file to have X-number of AA residue per line
def per_line():
  return 60


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
## CutSite templates for DFGmodel and N-/C-termini
def NCTermTempls():
  NCCutsites = ['3KQ7:VPE|RY/FAQ|YH',   # (c-in,  d-out, st) p38a
                '3NPC:VLK|RY/ITV|WY',   # (c-in,  d-out, st) jnk2
                '3UGC:EER|HL/QIR|DN',   # (c-in,  d-out, y)  jak2
                '3ETA:SRE|KI/DDL|HP',   # (c-in,  d-out, y)  igfr
                '3OCS:DPK|DL/LDV|MD',   # (c-out, d-in,  1)  BTK Q06187 (1t)
                '3VS6:PRE|SL/DDF|YT',   # (c-out, d-in,  1)  HCK P08631 (1t)
                '2R4B:KET|EL/SRM|AR',   # (c-out, d-in,  2)  ErbB4 Q15303 (mbum)
                '2H8H:PRE|SL/EDY|FT',   # (c-out, d-in,  2)  SRC P12931 (mbum)
                '4EHG:PDG|QI/ELL|AR',   # (c-out, d-in,  3)  BRaf P15056 (hbum)
                '4RFM:DPS|EL/AAI|AA',   # (c-out, d-in,  3)  ITK Q08881 (hbum)
                '1AD5:PRE|SL/DDF|YT',   # (c-out, d-in,  4)  HCK P08631 (norm)
                '3TV6:PDG|QI/ELL|AR',   # (c-out, d-in,  4)  BRaf P015056 (norm)                
                '4M15:DPS|EL/AEI|AE',   # (c-out, d-in,  s)  ITK Q08881 (mek)
                '3ZLW:KDD|DF/IKR|SD',   # (c-out, d-in,  s)  MEK1 Q02750 (mek)                
                '2WGJ:GPS|SL/SAI|FS',   # (c-out, d-in,  s)  cMet P08581 (met)
                '4XMO:GPS|SL/SAI|FS',   # (c-out, d-in,  s)  cMet P08581 (met)
                '4I22:KET|EF/SKM|AR',   # (c-out, d-in,  s)  EGFR P00533 (egfr)
                '2C0O:PRE|SL/DDF|YT',   # (c-out, d-in,  s)  HCK P08631 (egfr)

                '3EFJ:GPS|SL/AIF|ST',   # (c-out, d-out, y)   cmet
                '3NAX:RPE|DF/FES|VT' ]  # (c-out, d-out, s/t) pdk1
  return NCCutsites

def DFGTempls():
  DFGCutsites = [ '3KQ7:MGA|DL/CEL|KI/YVA|TR', # (c-in,  d-out, st) p38a
                  '3NPC:MDA|NL/CTL|KI/YVV|TR', # (c-in,  d-out, st) jnk2
                  '3UGC:PYG|SL/NRV|KI/GES|PI', # (c-in,  d-out, y)  jak2
                  '3ETA:AHG|DL/FTV|KI/GLL|PV', # (c-in,  d-out, y)  igfr
                  '3OCS:ANG|CL/GVV|KV/SKF|PV', # (c-out, d-in,  1)  BTK Q06187 (1t)
                  '3VS6:AKG|SL/LVC|KI/AKF|PI', # (c-out, d-in,  1)  HCK P08631 (1t)
                  '2R4B:PHG|CL/NHV|KI/GKM|PI', # (c-out, d-in,  2)  ErbB4 Q15303 (mbum)
                  '2H8H:SKG|SL/LVC|KV/AKF|PI', # (c-out, d-in,  2)  SRC P12931 (mbum)
                  '4EHG:EGS|SL/LTV|KI/LSG|SI', # (c-out, d-in,  3)  BRaf P15056 (hbum)
                  '4RFM:EHG|CL/QVI|KV/TKF|PV', # (c-out, d-in,  3)  ITK Q08881 (hbum)
                  '1AD5:AKG|SL/LVC|KI/AKF|PI', # (c-out, d-in,  4)  HCK P08631 (norm)
                  '3TV6:EGS|SL/LTV|KI/LSG|SI', # (c-out, d-in,  4)  BRaf P015056 (norm)                  
                  '4M15:EHG|CL/QVI|KV/TKF|PV', # (c-out, d-in,  s)  ITK Q08881 (mek)
                  '3ZLW:DGG|SL/GEI|KL/FVG|TR', # (c-out, d-in,  s)  MEK1 Q02750 (mek)                  
                  '2WGJ:KHG|DL/FTV|KV/AKL|PV', # (c-out, d-in,  s)  cMet P08581 (met)
                  '4XMO:KHG|DL/FTV|KV/AKL|PV', # (c-out, d-in,  s)  cMet P08581 (met)
                  '4I22:PFC|CL/QHV|KI/GKV|PI', # (c-out, d-in,  s)  EGFR P00533 (egfr)
                  '2C0O:AKG|SL/LVC|KI/AKF|PI', # (c-out, d-in,  s)  HCK P08631 (egfr)

                  '3EFJ:KHG|DL/FTV|KV/AKL|PV', # (c-out, d-out, y)   cMet
                  '3NAX:KNG|EL/MHI|QI/FVG|TA'] # (c-out, d-out, s/t) PDPK1
  return DFGCutsites


##########################################################################
## Default settings
def DefaultVariables():

  home_dir  = str(os.getcwd())
  Variables = {
    'ScriptDirectory'   : script_dir,
    'DatasetDirectory'  : dataset_dir,
    'PymolExecutable'   : pymol_dir,
    'POVMELocation'     : povme_dir,
    'HomeDirectory'     : home_dir,
    'ResultDirectory'   : home_dir+'/1_result',
    'WorkingDirectory'  : home_dir+'/0_tempfile',
    'POVMEDirectory'    : home_dir+'/0_tempfile/1_povme',

    'FastaDatabase'     : fasta_data,
    'KinaseProfile'     : kinase_profile,
    'PDBDirectory'      : stkinase_dir,
    'SuperposeRefResi'  : superpose_resi,
    'TemplatePDB'       : template_pdb,
    'SeqIdentThres'     : ident_thres,
    'AlignSwitchThres'  : align_thres,
    'BestMatchStruc'    : '',
    'SeqIdentity'       : '0',
    'CorrectFASTAFile'  : False,

    'OutputPrefix'      : 'cido.4L43_L',                    # User input
    'CHelixDFGModel'    : 'cido',
    'KinaseStructInput' : '???/4L43_A.1atp.pdb',
    'ModelKinaseFasta'  : 'e.g. 4L43_A',
    'KinaseFamily'      : 'serine',
    'NumberOfModel'     : '50',
    'NumberOfTopModel'  : '10',
    'NumberOfCPU'       : '10',

    'TemplateList'      : '',                               # variable
    'POVMETemplateFile' : '',                               # variable
    'ChimeraTemplList'  : 'chimera_xprefix.list',           # variable
    'ModifiedPIRFile'   : 'chimera_xprefix.pir',            # variable
    'PymolAlignPrefix'  : 'xprefix.align',                  # variable
    'AlignedModelList'  : 'xprefix.aligned_pdb.list',       # variable
    'POVMEStructure'    : 'xprefix.multi_frame.pdb',        # variable
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
#
