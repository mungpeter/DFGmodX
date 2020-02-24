#!/usr/bin/env python3

import os,sys,re,glob
from Bio import SeqIO
from CommonUtility import file_handle
from x_variables import TemplatePDBList

msg = '''\n\t> {0}
\t  [ fasta file with Kinase name ] [ conformation ] [ output prefix ]\n
\t  [ conformation: cidi -- cidi_s_templ or cidi_y_templ ]
\t  [               cido -- cido_s_templ or cido_y_templ ]
\t  [               codi -- codi_templ                   ]
\t  [               codo -- codo_templ                   ]
\t  [               wcd  -- wcd_templ                    ]\n'''.format(sys.argv[0])
if len(sys.argv) != 4: sys.exit(msg)

##########################################################################
## Get protein name from Fasta database
def main( fasta_inp, conformation, out_pref ):

  Names = [seq.id.split()[0].split('/')[1] for seq in SeqIO.parse(file_handle(fasta_inp), 'fasta')]
  print('\n# Number of kinase in input list: '+str(len(Names)))

  cluster = len(TemplatePDBList(conformation).split(','))
  CheckPIRColumn(Names, conformation, out_pref)
  CheckPIRFiles(Names, cluster, out_pref)

#  PrintOutput( Names, ExtractPIRLog(Names), out_pref )


##########################################################################
## Check if all pir(sub-cluster) for a conformation are generated
def CheckPIRFiles( Names, cluster, out_pref ):

  fo  = open(out_pref+'.pir-gen.txt', 'w')
  err = 0

  for name in Names:
    if not os.path.exists(name):
      print('# Protein directoy not found: '+name)
      continue

    Preps = glob.glob('{0}/1_result/*{0}*.pir*'.format(name))

    if not len(Preps):
      fo.write('xx no .pir generated: {0}\n'.format(name))
      err += 1
    elif len(Preps) == cluster:
#      fo.write('.. all .pir generated: {0}\n'.format(name))
      continue
    else:
      fo.write('## not all .pir are generated: {0}\n'.format(name))
      err += 1
      Run = []
      for prep in Preps:
        try:
          Run.append( int(prep.split('.')[-1]) )
        except TypeError:
          fo.write('** .pir may be wrong: {0}\n'.format(name))
          break
      Pir = sorted(Run)
      fo.write('  -- {0}\n'.format(Pir) )

      Tmpl = [x for x in range(1,cluster+1)]
      Diff = [x for x in Pir+Tmpl if x not in Tmpl or x not in Pir ]
      for miss in Diff:
        fo.write('  -- missing this: {0}\n'.format(miss))
  
  print('# Completion:\t'+str(len(Names)-err))
  print('# Error found:\t'+str(err)+'\n')

  fo.close()


##########################################################################
# Check column number of each seq. If target seq has more/fewer column than
# others, it is wrong and need fix
def CheckPIRColumn( Names, conformation, out_pref ):

  print('### Check _TEMP.***y1.fasta pir length ###')
  fo = open(out_pref+'.chk_column.txt', 'w')
  fo.write('## target sequence differs from aligned template sequences ##\n')

  for name in Names:
    for fasta in sorted(glob.glob('{0}/0_tempfile/_T*{0}*y1*fasta'.format(name))):

      if os.path.isfile(fasta):
        Ali  = SeqIO.parse(fasta, 'fasta')
        Leng = [ len(a.seq) for a in Ali ]
        diff = Leng[-1] - Leng[0]
        if diff:
          fo.write('{0}\t{1}\t{2}\n'.format(name, diff, fasta))


##########################################################################
## Look into .pir_log output to extract Template selection, seq identity
## of template to input fasta, and alignment data.
def ExtractPIRLog( Names ):
  Output = []
  for name in Names:
    if not os.path.exists(name):
      print('# Protein directory not found: '+name)
      continue
    with file_handle(name+'/'+name+'.pir_log') as fi:
      File = [l for l in fi]

    template, Identity, Alignment = None, [], []
    for idx, l in enumerate(File):
      if re.search(r'as structure template', l):
        template = l.split()[2]
        Identity  = File[idx:idx+3]
      if re.search(r'Check the Target FASTA', l):
        Alignment = File[idx+4:-1]
        break

    if len(Identity) == 0:
      Identity = ['None']
      pc_ident = 0.0
    else:
      for l in Identity:
        if re.search(r'Percent Identity:', l, re.IGNORECASE):
          pc_ident = float(l.split()[3])

    if len(Alignment) == 0:
      Alignment = ['None\n']

    Output.append([name, template, pc_ident, Identity, File[idx+4:-1]])

  print('# Number of kinase directory found: '+str(len(Output)))
  return Output


##########################################################################
## Output the kinase information
def PrintOutput( Names, Output, out_pref ):
  with open(sys.argv[2]+'.txt', 'w') as fo:
    count_l, count_m, count_h = 0, 0, 0
    Low_Idt, Mid_Idt, Hig_Idt = [], [], []
    for idx, itm in enumerate(Output):
      if itm[2] < 40.0:
        Low_Idt.append(itm)
        print('{0:8s}\t{1:16s}\t{2:5.1f}'.format(itm[0], itm[1], itm[2]))
        fo.write('  > ##ll '+itm[0]+'\n\n')
        count_l = count_l+1
      elif itm[2] < 50.0:
        Mid_Idt.append(itm)
        fo.write('  > ##mm '+itm[0]+'\n\n')
        count_m = count_m+1
      else:
        Hig_Idt.append(itm)
        fo.write('##hh '+itm[0]+'\n\n')
        count_h = count_h+1

      for l in itm[3]:
        fo.write(l)
      fo.write('\n')
      for l in itm[4]:
        fo.write(l)
      fo.write('\n######################################################\n\n')

  with open(out_pref+'.hig.txt', 'w') as fo:
    fo.write('# Number of kinase:\t'+str(len(Names))+'\n')
    fo.write('# Total mid-range (> 50%):\t'+str(count_h)+'\n')
    for itm in sorted(Hig_Idt, key=lambda tup:tup[2], reverse=True):
      fo.write('{0:8s}\t{1:16s}\t{2:5.1f}\n'.format(itm[0], itm[1], itm[2]))

  with open(out_pref+'.mid.txt', 'w') as fo:
    fo.write('# Number of kinase:\t'+str(len(Names))+'\n')
    fo.write('# Total mid-range (40-50):\t'+str(count_m)+'\n')
    for itm in sorted(Mid_Idt, key=lambda tup:tup[2], reverse=True):
      fo.write('{0:8s}\t{1:16s}\t{2:5.1f}\n'.format(itm[0], itm[1], itm[2]))

  with open(out_pref+'.low.txt', 'w') as fo:
    fo.write('# Number of kinase:\t'+str(len(Names))+'\n')
    fo.write('# Total low-range (< 40%):\t'+str(count_l)+'\n')
    for itm in sorted(Low_Idt, key=lambda tup:tup[2], reverse=True):
      fo.write('{0:8s}\t{1:16s}\t{2:5.1f}\n'.format(itm[0], itm[1], itm[2]))
  
  print('# Total low-range ident ( < 40%): '+str(count_l))
  print('# Total mid-range ident (40-50%): '+str(count_m))
  print('# Total hig-range ident ( > 50%): '+str(count_h))


##########################################################################
if __name__ == '__main__':

  main( sys.argv[1], sys.argv[2], sys.argv[3] )

##########################################################################
#
#   Peter M.U. Ung @ MSSM
#
#   v1.0  17.07.14
#   v2.0  18.03.30    extract the sequence identity and highlight it in result
#   v3.0  18.04.11    check pir completion for conformation with sub-clusters
#   v4.0  18.08.21    check if seq in _TEMP.***y1.fasta|_Temp.***y1.corr.fasta
#                     have same number of column
#
#   Purpose: Extract sequence identity information from directories for review
#
#
