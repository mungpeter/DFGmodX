#!/usr/bin/env python3

import sys
import os,re
import pandas as pd
## this script was used to realign a few sequences to be integrated back into 
## xtal database

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

os.getcwd()

def MuscleProfileAlign( fasta_database, fasta_file, temp_file ):
  os.system('muscle -profile -in1 "{0}" -in2 "{1}" -out "{2}" -maxiters 64 -seqtype protein -gapopen -5.0 -gapextend -2.0 -center 0.0 -quiet'.format(
              fasta_database, fasta_file, temp_file ) )
  Tget_List = list(SeqIO.parse(temp_file, 'fasta'))
  return Tget_List

Data = []
for fas in SeqIO.parse('x_dataset/x.fasta', 'fasta'):
  itm = SeqRecord(  id=fas.id.split('|')[0], seq=fas.seq,
                    description=fas.description, name=fas.id.split('|')[0] )
  Data.append(itm)

fo = open('templ_pdb.list', 'w')
for fas in Data:
  x = fas.id.split('_')
  fo.write('{0} {1}'.format(x[0], x[1]))
  fo.write('\n')
fo.close()

Rst = []
fasta_database = 'x_dataset/MD_human_kinome_alignment.dfgmod.2019.fasta'
for fas in Data:
  with open('temp.{0}.fasta'.format(fas.id), 'w') as fo:
    SeqIO.write(fas, fo, 'fasta')
  fasta_file = 'temp.{0}.fasta'.format(fas.id)
  temp_file  = 'x.{0}.fasta'.format(fas.id)
  Rst.append(MuscleProfileAlign(fasta_database, fasta_file, temp_file))

final = [item for sublist in Rst for item in sublist]
xxx = {}
for fas in final:
  if fas.id not in xxx:
    xxx[fas.id] = fas
aligned = [fas for fas in xxx.values()]
len(aligned)
xxx.keys()
aligned
with open('y.fasta', 'w') as fo:
  SeqIO.write(aligned, fo, 'fasta')

############################################################################
## add uniprot_id to kinome database
os.chdir('/Users/pmung/Dropbox (Schlessinger lab)/9_scripts/3_program/structures/3_DFGmodx/z_dataset')

uni_df = pd.read_csv('kinase_uniprot_id.txt',sep='\s+',comment='#',header=None)
uni_df.columns = ['Name','Uniprot']
uni_df.set_index(['Name'], inplace=True)
unips  = uni_df.to_dict('dict')['Uniprot']; unips['AKT1']

pwd
fastas = [fas for fas in SeqIO.parse('MD_human_kinome_alignment.2019.200324.fasta', 'fasta')]
fastas[5]
# search uniprot for kinase name, rebuild a new fasta file
new_fa = []
nogap_fa = []
failed = []
for fas in fastas:
  temp = fas.id.split('/')
  info = temp[0]
  name = info.split('|')[0]
  uni  = ''
  if name in unips:
    uni = unips[name]
  else:
    for uid in unips.keys():
      if re.search(name, uid, re.IGNORECASE):
        uni = unips[uid]
        break
  if not uni:
    failed.append(name)

  # rebuild a new fasta file with uniprot id
  try:
    new_id = '{0}|{1}/{2}'.format(info, uni, temp[1])
  except IndexError:
    new_id = '{0}|{1}'.format(info, uni)

  itm = SeqRecord(  id=new_id, seq=fas.seq,
                    description='uniprot={0}'.format(uni), name='' )
  new_fa.append(itm)

  nogap = SeqRecord(  id=new_id, seq=Seq(re.sub('-', '', str(fas.seq))),
                    description='uniprot={0}'.format(uni), name='' )
  nogap_fa.append(nogap)

new_fa[10:15]
type(new_fa[0].seq)
type(Seq(nogap_fa[0].seq))
with open('temp.new.fasta', 'w') as fo:
  for fa in new_fa:
    SeqIO.write(fa, fo, 'fasta')
with open('temp.nogap.fasta', 'w') as fo:
  for fa in nogap_fa:
    SeqIO.write(fa, fo, 'fasta')
