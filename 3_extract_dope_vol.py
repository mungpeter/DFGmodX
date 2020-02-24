#!/usr/bin/env python3

import re,sys,os
import numpy as np
import pandas as pd
from plotnine import *
from Bio import SeqIO

msg = '''> python3 {0}
\t[Fasta Database]\n\t[Prefix: cidi|cido|codi|codo|wcd]\n\t[Top Number of Volume]
\t[Output Prefix]
'''.format(sys.argv[0])
if len(sys.argv) != 5: sys.exit(msg)

##########################################################################
def main( fasta_database, prefix, top_modl, out_pref ):
  all_info = {}   # { kinase_name: [volumes], [zDOPE] }
  
  home_dir = os.getcwd()
  Info = ExtractData( home_dir, fasta_database, prefix, top_modl )
  os.chdir(home_dir)

  df_all = GenerateStatistics( Info, out_pref )

  GeneratePlots( df_all, out_pref )


##########################################################################
def GenerateStatistics( Info, out_pref ):
  all_vols, all_zdope, top_vols, top_zdope, Names = Info

  # [data] is equivalent to .T to df, keep kinase name as column name
  df_vols   = pd.DataFrame(all_vols)
  df_zdope  = pd.DataFrame(all_zdope)
  df_tvols  = pd.DataFrame(top_vols)
  df_tzdope = pd.DataFrame(top_zdope)

  vol_avg    = df_vols.median(axis=0)
  zdope_avg  = df_zdope.median(axis=0)
  tvol_avg   = df_tvols.median(axis=0)
  tzdope_avg = df_tzdope.median(axis=0)

  vol_std    = df_vols.std(axis=0)
  zdope_std  = df_zdope.std(axis=0)
  tvol_std   = df_tvols.std(axis=0)
  tzdope_std = df_tzdope.std(axis=0)

  df_all = pd.concat( [vol_avg, vol_std,  zdope_avg, zdope_std,
                       tvol_avg,tvol_std,tzdope_avg,tzdope_std],
                       axis=1, sort=False )
  df_all.columns = ['vol_avg','vol_std','zdope_avg','zdope_std',
                    'top_vol_avg','top_vol_std','top_zdope_avg','top_zdope_std']
  df_new = df_all.reset_index()
  df_new.rename(columns={'index':'kinase'}, inplace=True)

  df_new.to_csv(out_pref+'.csv', sep=',', header=True, index=False )
  df_new.to_excel(out_pref+'.xlsx', index=False)

  return(df_new)


###########################################################################
def GeneratePlots( df_all, out_pref ):

  df_kinase = list(df_all['kinase'])
  df_index  = list(df_all.index)

  df_vol    = df_all[['kinase', 'vol_avg', 'vol_std']].rename(columns={'vol_avg':'avg', 'vol_std':'std'})
  df_zdope  = df_all[['kinase', 'zdope_avg', 'zdope_std']].rename(columns={'zdope_avg':'avg', 'zdope_std':'std'})
  df_tvol   = df_all[['kinase', 'top_vol_avg', 'top_vol_std']].rename(columns={'top_vol_avg':'avg', 'top_vol_std':'std'})
  df_tzdope = df_all[['kinase', 'top_zdope_avg', 'top_zdope_std']].rename(columns={'top_zdope_avg':'avg', 'top_zdope_std':'std'})
  df_vol['group']    = 'all_vol'
  df_zdope['group']  = 'all_zdope'
  df_tvol['group']   = 'top_vol'
  df_tzdope['group'] = 'top_zdope'

  df = pd.concat( [df_vol, df_zdope, df_tvol, df_tzdope] )
  df_v = df[df['group'].str.contains('vol')].reset_index()
  df_z = df[df['group'].str.contains('zdope')].reset_index()

  a = ( ggplot(df_v, aes(x='index', y='avg', color='group')) +
        geom_line() + geom_point() + xlab('kinase') + ylab('Volume') +
        geom_errorbar(aes(ymin="avg-std", ymax="avg+std"), width=.1) +
        scale_x_continuous(breaks=df_index, labels=df_kinase) +
        theme(axis_text_x=element_text(rotation=45, hjust=1)) )  
  ggsave(plot=a, filename=out_pref+'.vol.png', dpi=300 )

  b = ( ggplot(df_z, aes(x='index', y='avg', color='group')) +
        geom_line() + geom_point() + xlab('kinase') + ylab('zDOPE')  +
        geom_errorbar(aes(ymin="avg-std", ymax="avg+std"), width=.1) +
        scale_x_continuous(breaks=df_index, labels=df_kinase) +
        theme(axis_text_x=element_text(rotation=45, hjust=1)) )  
  ggsave(b, filename='{0}.zdope.png'.format(out_pref), dpi=300 )


###########################################################################
def ExtractData( home_dir, fasta_database, prefix, top_modl ):
  all_vols, all_zdope, top_vols, top_zdope = {}, {}, {}, {}
  Names = []
  for idx, fasta in enumerate(SeqIO.parse(fasta_database, 'fasta')):
    family, name = fasta.id.split()[0].split('/')
    Names.append(name)

    # go in kinase directory
    if os.path.isdir(home_dir+'/'+name+'/1_result'):
      os.chdir(home_dir+'/'+name+'/1_result')
    else:
      print('  WARNING: {0} not found'.format(home_dir+'/'+name+'/1_result'))
      continue

    pname = '{0}.{1}'.format(prefix, name)
 
    if os.path.isfile('{0}.volumes.sorted.txt'.format(pname)):
      Vols  = [x.split() for x in open('{0}.volumes.sorted.txt'.format(pname), 'r')]
    else:
      print('  WARNING: {0} not found'.format('{0}.volumes.sorted.txt'.format(pname)))
      continue

    ZdopeDic, Zdope = {}, []
    if os.path.isfile('{0}.zDOPE.txt'.format(pname)):
      for itm in [x.split() for x in open('{0}.zDOPE.txt'.format(pname), 'r') if re.search(r'.pdb', x)]:
        x = int(re.sub('B9999', '', itm[0].split('.')[-2]))
        ZdopeDic[str(x)] =  float(itm[-1])
        Zdope.append(float(itm[-1]))
    else:
      print('  WARNING: {0} not found'.format('{0}.zDOPE.txt'.format(pname)))
      continue
    
    Top_vols  = Vols[0:top_modl]
    Top_zdope = [ ZdopeDic[model[0]] for model in Top_vols ]

    all_vols[name]  = np.array(list(zip(*Vols))[1]).astype(float)
    all_zdope[name] = Zdope 
    top_vols[name]  = np.array(list(zip(*Top_vols))[1]).astype(float)
    top_zdope[name] = Top_zdope

  return(all_vols, all_zdope, top_vols, top_zdope, Names)


##########################################################################
if __name__ == '__main__':
  main( sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4] )
