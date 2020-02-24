#!/usr/bin/env python3
import sys, re,os
import glob
from modeller import *              # Load standard Modeller classes
from modeller.automodel import *    # Load the automodel class
from modeller.parallel import *

## Add Chain_ID to the model 
## If the single chain has breaks, each 'fragment' is considered as a fragment
## and an additional 'self.chains[].name = '?' will be needed.
## The number of chain break is determined by the number of '/' found in the
## .pir alignment file, and the number is added to the AddChainID function 
## via python's 'super' function.
class AddChainID(automodel):

  ## determine how many segment in the final alignment by counting the 
  ## number of '/' 
  def __init__(self, env, alnfile, knowns, sequence):
    super(AddChainID, self).__init__(env, alnfile, knowns, sequence)
    self.seg_num = 1
    with open(alnfile, 'r') as f:
      readLine = False
      counter = 0
      with open('err_pir.output.txt', 'w') as fo:
        for line in f.readlines():
          if not re.search('^#', line):
            if readLine == False:
              if re.search(r'sequence:', line):
                readLine = True
                fo.write(line)
            else:
              self.seg_num += line.count('/')
              fo.write(line.strip()+' '+str(line.count('/'))+' '+str(self.seg_num)+"\n")
      if os.stat('err_pir.output.txt').st_size == 0:
        os.remove('err_pir.output.txt')
    
  def special_patches(self, aln):
    Seg = []
    for x in list(range(0,self.seg_num)):
      Seg.append('A')
    self.rename_segments(segment_ids=Seg)
    for idx in list(range(0,self.seg_num)):
      try:
        self.chains[idx].name = 'A'
      except IndexError:
        sys.exit('\n  -- Error: "AddChainID" class has wrong number of index --\n\n')

 

##########################################################################
# v1.0  -- 14.??.??     copied from sample in Modeller forum
# v2.0  -- 15.07.21     Ryan Smith added '/' checking and class manipulation
# v3.0  -- 16.08.02     remove 'err_pir.output.txt' if no error (filesize = 0b)
# v4.0  -- 16.10.20     ignore marked comment in .pir file
