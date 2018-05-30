#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys

# stolen from extract_2_3.py
def index2component_general(Nbands, N, ind):
  b=np.zeros((N),dtype=np.int_)
  s=np.zeros((N),dtype=np.int_)
  bs=np.zeros((N),dtype=np.int_)
  # the proposed back conversion assumes the indices are
  # given form 0 to max-1
  ind_tmp = ind - 1
  tmp=np.zeros((N+1),dtype=np.int_)
  for i in xrange(0,N+1):
    tmp[i] = (2*Nbands)**(N-i)

  for i in xrange(N):
    bs[i] = ind_tmp/tmp[i+1]
    s[i] = bs[i]%2
    b[i] = (bs[i]-s[i])//2
    ind_tmp = ind_tmp - tmp[i+1]*bs[i]

  return tuple(bs),tuple(b),tuple(s)

# compute a compound index from orbital indices only.
def component2index_band(Nbands, N, b):
  ind = 1
  for i in xrange(N):
     ind = ind + Nbands**(N-i-1)*b[i]
  return ind

def spinband2index(nbands, bands, spins):
  ''' we expect a 4 digit number for both bands and spins
      each band digit is element of [0,ndim-1]
      while each spin digit is lement of [0,1] = [spinup,spindown] '''
  # b1,b2,b3,b4 = int(str(bands)[0])+1,int(str(bands)[1])+1,int(str(bands)[2])+1,int(str(bands)[3])+1
  # s1,s2,s3,s4 = int(str(spins)[0])+1,int(str(spins)[1])+1,int(str(spins)[2])+1,int(str(spins)[3])+1
  index = 8*nbands**3*(2*(bands[0]+1)+(spins[0]+1)-3) + 4*nbands**2*(2*(bands[1]+1)+(spins[1]+1)-3)+2*nbands*(2*(bands[2]+1)+(spins[2]+1)-3)+2*(bands[3]+1)+(spins[3]+1)-2
  return index

f = h5py.File('threelegs_orig.hdf5','r')
g = h5py.File('threelegs_reversed','w')

a = f['ineq-001'].keys()

for gr in a:
  print(gr)
  bs, b, s = index2component_general(3,4,int(gr))
  b_rev = b[::-1]
  s_rev = s[::-1]
  new_index = spinband2index(3,b_rev, s_rev)
  g['ineq-001/{:05}'.format(new_index)] = f['ineq-001/{:05}'.format(int(gr))][()]

# a = f['ineq-002'].keys()

# for gr in a:
#   print(gr)
#   bs, b, s = index2component_general(3,4,int(gr))
#   b_rev = b[::-1]
#   s_rev = s[::-1]
#   new_index = spinband2index(3,b_rev, s_rev)
#   g['ineq-002/{:05}'.format(new_index)] = f['ineq-002/{:05}'.format(int(gr))][()]
