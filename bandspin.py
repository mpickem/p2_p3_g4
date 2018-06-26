#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import h5py
import sys
from optparse import OptionParser

def spinband2index(nbands, bands, spins):
    b1,b2,b3,b4=bands[0]+1,bands[1]+1,bands[2]+1,bands[3]+1
    s1,s2,s3,s4=spins[0]+1,spins[1]+1,spins[2]+1,spins[3]+1
    index = 8*nbands**3*(2*b1+s1-3) + 4*nbands**2*(2*b2+s2-3)+2*nbands*(2*b3+s3-3)+2*b4+s4-2
    return index

""" Script to create python list of bandspin indices
for w2dynamics component sympling """

parser = OptionParser()
(options, args) = parser.parse_args()

if (len(args) != 1):
  print("Script Usage: <script> #bands")
  sys.exit()

ndim = int(args[0])

mylist = []

for i in xrange(ndim):
  for j in xrange(ndim):
    for k in xrange(ndim):
      for l in xrange(ndim):

        # kanamori interaction; xxxx xxyy xyyx xyxy
        if (i==j and k==l): #xxyy
          pass
        elif (i==l and j==k): #xyyx
          pass
        elif (i==k and j==l): #xyxy
          pass
        else:
          continue

        for s1 in xrange(2):
          for s2 in xrange(2):
            for s3 in xrange(2):
              for s4 in xrange(2):

                if (s1==s2 and s3==s4): # sigma sigma'
                  pass
                elif (s1==s4 and s2==s3): # sigma sigma' quer
                  pass
                else:
                  continue

                mylist.append(spinband2index(ndim,tuple((i,j,k,l)),tuple((s1,s2,s3,s4))))


mylist.sort()
np.savetxt('indices.dat', mylist, fmt='%i')
print(mylist)
print('Kanamori max: ',2**4*ndim**4)
