#! /usr/bin/env python
"""
This script is meant for the extraction of the full vertex __F__
from a component sampled two-particle Green's function
within w2dynamics.
"""

from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
from optparse import OptionParser

__version__ = 0.3
__author__ = "Matthias Pickem"

parser = OptionParser(usage = "usage: python %prog [Options]")
parser.add_option("--1pg", dest="g1", type="string",
                  help="DMFT file - one particle Green's function")
parser.add_option("--2pg", dest="g2", type="string",
                  help="DMFT file - two particle Green's function")
parser.add_option("--type", dest="g2type", type=int, default=2,
                  help="Type of the two particle Green's function (0=Gcon, 1=Chi, 2=G2)")
parser.add_option("--ineq", dest="ineq", type=int, default=1,
                  help="Inequivalent atom from the DMFT calculation")
parser.add_option("--nbands", dest="nbands", type=int, default=1,
                  help="Number of correlated bands on the impurity")
parser.add_option("--bands", dest="bgroup", type=int, default=1111,
                  help="Band combination we want to consider")
parser.add_option("--spins", dest="sgroup", type=int, default=1111,
                  help="Spin combination we want to consider: gets ignored if channel options is activated")
parser.add_option("--channel", dest="channel", type=int, default=0,
                  help="Channel for the sampled group (0=spinband combination, 1=dens/magn)")
parser.add_option("--gamma", dest="gamma", type="string",
                  help="Directly measured gamma file for comparison")
parser.add_option("--offset", dest="offset", type=int,
                  help="the fuck do i know")

def spinband2index(nbands, bands, spins):
  ''' we expect a 4 digit number for both bands and spins
      each band digit is element of [1,ndim]
      while each spin digit is lement of [1,2] = [spinup,spindown] '''
  if bands >= 1111 and bands <= (nbands*1111) and spins >= 1111 and spins <= 2222:
    b1,b2,b3,b4 = int(str(bands)[0]),int(str(bands)[1]),int(str(bands)[2]),int(str(bands)[3])
    s1,s2,s3,s4 = int(str(spins)[0]),int(str(spins)[1]),int(str(spins)[2]),int(str(spins)[3])
    index = 8*nbands**3*(2*b1+s1-3) + 4*nbands*+3*(2*b2+s2-3)+2*nbands*(2*b3+s3-3)+2*b4+s4-2
    return index
  else:
    print('Wrong bands or spin combination')
    sys.exit()

options, args = parser.parse_args()

input_g1 = h5py.File(options.g1,'r')
input_g2 = h5py.File(options.g2,'r')

beta = input_g1[".config"].attrs["general.beta"] # inverse temperature
ineq = options.ineq

giw_file = input_g1["dmft-last/ineq-{:03}/giw/value".format(ineq)][()]
# shape : nband, 2, 2*niw
niw = giw_file.shape[-1]//2
# paramagnetic solution only, average over spin component
# giw = np.zeros((nbands,2*niw), dtype=np.complex128)
giw = np.mean(giw_file, axis=1)

if (options.channel == 0):
  bandspin = spinband2index(options.nbands, options.bgroup, options.sgroup)
  bandspin_rev = spinband2index(options.nbands, int(str(options.bgroup)[::-1]), int(str(options.sgroup)[::-1]))
else:
  bandspin1 = spinband2index(options.nbands, options.bgroup, 1111)
  bandspin2 = spinband2index(options.nbands, options.bgroup, 2222)
  bandspin3 = spinband2index(options.nbands, options.bgroup, 1122)
  bandspin4 = spinband2index(options.nbands, options.bgroup, 2211)
  bandspin5 = spinband2index(options.nbands, options.bgroup, 1221)
  bandspin6 = spinband2index(options.nbands, options.bgroup, 2112)
  bandspin1_rev = spinband2index(options.nbands, int(str(options.bgroup)[::-1]), 1111)
  bandspin2_rev = spinband2index(options.nbands, int(str(options.bgroup)[::-1]), 2222)
  bandspin3_rev = spinband2index(options.nbands, int(str(options.bgroup)[::-1]), 2211)
  bandspin4_rev = spinband2index(options.nbands, int(str(options.bgroup)[::-1]), 1122)
  bandspin5_rev = spinband2index(options.nbands, int(str(options.bgroup)[::-1]), 1221)
  bandspin6_rev = spinband2index(options.nbands, int(str(options.bgroup)[::-1]), 2112)

# BETA HERE
if (options.channel == 0):
  g2_file = beta*input_g2["worm-001/ineq-{:03}/g4iw-worm/{:05}/value".format(ineq,bandspin)][()]
  n4iwf = g2_file.shape[0]//2
  n4iwb = g2_file.shape[-1]//2
else:
  g2_file_1 = beta*input_g2["worm-001/ineq-{:03}/g4iw-worm/{:05}/value".format(ineq,bandspin1)][()]
  g2_file_2 = beta*input_g2["worm-001/ineq-{:03}/g4iw-worm/{:05}/value".format(ineq,bandspin2)][()]
  g2_file_3 = beta*input_g2["worm-001/ineq-{:03}/g4iw-worm/{:05}/value".format(ineq,bandspin3)][()]
  g2_file_4 = beta*input_g2["worm-001/ineq-{:03}/g4iw-worm/{:05}/value".format(ineq,bandspin4)][()]
  g2_file_5 = beta*input_g2["worm-001/ineq-{:03}/g4iw-worm/{:05}/value".format(ineq,bandspin5)][()]
  g2_file_6 = beta*input_g2["worm-001/ineq-{:03}/g4iw-worm/{:05}/value".format(ineq,bandspin6)][()]
  n4iwf = g2_file_1.shape[0]//2
  n4iwb = g2_file_1.shape[-1]//2

# this comes from user input where we start counting from 1
b1 = int(str(options.bgroup)[0])-1
b2 = int(str(options.bgroup)[1])-1
b3 = int(str(options.bgroup)[2])-1
b4 = int(str(options.bgroup)[3])-1

if (options.channel == 0):
  s1 = int(str(options.sgroup)[0])-1
  s2 = int(str(options.sgroup)[1])-1
  s3 = int(str(options.sgroup)[2])-1
  s4 = int(str(options.sgroup)[3])-1

chi0 = np.zeros((2*n4iwf,2*n4iwf,2*n4iwb+1), dtype=np.complex128)
chi0_inv = np.zeros((2*n4iwf,2*n4iwf,2*n4iwb+1), dtype=np.complex128)
# vertical
if (options.channel == 0):
  if b1 == b2 and b3 == b4 and s1 == s2 and s3 == s4:
    g2_file[:,:,n4iwb] -= beta*giw[b1,niw-n4iwf:niw+n4iwf,None]*giw[b3,None,niw-n4iwf:niw+n4iwf]
  for nbos in xrange(2*n4iwb+1):
    chi0[:,:,nbos] = np.diag(beta*giw[b4,niw-n4iwf:niw+n4iwf]*giw[b3,niw-n4iwf+n4iwb-nbos:niw+n4iwf+n4iwb-nbos])
    chi0_inv[:,:,nbos] = np.diag(-1.0/(beta*giw[b4,niw-n4iwf:niw+n4iwf]*giw[b3,niw-n4iwf+n4iwb-nbos:niw+n4iwf+n4iwb-nbos]))
    if b1 == b4 and b2 == b3 and s1 == s4 and s2 == s3:
      g2_file[:,:,nbos] += chi0[:,:,nbos]

    g2_file[:,:,nbos] = g2_file[:,:,nbos].dot(chi0_inv[:,:,nbos])

  gamma_sum = np.sum(g2_file, axis=0)
else:
  g2_dens = (g2_file_1 + g2_file_2 + g2_file_3 + g2_file_4)/2.0
  g2_magn = (g2_file_1 + g2_file_2 - g2_file_3 - g2_file_4 + g2_file_5 + g2_file_6 )/4.0
  if b1 == b2 and b3 == b4:
    g2_dens[:,:,n4iwb] -= 2*beta*giw[b1,niw-n4iwf:niw+n4iwf,None]*giw[b3,None,niw-n4iwf:niw+n4iwf]
  for nbos in xrange(2*n4iwb+1):
  #   pass
    # this diagonal routines makes sure that nu = nu'
    # further we define the 2 right side legs
    chi0[:,:,nbos] = np.diag(beta*giw[b4,niw-n4iwf:niw+n4iwf]*giw[b3,niw-n4iwf+n4iwb-nbos:niw+n4iwf+n4iwb-nbos])
    chi0_inv[:,:,nbos] = np.diag(-1.0/(beta*giw[b4,niw-n4iwf:niw+n4iwf]*giw[b3,niw-n4iwf+n4iwb-nbos:niw+n4iwf+n4iwb-nbos]))
    if b1 == b4 and b2 == b3: # horizontal only for uuuu
      g2_dens[:,:,nbos] += chi0[:,:,nbos] # this has the spin configuration for this part
      g2_magn[:,:,nbos] += chi0[:,:,nbos] # this has the spin configuration for this part

    g2_dens[:,:,nbos] = g2_dens[:,:,nbos].dot(chi0_inv[:,:,nbos])
    g2_magn[:,:,nbos] = g2_magn[:,:,nbos].dot(chi0_inv[:,:,nbos])

  gamma_sum_d = np.sum(g2_dens, axis=0)
  gamma_sum_m = np.sum(g2_magn, axis=0)


# COMPARISON WITH INTERNAL SUMMED CHI0F object

input_internal = h5py.File('/home/lv70585/pickemERC/SRVO-VQ-TESTING/1imp-test/output-internal/adga-20180525-161024.223-output.hdf5','r')
# input_internal = h5py.File('/home/lv70585/pickemERC/SRVO-VQ-TESTING/1imp-test/output-gamma+chi0/adga-20180524-152953.147-output.hdf5','r')
# input_internal = h5py.File('/home/lv70585/pickemERC/SRVO-VQ-TESTING/1imp-test/output-gammachi0/adga-20180524-154517.677-output.hdf5','r')






# COMPARISON WITH DIRECTLY MEASURED GAMMA
input_gamma = h5py.File(options.gamma,'r')
if (options.channel == 0):
  gamma = beta*input_gamma["worm-001/ineq-{:03}/p3iw-worm/{:05}/value".format(ineq,bandspin_rev)][()]
  n3iwf = gamma.shape[0]//2
  n3iwb = gamma.shape[-1]//2
else:
  gamma_1 = beta*input_gamma["worm-001/ineq-{:03}/p3iw-worm/{:05}/value".format(ineq,bandspin1_rev)][()]
  gamma_2 = beta*input_gamma["worm-001/ineq-{:03}/p3iw-worm/{:05}/value".format(ineq,bandspin2_rev)][()]
  gamma_3 = beta*input_gamma["worm-001/ineq-{:03}/p3iw-worm/{:05}/value".format(ineq,bandspin3_rev)][()]
  gamma_4 = beta*input_gamma["worm-001/ineq-{:03}/p3iw-worm/{:05}/value".format(ineq,bandspin4_rev)][()]
  gamma_5 = beta*input_gamma["worm-001/ineq-{:03}/p3iw-worm/{:05}/value".format(ineq,bandspin5_rev)][()]
  gamma_6 = beta*input_gamma["worm-001/ineq-{:03}/p3iw-worm/{:05}/value".format(ineq,bandspin6_rev)][()]
  n3iwf = gamma_1.shape[0]//2
  n3iwb = gamma_2.shape[-1]//2

occ = input_g1["dmft-last/ineq-{:03}/occ/value".format(ineq)][()]
dens = np.diagonal(np.diagonal(occ,0,-4,-2),0,-3,-2)

print(dens)
if (options.channel == 0):
  if b1 == b2 and b3 == b4 and s1 == s2 and s3 == s4:
    gamma[:,n3iwb] += beta*giw[b3,niw-n3iwf:niw+n3iwf]*(1.0-dens[b1,s1])

  #horizontal
  if b1 == b4 and b3 == b2 and s1 == s4 and s3 == s2:
    for nbos in xrange(2*n3iwb+1):
      gamma[:,nbos] += giw[b1,niw-n3iwf:niw+n3iwf]*giw[b2,niw-n3iwf+n3iwb-nbos:niw+n3iwf+n3iwb-nbos]

  for nbos in xrange(2*n3iwb+1):
    gamma[:,nbos] = gamma[:,nbos] * (-1.0/(beta*giw[b1,niw-n3iwf:niw+n3iwf]*giw[b2,niw-n3iwf+n3iwb-nbos:niw+n3iwf+n3iwb-nbos]))
else:
  gamma_d = (gamma_1 + gamma_2 + gamma_3 + gamma_4)/2.0
  gamma_m = (gamma_1 + gamma_2 - gamma_3 - gamma_4 + gamma_5 + gamma_6)/4.0
  # gamma_m = (gamma_5 + gamma_6)/2.0

  # vertical
  if b1 == b2 and b3 == b4: #and s1 == s2 and s3 == s4:
    gamma_d[:,n3iwb] += 2*beta*giw[b3,niw-n3iwf:niw+n3iwf]*(1.0-dens[b1,0])

  #horizontal
  if b1 == b4 and b3 == b2:
    for nbos in xrange(2*n3iwb+1):
      gamma_d[:,nbos] += giw[b1,niw-n3iwf:niw+n3iwf]*giw[b2,niw-n3iwf+n3iwb-nbos:niw+n3iwf+n3iwb-nbos]
      gamma_m[:,nbos] += giw[b1,niw-n3iwf:niw+n3iwf]*giw[b2,niw-n3iwf+n3iwb-nbos:niw+n3iwf+n3iwb-nbos]

  for nbos in xrange(2*n3iwb+1):
    gamma_d[:,nbos] = gamma_d[:,nbos] * (-1.0/(beta*giw[b1,niw-n3iwf:niw+n3iwf]*giw[b2,niw-n3iwf+n3iwb-nbos:niw+n3iwf+n3iwb-nbos]))
    gamma_m[:,nbos] = gamma_m[:,nbos] * (-1.0/(beta*giw[b1,niw-n3iwf:niw+n3iwf]*giw[b2,niw-n3iwf+n3iwb-nbos:niw+n3iwf+n3iwb-nbos]))


# COMPARISON WITH EXTRACT 2 3 SCRIPT

# extract = h5py.File('threelegs_reversed.hdf5','r')
# extract_1 = extract["ineq-{:03}/{:05}".format(ineq,bandspin1)][()]
# extract_2 = extract["ineq-{:03}/{:05}".format(ineq,bandspin2)][()]
# extract_3 = extract["ineq-{:03}/{:05}".format(ineq,bandspin3)][()]
# extract_4 = extract["ineq-{:03}/{:05}".format(ineq,bandspin4)][()]
# extract_5 = extract["ineq-{:03}/{:05}".format(ineq,bandspin5)][()]
# extract_6 = extract["ineq-{:03}/{:05}".format(ineq,bandspin6)][()]

# extract_d = (extract_1 + extract_2 + extract_3 + extract_4)/2.0
# extract_m = (extract_1 + extract_2 - extract_3 - extract_4 + extract_5 + extract_6)/4.0


# threelegs = h5py.File('/home/lv70585/pickemERC/SRVO-VQ-TESTING/1imp-test/output-p3/adga-20180524-191731.992-output.hdf5','r')
# threelegs = h5py.File('/home/lv70585/pickemERC/SRVO-VQ-TESTING/1imp-test/output-p3/adga-20180524-194534.162-output.hdf5','r')
# threelegs = h5py.File('/home/lv70585/pickemERC/SRVO-VQ-TESTING/1imp-test/output-external/adga-20180525-113535.105-output.hdf5','r')
# threelegs = h5py.File('/home/lv70585/pickemERC/SRVO-VQ-TESTING/1imp-test/output-external/adga-20180525-120011.798-output.hdf5','r')
# threelegs = h5py.File('/home/lv70585/pickemERC/SRVO-VQ-TESTING/1imp-test/output-external/adga-20180525-131753.186-output.hdf5','r')
# threelegs = h5py.File('/home/lv70585/pickemERC/SRVO-VQ-TESTING/1imp-test/output-external/adga-20180525-133845.887-output.hdf5','r')
# threelegs = h5py.File('/home/lv70585/pickemERC/SRVO-VQ-TESTING/1imp-test/output-external/adga-20180525-141305.583-output.hdf5','r')
# threelegs = h5py.File('/home/lv70585/pickemERC/SRVO-VQ-TESTING/1imp-test/output-external/adga-20180525-142357.148-output.hdf5','r')
# threelegs = h5py.File('/home/lv70585/pickemERC/SRVO-VQ-TESTING/1imp-test/output-external/adga-20180525-144132.233-output.hdf5','r')
# threelegs = h5py.File('/home/lv70585/pickemERC/SRVO-VQ-TESTING/1imp-test/output-external/adga-20180525-145447.696-output.hdf5','r')
threelegs = h5py.File('/home/lv70585/pickemERC/SRVO-VQ-TESTING/1imp-test/output-external/adga-20180525-161017.908-output.hdf5','r')


# g = plt.figure()
# # plt.plot(beta*gamma_m[n3iwf-n4iwf:n3iwf+n4iwf,n3iwb].real, label='internal')
# plt.plot(threelegs['gamma/magn'][b1,b2,b3,b4][:,n4iwb+options.offset].real, label='external_m')
# plt.plot(input_internal['gamma/magn'][b1,b2,b3,b4,:,n4iwb+options.offset].real, label='internal_m')
# plt.plot(threelegs['gamma/dens'][b1,b2,b3,b4][:,n4iwb+options.offset].real, label='external_d')
# plt.plot(input_internal['gamma/dens'][b1,b2,b3,b4,:,n4iwb+options.offset].real, label='internal_d')
# plt.legend(loc='best')
# plt.show()
# sys.exit()


if True:
  g = plt.figure(figsize=(15,7))
  ax1 = g.add_subplot(121)
  ax1.set_title(r'$P3$')
  plt.pcolormesh(beta*gamma[n3iwf-n4iwf:n3iwf+n4iwf,n3iwb-n4iwb:n3iwb+n4iwb+1].real)
  vmin, vmax = plt.gci().get_clim()
  plt.colorbar()

  ax2 = g.add_subplot(122)
  ax2.set_title(r'$G4$')
  plt.pcolormesh(gamma_sum[:,:].real, vmin=vmin, vmax=vmax)
  plt.colorbar()

  plt.show()

if False:
# above here ... everything works properly
  g=plt.figure(figsize=(25,10))

  ax5 = g.add_subplot(243)
  ax5.set_title(r'$\gamma_d$')
  plt.pcolormesh(beta*gamma_d[n3iwf-n4iwf:n3iwf+n4iwf,n3iwb-n4iwb:n3iwb+n4iwb+1].real)
  vmind, vmaxd = plt.gci().get_clim()
  plt.colorbar()

  ax6 = g.add_subplot(247)
  ax6.set_title(r'$\gamma_m$')
  plt.pcolormesh(beta*gamma_m[n3iwf-n4iwf:n3iwf+n4iwf,n3iwb-n4iwb:n3iwb+n4iwb+1].real)
  vminm, vmaxm = plt.gci().get_clim()
  plt.colorbar()

  ax1 = g.add_subplot(241)
  ax1.set_title(r'$\gamma_d$')
  plt.pcolormesh(gamma_sum_d.real,vmin=vmind,vmax=vmaxd)
  plt.colorbar()

  ax2 = g.add_subplot(245, sharex=ax1, sharey=ax1)
  ax2.set_title(r'$\gamma_m$')
  plt.pcolormesh(gamma_sum_m.real,vmin=vminm,vmax=vmaxm)
  plt.colorbar()

  ax3 = g.add_subplot(242)
  ax3.set_title(r'$\gamma_d$')
  plt.pcolormesh(input_internal['gamma/dens'][b1,b2,b3,b4,:,:].real, vmin=vmind, vmax=vmaxd)
  plt.colorbar()

  ax4 = g.add_subplot(246)
  ax4.set_title(r'$\gamma_m$')
  plt.pcolormesh(input_internal['gamma/magn'][b1,b2,b3,b4,:,:].real, vmin=vminm, vmax=vmaxm)
  plt.colorbar()

  ax5 = g.add_subplot(244)
  ax5.set_title(r'$\gamma_d$')
  plt.pcolormesh(threelegs['gamma/dens'][b1,b2,b3,b4].real)
  plt.colorbar()

  ax2 = g.add_subplot(248)
  ax2.set_title(r'$\gamma_m$')
  plt.pcolormesh(threelegs['gamma/magn'][b1,b2,b3,b4].real)
  plt.colorbar()

  plt.show()
