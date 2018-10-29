import numpy as np
import NFWtriaxMCmodel
import sys
#sys.path.append('/Users/brandyn/astro/')
import dataHandle
import cPickle
#sys.path.append('/Users/brandyn/astro/baryonimpact')
#import Clusters

r200true = 2.
ctrue = 4.
atrue = .99
btrue = 1.
thetatrue = 0.
phitrue = 0.

#zl = 0.25
#zs=1.
ndens=10
sige = 0.25
fov = 10.

nwalkers = 20
lenburn = 200
lenchain = 500
npix = 100

#triax_file = open('./testcats/triax_spherical_bf.txt','r')
#triax_file = open("/Users/brandyn/astro/triax/testcats/triax_prolate_unequal_los.txt", 'r')
#triax_file = open("/Users/brandyn/astro/triax/testcats/triax_sphere_pos.txt", 'r')
#xx = np.load("./testcats/x.npy")
#yy = np.load("./testcats/y.npy")

#triax_file = open('/work/04144/blee/lonestar/triaxpy/triaxcats/triax_spherical_bf','r')
#xx = np.load("/work/04144/blee/lonestar/triaxpy/triaxcats/x.npy")
#yy = np.load("/work/04144/blee/lonestar/triaxpy/triaxcats/y.npy")

clusterName = 'A209.W-S-I+'

triax_file = open("/Users/brandyn/astro/CLASHtxtFiles/"+clusterName+".good.cc.txt",'r')
xx = np.load('/Users/brandyn/astro/CLASHtxtFiles/x.npy')
yy = np.load('/Users/brandyn/astro/CLASHtxtFiles/y.npy')

orient,orient_err,axis_ratio,axis_ratio_err = dataHandle.return_clash_orient(clusterName)

g1 = []
g2 = []
for line in triax_file.readlines():
	fields=line.split(',')
	g1.append(float(fields[0]))
	g2.append(float(fields[1]))
	zs = float(fields[5])

g1 = np.array(g1)
g2 = np.array(g2)

center,zl,r500 = dataHandle.return_cluster_attributes(clusterName+".good.cc.cat")

r200chain, cchain, achain, bchain, thetachain, phichain = NFWtriaxMCmodel.MCfit(nwalkers, lenburn, lenchain, xgal=xx, ygal=yy, eps1=g1, eps2=g2,orient=orient, orient_err=orient_err, axis_ratio=axis_ratio, axis_ratio_err=axis_ratio_err, zl=zl,zs=zs)

print r200chain[-1]
print cchain[-1]
print achain[-1]
print bchain[-1]
print thetachain[-1]
print phichain[-1]

print "avg r200: "+str(np.average(r200chain))
print "avg c: "+str(np.average(cchain))
print "avg a: "+str(np.average(achain))
print "avg b: "+str(np.average(bchain))
print "avg theta:"+str(np.average(thetachain))
print "avg phi: "+str(np.average(phichain))