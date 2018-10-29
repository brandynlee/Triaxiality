import numpy as np
import cPickle
from astropy.coordinates import SkyCoord
from astropy import units as u
import sys
#sys.path.append('/Users/brandyn/astro/baryonimpact')
import colossus.cosmology.cosmology as Cosmology
#import pzutil as util
from astropy.io import fits

def returnTrueInfo():

	catPath = '/Users/brandyn/astro/cosmoOWLS/SIMS/matched_cluster_catalog_z_0.250.dat'
	cat = open(catPath, 'r')
	
	fof_dm = []
	fof_agn80 = []
	fof_agn87 = []
	lm_dm = []
	lm_agn80 = []
	lm_agn87 = []
	r_dm = []
	r_agn80 = []
	r_agn87 = []
	for line in cat.readlines():
		fields = line.split()
		fof_dm.append(int(fields[0]))
		lm_dm.append(float(fields[1]))
		r_dm.append(float(fields[2]))
		fof_agn80.append(int(fields[3]))
		lm_agn80.append(float(fields[4]))
		r_agn80.append(float(fields[5]))
		fof_agn87.append(int(fields[6]))
		lm_agn87.append(float(fields[7]))
		r_agn87.append(float(fields[8]))

	info_dm = [fof_dm, lm_dm, r_dm]
	info_agn80 = [fof_agn80, lm_agn80, r_agn80]
	info_agn87 = [fof_agn87, lm_agn87, r_agn87]
	return (info_dm, info_agn80, info_agn87)

def returnTrueInfo_wref(logM=True, mdef='m200c', missing=False):

	if mdef=='m200c':
		catPath = '/Users/brandyn/astro/cosmoOWLS/SIMS/matched_cluster_catalog_z_0.250_wref.dat'
	elif mdef=='m500c':
		catPath = '/Users/brandyn/astro/cosmoOWLS/SIMS/matched_cluster_catalog_z_0.250_m500_wref.dat'
	elif mdef=='m200m':
		catPath = '/Users/brandyn/astro/cosmoOWLS/SIMS/matched_cluster_catalog_z_0.250_m200m_wref.dat'
	else:
		raise Exception("Invalid mass definition")

	cat = open(catPath, 'r')
	
	fof_dm = []
	fof_agn80 = []
	fof_agn87 = []
	fof_ref = []
	m_dm = []
	m_agn80 = []
	m_agn87 = []
	m_ref = []
	r_dm = []
	r_agn80 = []
	r_agn87 = []
	r_ref = []
	missingFlag = []
	for line in cat.readlines():
		fields = line.split()
		fof_dm.append(int(fields[0]))
		r_dm.append(float(fields[2]))
		fof_agn80.append(int(fields[3]))
		r_agn80.append(float(fields[5]))
		fof_agn87.append(int(fields[6]))
		r_agn87.append(float(fields[8]))
		fof_ref.append(int(fields[9]))
		r_ref.append(float(fields[11]))
		if missing:
			missingFlag.append(int(fields[12]))

		if logM:
			m_dm.append(float(fields[1]))
			m_agn80.append(float(fields[4]))
			m_agn87.append(float(fields[7]))
			m_ref.append(float(fields[10]))
		else:
			m_dm.append(10.**float(fields[1]))
			m_agn80.append(10.**float(fields[4]))
			m_agn87.append(10.**float(fields[7]))
			m_ref.append(10.**float(fields[10]))

	info_dm = [fof_dm, m_dm, r_dm]
	info_agn80 = [fof_agn80, m_agn80, r_agn80]
	info_agn87 = [fof_agn87, m_agn87, r_agn87]
	info_ref = [fof_ref, m_ref, r_ref]
	if missing:
		return (info_dm, info_agn80, info_agn87, info_ref, missingFlag)
	return (info_dm, info_agn80, info_agn87, info_ref)

def return_cluster_centers(): #for use on lonestar 5 with AGN 8.0 sims, current version is in generate photonlist
	catPath = '/work/04144/blee/lonestar/cosmoOWLS/SIMS/AGN_L400N1024_WMAP7/z_0.250/xray/cluster_catalog_AGN80_for_xray.dat'
	cat = open(catPath,'r')

	fof = []
	posgrp = []
	r200 = []
	r500 = []

	for line in cat.readlines():
		fields = line.split()
		fof.append(int(fields[0]))
		posgrp.append((float(fields[1]),float(fields[2]),float(fields[3])))
		r200.append(float(fields[4]))
		r500.append(float(fields[5]))
	return (fof,posgrp,r200,r500)

#input is a list of m200s, arrays of bins that are not log10'd beforehand
def bin_masses(masses, bins, logbin = False, lists=False):
	if lists:
		masses = np.array(masses)

	if logbin:
		masses = np.log10(masses)
		#bins = np.log10(bins)

	return np.digitize(masses, bins)

def cluster_like_inputs(mcmc_path, sim_suffix, fof_bin_list, bin_m200_true, num_samples=2000, c_cap=None, projection='all', binnum = None):
	nclusters_bin = len(fof_bin_list)
	projections = ['x','y','z']

		#m_sample = np.zeros((nbinclusters*3,nsamples))

	m_sample = []
	delta_lnm = []
	clusters_excluded = 0	
	for i in range(nclusters_bin):
		fof = fof_bin_list[i]
		m200_true = bin_m200_true[i]

		#if binnum==0:
			#samplefile = mcmc_path+'cluster_'+str(fof)+'_'+sim_suffix+'.pkl'
		#elif binnum == 2 or binnum == 3:
		samplefile = mcmc_path+'cluster_'+str(fof)+'_'+sim_suffix+'_m500.pkl'
		with open(samplefile,'rb') as input:
			cluster_chains = cPickle.load(input)
			for proj in projections:
				if projection <> 'all' and proj <> projection: continue

				m200chain = cluster_chains[proj]['m500']
				cchain = cluster_chains[proj]['c500']

				if c_cap is not None:
					m200chain = m200chain[cchain<c_cap]
					cchain = cchain[cchain<c_cap]
				if len(m200chain)<num_samples:
					print fof
					print proj
					clusters_excluded = clusters_excluded + 1
					continue
				idx = np.random.choice(np.arange(len(m200chain)),num_samples,replace=False)
				m200chain = m200chain[idx]
				cchain = cchain[idx]
				delta_lnm.append(np.log(m200chain) - np.log(m200_true))
				m_sample.append(m200chain)

	print("num of clusters excluded: " + str(clusters_excluded))
	m_sample = np.array(m_sample)
	delta_lnm = np.array(delta_lnm)

	return (m_sample, delta_lnm)

def cluster_like_inputs_analyticNFW(mcmc_path,id_list, noise_suffix = '25', c_cap = None, m_cap=None, num_samples=2000):
	nclusters_bin = len(id_list)
	m_sample = []
	delta_lnm = []

	for i in range(nclusters_bin):
		cid = id_list[i]
		samplefile = mcmc_path+'nfwcluster_'+str(cid)+'_noise'+noise_suffix+'.pkl'
		with open(samplefile,'rb') as input:
			cluster_chains = cPickle.load(input)

			m200chain = cluster_chains['m200']
			cchain = cluster_chains['c200']
			m200_true = cluster_chains['m200_true']
			if m_cap is not None:
				if m200_true >= m_cap: continue

			if c_cap is not None:
				m200chain = m200chain[cchain<c_cap]
				cchain = cchain[cchain<c_cap]

			if len(m200chain)<num_samples:
				print cid
				continue

			idx = np.random.choice(np.arange(len(m200chain)),num_samples,replace=False)
			m200chain = m200chain[idx]
			cchain = cchain[idx]

			delta_lnm.append(np.log(m200chain) - np.log(m200_true))
			m_sample.append(m200chain)

	m_sample = np.array(m_sample)
	delta_lnm = np.array(delta_lnm)

	return (m_sample, delta_lnm)

def read_bias_scatter(bias_file, use_bins):
	with open(bias_file, 'rb') as input:
		chains_in = cPickle.load(input)
		mumean = []
		sigmean = []
		for bin_idx in use_bins:

			sigchain = chains_in['sigma'][bin_idx]
			muchain = chains_in['mu'][bin_idx]

			mumean.append(np.percentile(muchain,50))
			sigmean.append(np.percentile(sigchain,50))

	return (mumean, sigmean)

def return_wtg_catalog(path, fileName):
	cat = fits.open(path+fileName)

	RA = cat[1].data['ALPHA_J2000']
	Dec = cat[1].data['DELTA_J2000']
	eps1 = cat[1].data['gs1corr']  # gs1
	eps2 = cat[1].data['gs2corr']  # gs2
	# Extract redshift information
	zS = cat[2].data['z_s']
	mean_zS = np.mean(zS)

	return (RA, Dec, eps1, eps2, mean_zS)

def return_cluster_attributes(fileName):
	if fileName == 'A209.W-S-I+.good.cc.cat':
    	#center[0] in u.hourangle, center[1] in u.degrees
  		center = ('01:31:53.139','-13:36:48.35')
		zL = 0.206
		r500 = 1.53
	#if fileName == 'A383.W-S-I+.good.cc_t1.cat':
    	#center = ('02:48:03.268', '-03:31:46.43')
    	#zL = 0.188
	#if fileName == 'A611.W-C-IC.good.cc_t1.cat':
    	#center = ('08:00:56.818', '36:03:25.52')
    	#zL = 0.288
	#if fileName == 'A2261.W-C-RC.good.cc_t1.cat':
    	#center = ('17:22:26.986', '32:07:57.89')
    	#zL = 0.224
	#if fileName == 'MACS0329-02.W-J-V.good.cc_t1.cat':
    	#center = ('03:29:41.459','-02:11:45.52')
    	#zL = 0.450
	if fileName == 'MACS0429-02.W-C-RC.good.cc.cat':
		center = ('04:29:36.001','-02:53:05.63')
		zL = 0.399
		r500=1.10
	if fileName== 'MACS0647+70.W-C-RC.good.cc.cat':
		center = ('06:47:49.682','70:14:56.05 ')
		zL = 0.592
		r500=1.26
	if fileName == 'MACS0717+37.W-J-V.good.cc.cat':
		center = ('07:17:32.088','37:45:20.94')
		zL = 0.546
		r500=1.69
	if fileName == 'MACS0744+39.W-C-RC.good.cc.cat':
		center = ('07:44:52.310','39:27:26.80')
		zL = 0.698
		r500=1.69  
	if fileName == 'MACS1115+01.W-C-RC.good.cc.cat':
		center = ('11:15:51.881 ','01:29:54.98')
		zL = 0.355
		r500 = 1.28
	if fileName == 'MACS1149+22.W-J-V.good.cc.cat':
		center = ('11:49:35.426','22:24:03.62')
		zL = 0.544 
		r500=1.53 
	if fileName == 'MACS1206-08.W-C-RC.good.cc.cat':
		center = ('12:06:12.293','-08:48:06.22')
		zL = 0.439
		r500=1.61  
	if fileName == 'MACS1347-11.W-C-RC.good.cc.cat':
    	#ID for this one on paper is 'RXJ1347.5-1144'
		center = ('13:47:30.778 ','-11:45:09.43')
		zL = 0.451
		r500=1.67 
	if fileName == 'MACS1423+24.W-C-IC.good.cc.cat':
		center = ('14:23:47.923','24:04:42.77')
		zL = 0.543
		r500=1.35  
	if fileName == 'MACS1532+30.W-C-RC.good.cc.cat':
    #ID on paper 'RXJ1532.8+3021'
		center = ('15:32:53.830 ','30:20:59.38')
		zL =  0.363
		r500=1.31
	if fileName == 'MACS1720+35.W-J-V.good.cc.cat':
		center = ('17:20:16.666','35:36:23.35')
		zL =  0.387 
		r500=1.14 
	if fileName == 'MACS2129-07.W-C-RC.good.cc.cat':
		center = ('21:29:25.723','-07:41:30.84')
		zL =  0.588
		r500=1.28 
	if fileName == 'RXJ2129.W-J-V.good.cc.cat':
		center = ('21:29:39.727','00:05:18.15')
		zL =  0.235  
		r500=1.28
	return (center,zL, r500)

def process_wtg(fileName):
	
	params = {'flat': True, 'H0': 70.0, 'Om0': 0.3, 'Ob0': 0.049, 'sigma8': 0.81, 'ns': 0.95}
	Cosmology.addCosmology('myCosmo', params)
	cosmo = Cosmology.setCosmology('myCosmo')

	#chooseCosmology = 'WMAP7-ML'
	#cosmo=Cosmology.setCosmology(chooseCosmology)
	dataPath = '/Users/brandyn/astro/clash_cat/'
	center, zL, r500 = return_cluster_attributes(fileName)
	Dl = cosmo.angularDiameterDistance(zL)/cosmo.h #Dl = Mpc from Mpc/h

	RA, Dec, eps1, eps2, zS = return_wtg_catalog(dataPath,fileName)
	c = SkyCoord(ra=RA,dec=Dec, distance=Dl, unit=(u.degree,u.degree,u.Mpc),frame='icrs')
	xcen = SkyCoord(center[0], center[1],distance = Dl, unit = (u.hourangle,u.degree, u.Mpc),frame='icrs')

	aframe=xcen.skyoffset_frame() #create new frame centered at cluster center
	coord = c.transform_to(aframe) #transform coordinates to new cluster frame

	x= np.array(-coord.cartesian.y) #negative cuz RA points opposite way
	y = np.array(coord.cartesian.z)

	f500=.073

	psi = np.arctan2(y,x)
	rgal = (x**2 + y**2)**0.5
	epstan = -(eps1*np.cos(2*psi)+eps2*np.sin(2*psi))
	epstancoor = epstan/(1-f500*np.exp(1-(rgal/r500)))
	eps1 = eps1/(1-f500*np.exp(1-rgal/r500))
	eps2 = eps2/(1-f500*np.exp(1-rgal/r500))

	#return (x,y, epstancorr, zS)
	return(x,y,eps1,eps2, zS)

def return_clash_orient(fileName):
	orient_file = open('/Users/brandyn/astro/CLASHtxtFiles/cluster_orientations.txt','r')

	for line in orient_file.readlines()[1:]:
		fields = line.split()
		if fileName == fields[0]:
			orient = float(fields[1])
			orient_err = float(fields[2])
			if float(fields[3])>1:
				axis_ratio = 1./float(fields[3])
			else:
				axis_ratio = float(fields[3])
			axis_ratio_err = float(fields[4])

	return (orient,orient_err,axis_ratio,axis_ratio_err)

#return rs fits from color cut method from overlapping WTG clusters
def return_wtg_fits(clusterName):
	#cosmo = Cosmology.setCosmology('WMAP7-ML')
	params = {'flat': True, 'H0': 70.0, 'Om0': 0.3, 'Ob0': 0.049, 'sigma8': 0.81, 'ns': 0.95}
	Cosmology.addCosmology('myCosmo', params)
	cosmo = Cosmology.setCosmology('myCosmo')
	center, zl, zs = return_redshift_center(clusterName+'.txt')

	if clusterName == 'A383.W-S-I+.good.cc_t1':
		rs = 0.45
 
	if clusterName == 'A209.W-S-I+.good.cc_t1':
		rs = 0.56
 
	if clusterName == 'A611.W-C-IC.good.cc_t1':
		rs = 0.50
 
	if clusterName == 'A2261.W-C-RC.good.cc_t1':
		rs = 0.63
 
	if clusterName == 'MACS0329-02.W-J-V.good.cc_t1':
		rs = 0.52
 
	if clusterName == 'MACS0429-02.W-C-RC.good.cc_t1':
		rs = 0.51
 
	if clusterName == 'MACS0647+70.W-C-RC.good.cc_t1':
		rs = 0.52
 
	if clusterName == 'MACS0717+37.W-J-V.good.cc_t1':
		rs = 0.66
 
	if clusterName == 'MACS0744+39.W-C-RC.good.cc_t1':
		rs = 0.56
 
	if clusterName == 'MACS1115+01.W-C-RC.good.cc_t1':
		rs = 0.51
 
	if clusterName == 'MACS1149+22.W-J-V.good.cc_t1':
		rs = 0.51
 
	if clusterName == 'MACS1206-08.W-C-RC.good.cc_t1':
		rs = 0.5
 
	if clusterName == 'MACS1347-11.W-C-RC.good.cc_t1':
		rs = 0.53
 
	if clusterName == 'MACS1423+24.W-C-IC.good.cc_t1':
		rs = 0.42
 
	if clusterName == 'MACS1532+30.W-C-RC.good.cc_t1':
		rs = 0.48
 
	if clusterName == 'MACS1720+35.W-J-V.good.cc_t1':
		rs = 0.4
 
	if clusterName == 'MACS2129-07.W-C-RC.good.cc_t1':
		rs = 0.52
 
	if clusterName == 'RXJ2129.W-J-V.good.cc_t1':
		rs = 0.39

	r200 = rs*4. #using concentration =4
	M200 = 200.*(4./3.)*np.pi*(r200)**3.*cosmo.rho_c(zl)*1000**3*cosmo.h**2

	return M200

def return_miyoung_fits():
	fits_dir = '/Users/brandyn/astro/CLASH/fits/'

	fits_file = open(fits_dir+'miyoung_c4_fits.txt', 'r')

	params = {'flat': True, 'H0': 70.0, 'Om0': 0.3, 'Ob0': 0.049, 'sigma8': 0.81, 'ns': 0.95}
	Cosmology.addCosmology('myCosmo', params)

	cluster_name = []
	M200 = []
	for line in fits_file.readlines():
		#cosmo = Cosmology.setCosmology('WMAP7-ML')
		cosmo = Cosmology.setCosmology('myCosmo')
		fields = line.split()
		cluster_name.append(fields[0].strip(',')[:-4])
		center, zl, zs = return_redshift_center(cluster_name[-1]+'.txt')
		r200 = float(fields[2])
		M200.append(200.*(4./3.)*np.pi*(r200)**3.*cosmo.rho_c(zl)*1000**3*cosmo.h**2)

	return (cluster_name,M200)

def return_matt_fits():
	fits_dir = '/Users/brandyn/astro/CLASH/fits/'
	fits_file = open(fits_dir+'matt_c4_fits.txt','r')

	M200 = []
	for line in fits_file.readlines():
		fields = line.split(',')
		M200.append(float(fields[0]))

	return M200