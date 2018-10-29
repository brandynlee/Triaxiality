import dataHandle
import numpy as np

fileNames=['A209.W-S-I+.good.cc.cat', 'MACS0429-02.W-C-RC.good.cc.cat', 'MACS0647+70.W-C-RC.good.cc.cat',
	'MACS0717+37.W-J-V.good.cc.cat', 'MACS0744+39.W-C-RC.good.cc.cat',
	'MACS1115+01.W-C-RC.good.cc.cat', 'MACS1149+22.W-J-V.good.cc.cat', 'MACS1206-08.W-C-RC.good.cc.cat', 
	'MACS1347-11.W-C-RC.good.cc.cat', 'MACS1423+24.W-C-IC.good.cc.cat', 'MACS1532+30.W-C-RC.good.cc.cat',
	'MACS1720+35.W-J-V.good.cc.cat', 'MACS2129-07.W-C-RC.good.cc.cat', 'RXJ2129.W-J-V.good.cc.cat']
outPath = '/Users/brandyn/astro/CLASHtxtFiles/'

rmax = 3.
npix = 100

for fileName in fileNames:
	print fileName[:-4]

	x,y,eps1,eps2,zS = dataHandle.process_wtg(fileName)

	xgrid = np.linspace(-3.,3.,npix)
	ygrid = np.linspace(-3.,3.,npix)

	xxgrid, yygrid = np.meshgrid(xgrid,ygrid)

	x_list = xgrid.tolist()
	y_list = ygrid.tolist()

	eps1_avg = []
	eps2_avg = []
	eps1_std = []
	eps2_std = []
	ngal = []
	x_out = []
	y_out = []
	zs_list = []
	Ngal = 0
	for i in range(len(x_list)-1):
		keep_x = np.logical_and(x_list[i]<x,x<x_list[i+1])
		x_temp_0 = x[keep_x]
		y_temp_0 = y[keep_x]
		eps1_temp_0 = eps1[keep_x]
		eps2_temp_0 = eps2[keep_x]
		for j in range(len(y_list)-1):
			keep_y = np.logical_and(y_list[j]<y_temp_0,y_temp_0<y_list[j+1])
			x_temp = x_temp_0[keep_y]
			y_temp = y_temp_0[keep_y]
			eps1_temp = eps1_temp_0[keep_y]
			eps2_temp = eps2_temp_0[keep_y]

			x_out.append(x_list[i])
			y_out.append(y_list[j])
			eps1_avg.append(np.average(eps1_temp)) 
			eps2_avg.append(np.average(eps2_temp))
			eps1_std.append(np.std(eps1_temp))
			eps2_std.append(np.std(eps2_temp))
			ngal.append(eps1_temp.shape[0])
			Ngal = Ngal+ngal[-1]
			zs_list.append(zS)

	print Ngal
	np.save(outPath+"x.npy",np.array(x_out))
	np.save(outPath+"y.npy",np.array(y_out))
	#print zs_list
	np.savetxt(outPath+fileName[:-4]+".txt", np.transpose([np.array(eps1_avg), np.array(eps2_avg), np.array(eps1_std), np.array(eps2_std), np.array(ngal), np.array(zs_list)]),delimiter=',', fmt='%s')


