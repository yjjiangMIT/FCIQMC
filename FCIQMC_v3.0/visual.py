import matplotlib
matplotlib.rcParams['backend'] = 'Qt4Agg'
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

def draw_H():
	
	import integrals
	import ctrl_panel
	
	dim = len(ctrl_panel.key_list)
	H = []
	for key_i in ctrl_panel.key_list:
		vec = []
		for key_j in ctrl_panel.key_list:
			if key_i != key_j:
				entry = integrals.mat_element(key_i, key_j, integrals.ref_energy)
				vec.append(abs(entry))
			else:
				vec.append(0)
		vec.append(0)
		H.append(vec)
	H.append([0] * (dim+1))
	plot_once_2D(dim, H)

def new_spawn_map(dim):
	spawn_map = []
	for i in range(dim+1):
		spawn_map.append([1e-1] * (dim+1))
	return spawn_map

def plot_update_1D(fig_index, data, color, ax_range, label):
	plt.figure(fig_index)
	plt.clf()
	plt.axis(ax_range)
	x = data[0]
	for i in range(1,len(data)):
		plt.plot(x, data[i], color[i-1])
	plt.xlabel(label[0])
	plt.ylabel(label[1])
	plt.pause(0.01)

def plot_update_2D(dim, data):
	
	side = range(dim+1)
	X,Y = np.meshgrid(side, side)
	Z = np.array(data)
	
	plt.figure(4)
	plt.clf()
	plt.pcolormesh(X, Y, Z, norm = LogNorm())
	plt.colorbar()
	plt.axis('equal')
	plt.xlim([0, dim])
	plt.ylim([dim, 0])
	plt.pause(0.01)

def plot_once_2D(dim, data):
	
	side = range(dim+1)
	X,Y = np.meshgrid(side, side)
	Z = np.array(data)
	
	f = plt.figure(5)
	plt.ion()
	plt.pcolormesh(X, Y, Z, norm = LogNorm())
	plt.colorbar()
	plt.axis('equal')
	plt.xlim([0, dim])
	plt.ylim([dim, 0])
	plt.show()