import wx
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.signal import argrelextrema

def get_path(wildcard):
	app = wx.App(None)
	style = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST
	dialog = wx.FileDialog(None, 'Open', wildcard=wildcard, style=style)
	dialog.SetDirectory("D:\\LIONEL\\PROGRAMMING\\CODE\\CORE\\IMAGEJ\\OPT_OPTIMISE_TRACE\\")
	if dialog.ShowModal() == wx.ID_OK:
		path = dialog.GetPath()
	else:
		path = None
	dialog.Destroy()
	return path

file = get_path('*.txt')

# import file
X = []
Y = []
Z = []
f = open(file,'r')
for line in f:
	row = line.rstrip().split(',')
	if len(row)>1:
		X.append(float(row[0]))
		Y.append(float(row[1]))
		Z.append(float(row[2]))
f.close()

X = np.array(X)
Y = np.array(Y)
Z = np.array(Z)

# Determine Curvature
CURV = []
for i in range(len(X)-2):
	X0 = X[i+0]
	X1 = X[i+1]
	X2 = X[i+2]
	Y0 = Y[i+0]
	Y1 = Y[i+1]
	Y2 = Y[i+2]
	Z0 = Z[i+0]
	Z1 = Z[i+1]
	Z2 = Z[i+2]
	
	V1 = np.array([X1-X0, Y1-Y0, Y1-Y0])
	V2 = np.array([X2-X1, Y2-Y1, Y2-Y1])
	L1 = np.sqrt(np.sum(V1*V1))
	L2 = np.sqrt(np.sum(V2*V2))
	
	cross = np.cross(V1,V2) / (L1*L2)
	curvature = np.sqrt(np.sum(cross*cross)) / (L1+L2)
	CURV.append(curvature)

CURV = np.array(CURV)	
P = argrelextrema(CURV, np.greater, order = 3)
POS = np.array(P[0], int)

# Plot curve
fig = plt.figure(1)
plt.plot(CURV)
plt.plot(POS, CURV[POS], 'ro')

# plot 3D
fig = plt.figure(2)
ax = fig.add_subplot(111, projection='3d')
ax.plot(X, Y, Z, label='root')
ax.plot(X[POS+1], Y[POS+1], Z[POS+1], 'ro')
#ax.invert_zaxis()
#max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max() / 2.0
#mean_x = X.mean()
#mean_y = Y.mean()
#mean_z = Z.mean()
#sxy = 2.
#ax.set_xlim(mean_x - max_range/2./sxy, mean_x + max_range/2./sxy)
#ax.set_ylim(mean_y - max_range/2./sxy, mean_y + max_range/2./sxy)
#ax.set_zlim(Z.min() , Z.max() )
fig.tight_layout(pad=0.4, w_pad=0.0, h_pad=1.0)
	
plt.show()