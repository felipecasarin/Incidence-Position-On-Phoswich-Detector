# -*- coding: cp1252 -*-
from __future__ import division 
from math import *
import numpy as np
import matplotlib.pyplot as plt
from numpy import cross, eye, dot
import ExpPat as EX
import scipy as sp
from scipy.optimize import least_squares

#list of experimental data point positions
xtab=EX.xtab
ytab=EX.ytab
Nx=xtab.size
Ny=ytab.size

ExDataTab=EX.ExDataTab #data: ExDataTab=intable.reshape(Nx,Ny,2,4)

def recaldata(Tab,b,c):
	rec=np.empty((Nx,Ny,2,4))
	for ix in range(0,Nx):
		for iy in range(0,Ny):
			L=(Tab[ix][iy][0]+c[0])*b[0]
			rec[ix][iy][0][0:4]=L[0:4] #lines
			C=(Tab[ix][iy][1]+c[1])*b[1]
			rec[ix][iy][1][0:4]=C[0:4] #columns
	return rec
	

def normdata(Tab):
	S=np.sum(Tab,axis=3)
	D=np.zeros((Nx,Ny,2,4))
	for ix in range(0,Nx):
		for iy in range(0,Ny):
			D[ix][iy][0][0:4]=Tab[ix][iy][0][0:4]/S[ix][iy][0] #lines
			D[ix][iy][1][0:4]=Tab[ix][iy][1][0:4]/S[ix][iy][1] #columns
	return D

