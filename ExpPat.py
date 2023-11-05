# -*- coding: cp1252 -*-
from __future__ import division 
from math import *
import numpy as np
import matplotlib.pyplot as plt
from numpy import cross, eye, dot
from mpl_toolkits.mplot3d import Axes3D

tabname="./Mean.txt"

intable=np.loadtxt(tabname)


xtab=np.array([-5.3,-3.3,-0.3,1.7,4.7])
ytab=np.array([-5.3,-3.3,-0.3,1.7,4.7,6.7])






Nx=xtab.size
Ny=ytab.size

ExDataTab=intable.reshape(Nx,Ny,2,4)
def SwapLineData(Tab):
	for ix in range(0,Nx):
		for iy in range(0,Ny):
			L=Tab[ix][iy][0][0:4]
			L=np.flip(L,0)
			Tab[ix][iy][0][0:4]=L

SwapLineData(ExDataTab) # lines swapped (L ... L4 -> L4 ... L1) to conform to model line order 

errortabname="./Error.txt"
errortab=np.loadtxt(errortabname)
ErrorDat=errortab.reshape(Nx,Ny,2,4)
SwapLineData(ErrorDat) # lines swap ...

RelErrorDat=(ErrorDat/ExDataTab)

lcal=np.array([1.,1.,1.,1.])
ccal=np.array([1.,1.,1.,1.])

cal=np.array([lcal,ccal])


def SetUnCal():
	lcal=np.ones(4)
	ccal=np.ones(4)
	cal=np.array([lcal,ccal])
	
def SetCal(al,ac):
	lcal=al
	ccal=ac
	cal=np.array([lcal,ccal])

def FunX(iy,il,ilc):
	a=np.array([])
	for ix in range(0,7):
		a=np.append(a,ExDataTab[ix][iy][ilc][il]/cal[ilc][il])
	return a

def FunY(ix,il,ilc):
	a=np.array([])
	for iy in range(0,6):
		a=np.append(a,ExDataTab[ix][iy][ilc][il]/cal[ilc][il])
	return a
	
def FunXN(iy,il,ilc):
	a=np.array([])
	for ix in range(0,7):
		norm=np.sum(ExDataTab[ix][iy][ilc][0:4]/cal[ilc][0:4])
		a=np.append(a,ExDataTab[ix][iy][ilc][il]/cal[ilc][il]/norm)
	return a
	
def FunYN(ix,il,ilc):
	a=np.array([])
	for iy in range(0,6):
		norm=np.sum(ExDataTab[ix][iy][ilc][0:4]/cal[ilc][0:4])
		a=np.append(a,ExDataTab[ix][iy][ilc][il]/cal[ilc][il]/norm)
	return a
	

	
def exp_xy(ixtab,iytab):
	return np.array([xtab[ixtab],ytab[iytab]])

def patL(ixtab,iytab): # returns line pattern of exp. point ixtab,iytab
	tablin=ixtab
	tabcol_i=iytab*8
	tabcol_f=tabcol_i+4
	return (intable[tablin][tabcol_i:tabcol_f])/lcal

def patC(ixtab,iytab):
	tablin=ixtab
	tabcol_i=iytab*8+4
	tabcol_f=tabcol_i+4
	return (intable[tablin][tabcol_i:tabcol_f])/ccal
	
def sumL(ixtab,iytab):
	return np.sum(patL(ixtab,iytab))
	
def sumC(ixtab,iytab):
	return np.sum(patC(ixtab,iytab))

def NpatL(ixtab,iytab):
	return patL(ixtab,iytab)/sumL(ixtab,iytab)

def NpatC(ixtab,iytab):
	return patC(ixtab,iytab)/sumC(ixtab,iytab)

def DataSetSumL():
	s=0.
	for ixtab in range(0,Nx):
		for iytab in range(0,Ny):
			s=s+sumL(ixtab,iytab)
	return s  #39160.82
	
def DataSetSumC():
	s=0.
	for ixtab in range(0,Nx):
		for iytab in range(0,Ny):
			s=s+sumC(ixtab,iytab)
	return s #42251.22



