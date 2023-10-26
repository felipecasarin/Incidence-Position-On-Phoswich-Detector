# -*- coding: cp1252 -*-
from __future__ import division 
from math import *
import numpy as np
import matplotlib.pyplot as plt
from numpy import cross, eye, dot
import ImageModelClass as IMC
import ExpPat as EX
import ExpReCal as ERC
import scipy as sp
from scipy.optimize import least_squares
import sys

print("test commit")


#list of experimental data point positions
xtab=EX.xtab
ytab=EX.ytab
Nx=xtab.size
Ny=ytab.size

ExDataTab=EX.ExDataTab #data: ExDataTab=intable.reshape(Nx,Ny,2,4)
RelErrorDat=EX.RelErrorDat

a=np.ones((4,4))
b=np.ones((2,4))
c=np.zeros((2,4))
ExN=ERC.normdata(ERC.recaldata(ExDataTab,b,c))
#print (ExN)

eps0=0.5
#Ybg0=1.E-4

def SaveParameters():
	np.savetxt("par_eps.txt",np.array([eps0]))
	#np.savetxt("par_Ybg0.txt",np.array([Ybg0]))
	np.savetxt("par_a.txt", np.ravel(a))
	np.savetxt("par_b.txt", np.ravel(b))
	np.savetxt("par_c.txt", np.ravel(c))
	return

IM=IMC.ImageModel(eps0)

YpixModel=IM.YpixModel #model: YpixPattern(x,y,ep=eps, a=np.ones((4,4)),istonormalize=1):
YPattern=IM.YPattern	#model: YPattern(Ypm,direction,a=np.ones((4,4)),b=np.ones((2,4)),istonormalize=1)


def Ytable(xoff=0.,yoff=0.,istonorm=1,a=np.ones((4,4))): # theoretical yields
	YT=np.zeros((Nx,Ny,2,4))
	for ix in range(0,Nx):
		for iy in range(0,Ny):
			x=xtab[ix]-xoff
			y=ytab[iy]-yoff
			Ypm=YpixModel(x,y)*a
			YT[ix][iy][0][0:4]=YPattern(Ypm,0,istonorm) #lines
			YT[ix][iy][1][0:4]=YPattern(Ypm,1,istonorm) #columns       
            
	return YT

def errorfunc(xoff,yoff,dTab): #dTab is the data table after calib.+ renorm. for example
	return Ytable(xoff,yoff)-dTab
	
def errorfunc(xoff,yoff,dTab,YTab): #dTab is the data table after calib.+ renorm. for example
	return YTab-dTab

def PatternError(x,y,LCdat): 
	Ypm=YpixModel(x,y) # theory
	YE=YPattern(Ypm,0)-LCdat[0]
	YE=np.append(YE,YPattern(Ypm,1)-LCdat[1])
	return YE

# fitting functions ****
def fabc(X,Y_dat,Y_sig):
	IM.eps=X[0] # sets model parameters for errorfunc eval
	IM.Ybg=X[1]
	#a=np.ones((4,4)) # not going to fit a
	b=np.reshape(X[2:10],(2,4))
	c=np.reshape(X[10:18],(2,4))
	a=np.reshape(X[18:34],(4,4))
	ExN=ERC.normdata(ERC.recaldata(Y_dat,b,c))# 
	YTab=Ytable(0,0,1,a)
	#Tab=np.reshape(ExN,(Nx,Ny,2,4))
	f=np.ravel(errorfunc(0.,0.,ExN,YTab)/abs(Y_sig))[0:328] # drop last experimental point (4l,4c, [328:336])
	return f
	
def fbc(X,Y_dat,Y_sig):
	IM.eps=X[0] # sets model parameters for errorfunc eval
	#IM.Ybg=X[1]
	#a=np.ones((4,4)) # not going to fit a
	b=np.reshape(X[1:9],(2,4))
	c=np.reshape(X[9:17],(2,4))
	ExN=ERC.normdata(ERC.recaldata(Y_dat,b,c))# 
	YTab=Ytable(0,0)
	#Tab=np.reshape(ExN,(Nx,Ny,2,4))
	f=np.ravel(errorfunc(0.,0.,ExN,YTab)/abs(Y_sig))[0:328] # drop last experimental point (4l,4c, [328:336])
	return f
	
def fb(X,Y_dat,Y_sig):
	IM.eps=X[0] # sets model parameters for errorfunc eval
	IM.Ybg=X[1]
	#a=np.ones((4,4)) # not going to fit a
	b=np.reshape(X[2:10],(2,4))
	ExN=ERC.normdata(ERC.recaldata(Y_dat,b,np.zeros((2,4))))# 
	YTab=Ytable(0,0)
	#Tab=np.reshape(ExN,(Nx,Ny,2,4))
	f=np.ravel(errorfunc(0.,0.,ExN,YTab)/abs(Y_sig))[0:328] # drop last experimental point (4l,4c, [328:336])
	return f



def fpat(X,Y_dat): # pattern error function
	x=X[0]
	y=X[1]
	f=PatternError(x,y,Y_dat)
	return f

def f0(X,Y_dat): # for fitting sum L - sum C =0
	b=np.reshape(X[0:8],(2,4))
	c=np.reshape(X[8:16],(2,4))
	ExD=(ERC.recaldata(Y_dat,b,c))
	F=np.array([])
	for ixtab in range(0,Nx):
		for iytab in range(0,Ny):
			if not ( (ixtab==6) and (iytab==5) ):
				L=np.sum(ExD[ixtab][iytab][0][0:4])
				C=np.sum(ExD[ixtab][iytab][1][0:4])
				F=np.append(F,(L-C))
	return F
				
#****
# fitting procedures ****


def fitpat(ixtab,iytab): # pattern fit of x,y
	# loads model and calib. fit parameters
	IM.eps=np.loadtxt("par_eps.txt")
	IM.Ybg=np.loadtxt("par_Ybg0.txt")
	a=np.loadtxt("par_a.txt")
	b=np.loadtxt("par_b.txt")
	c=np.loadtxt("par_c.txt")
	ExN=ERC.normdata(ERC.recaldata(ExDataTab,b,c))
	xguess=np.sum(IM.ruler*ExN[ixtab][iytab][0])*3.
	yguess=np.sum(IM.ruler*ExN[ixtab][iytab][1])*3.
	#print("x,y guess:",xguess,yguess)
	X=np.array([xguess,yguess]) # initial x,y guess
	Lmax=10.
	bdi=np.array([-Lmax,-Lmax]) # lower limits
	bds=np.array([Lmax,Lmax]) # upper limits	
	result = sp.optimize.least_squares(fpat,X,bounds=(bdi,bds),args = ([ExN[ixtab][iytab][0:2]]))
	#print("hole x,y=",xtab[ixtab],ytab[iytab], "\n x,y fit result:")
	print(result.x)
	#fil=open("xy_result.txt",'w')
	#print(result,file=fil)
	#fil.close()	
	return result.x

def fitabc(eps0=0.5,Ybg0=0.0,a=np.ones((4,4)),b=np.ones((2,4)),c=np.zeros((2,4)),blim=2.,clim=25.): #blim minimum, xylim min/max both x/y
	#initial values
	X=np.array([])
	X=np.append(X,eps0) # initial eps
	X=np.append(X,Ybg0) # initial Ybg
	X=np.append(X,np.ravel(b))  # initial b's
	X=np.append(X,np.ravel(c))  # initial c's
	X=np.append(X,np.ravel(a))  # initial a's
	#lower limits:
	bdi=np.array([])
	bdi=np.append(bdi,0.)
	bdi=np.append(bdi,0.) 
	bdi=np.append(bdi,np.ones((8))/blim) 
	bdi=np.append(bdi,np.ones((8))*(-clim)) 	
	bdi=np.append(bdi,np.ones((16))*0.1) 	#a min
	bdi[3]=0.999
	bdi[11]=-0.01 # effectively fixes b,c for this line
	#upper limits:
	bds=np.array([])
	bds=np.append(bds,1.)
	bds=np.append(bds,0.0001)
	bds=np.append(bds,np.ones((8))*blim) 
	bds=np.append(bds,np.ones((8))*(clim))
	bds=np.append(bds,np.ones((16)))
	bds[3]=1.001
	bds[11]=0.01

	result = sp.optimize.least_squares(fabc,X,bounds=(bdi,bds),args = ([ExDataTab,ExN*RelErrorDat]))
	print(result)
	np.savetxt("par_eps.txt",np.array([result.x[0]]))
	np.savetxt("par_Ybg0.txt",np.array([result.x[1]]))
	np.savetxt("par_b.txt", result.x[2:10])
	np.savetxt("par_c.txt", result.x[10:18])
	np.savetxt("par_a.txt", result.x[18:34])
	return result.x


def fitbc(eps0=0.5,Ybg0=0.,b=np.ones((2,4)),c=np.zeros((2,4)),blim=10.,clim=0.01): #blim minimum, xylim min/max both x/y
	#initial values
	X=np.array([])
	X=np.append(X,eps0) # initial eps
	X=np.append(X,Ybg0) # initial Ybg
	X=np.append(X,np.ravel(b))  # initial b's
	X=np.append(X,np.ravel(c))  # initial c's
	#lower limits:
	bdi=np.array([0.])
	bdi=np.append(bdi,0.) 
	bdi=np.append(bdi,np.ones((8))/blim) 
	bdi=np.append(bdi,np.ones((8))*(-clim)) 	
	bdi[3]=0.999
	bdi[11]=-0.01 # effectively fixes b,c for this line
	#upper limits:
	bds=np.array([])
	bds=np.append(bds,1.)
	bds=np.append(bds,0.0001)
	bds=np.append(bds,np.ones((8))*blim) 
	bds=np.append(bds,np.ones((8))*(clim))
	bds[3]=1.001
	bds[11]=0.01

	result = sp.optimize.least_squares(fbc,X,bounds=(bdi,bds),args = ([ExDataTab,ExN*RelErrorDat]))
	print(result)
	np.savetxt("par_eps.txt",np.array([result.x[0]]))
	np.savetxt("par_Ybg0.txt",np.array([result.x[1]]))
	np.savetxt("par_b.txt", result.x[2:10])
	np.savetxt("par_c.txt", result.x[10:18])
	return result.x

def newfitbc(eps0=0.5,Ybg0=0.005,b=np.ones((2,4)),c=np.zeros((2,4)),blim=10.,clim=1000.): #blim minimum, xylim min/max both x/y
#initial values
	print("Inicializando")
	sys.stdout.flush()
	c.fill(100.)
	X=np.array([])
	X=np.append(X,eps0) # initial eps
	#X=np.append(X,Ybg0) # initial Ybg
	X=np.append(X,np.ravel(b)) # initial b’s
	X=np.append(X,np.ravel(c)) # initial c’s
#lower limits:
	bdi=np.array([0.])
	#bdi=np.append(bdi,0.)
	bdi=np.append(bdi,np.ones((8))/blim)
	bdi=np.append(bdi,np.ones((8))*(-clim))
	bdi[3]=0.999
	bdi[11]=100. # effectively fixes b,c for this line
#upper limits:
	bds=np.array([])
	bds=np.append(bds,1.)
	#bds=np.append(bds,0.01)
	bds=np.append(bds,np.ones((8))*blim)
	bds=np.append(bds,np.ones((8))*(clim))
	bds[3]=1.001
	bds[11]=110.
#bds[1]=0.01
	result = sp.optimize.least_squares(fbc,X,bounds=(bdi,bds),args = ([ExDataTab,ExN*RelErrorDat]),loss='cauchy',tr_solver='lsmr', verbose=2)
##Configurar f_scale e explorar outras configurações##
##loss arctan?##
##Tentar methods diferentes##
	print(result)
	np.savetxt("par_eps_tryNOYBG.txt",np.array([result.x[0]]))
	#np.savetxt("par_Ybg0_try.txt",np.array([result.x[1]]))
	np.savetxt("par_b_tryNOYBG.txt", result.x[1:9])
	np.savetxt("par_c_tryNOYBG.txt", result.x[9:17])
	return result.x

#****     

def fitb(eps0=0.5,Ybg0=0.,b=np.ones((2,4)),blim=10.): #blim minimum, xylim min/max both x/y
	#initial values
	X=np.array([])
	X=np.append(X,eps0) # initial eps
	X=np.append(X,Ybg0) # initial Ybg
	X=np.append(X,np.ravel(b))  # initial b's
	#lower limits:
	bdi=np.array([0.])
	bdi=np.append(bdi,0.) 
	bdi=np.append(bdi,np.ones((8))/blim) 	
	bdi[3]=0.999
	#upper limits:
	bds=np.array([])
	bds=np.append(bds,1.)
	bds=np.append(bds,0.0001)
	bds=np.append(bds,np.ones((8))*blim) 
	bds[3]=1.001

	result = sp.optimize.least_squares(fb,X,bounds=(bdi,bds),args = ([ExDataTab,ExN*RelErrorDat]), loss = 'cauchy')
	print(result)
	np.savetxt("par_eps.txt",np.array([result.x[0]]))
	np.savetxt("par_Ybg0.txt",np.array([result.x[1]]))
	np.savetxt("par_b.txt", result.x[2:10])
	return result.x




#fitbc(0.5,0.0,np.ones((2,4)),np.zeros((2,4)),blim=2.,clim=0.1)


def fitzero(b=np.ones((2,4)),c=np.zeros((2,4))): #
	#initial values
	X=np.array([])
	X=np.append(X,np.ravel(b))  # initial b's
	X=np.append(X,np.ravel(c))  # initial c's
	#lower limits:
	bdi=np.array([])
	bdi=np.append(bdi,np.ones((8))*0.1) 
	bdi=np.append(bdi,np.ones((8))*(-50.))	

	#upper limits:
	bds=np.array([])
	bds=np.append(bds,np.ones((8))*10.)
	bds=np.append(bds,np.ones((8))*50.) 	
	
	#effectively fix b[1] and c[1]
	bdi[1]=0.999999
	bds[1]=1.000001
	bdi[8]=-0.000001
	bds[8]=0.000001
	F=f0(X,ExDataTab)
	print()
	print(F)
	print("rms-i:",sqrt(np.sum(np.power(F,2)))/F.size)
	result = sp.optimize.least_squares(f0,X,bounds=(bdi,bds),args = ([ExDataTab]))
	print(result)
	np.savetxt("par_b.txt", result.x[0:8])
	np.savetxt("par_c.txt", result.x[8:16])
	F=f0(result.x,ExDataTab)	
	print()
	print(F)
	print("rms-f:",sqrt(np.sum(np.power(F,2)))/F.size)

	return result.x


#fitzero()
#fitbc(0.5,0.,b=np.ones((2,4)),c=np.zeros((2,4)),blim=2.,clim=25.)
#fitabc()
