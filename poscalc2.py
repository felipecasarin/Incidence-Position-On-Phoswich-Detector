from __future__ import division 
import ROOT
from ROOT import TNtuple
import sys
import numpy as np
from math import *
import matplotlib.pyplot as plt
from numpy import cross, eye, dot
import ImageModelClass as IMC
import ExpPat as EX
import ExpReCal as ERC
import scipy as sp
import Fits as ft
from scipy.optimize import least_squares
from ROOT import TFile, TTree, gROOT, addressof
from array import array
from scipy import linalg



#eps0=1.
#Ybg0=1.E-4


tab_eps0 = "./par_eps_tryNOYBG.txt"
eps0 = np.loadtxt(tab_eps0)

#tab_Ybg0 = "./par_Ybg0_2.txt"
#Ybg0 = np.loadtxt(tab_Ybg0)

#Import of calibration parameters b calculated by fitb() from Fits.py
tab_b="./par_b_tryNOYBG.txt"

b=np.loadtxt(tab_b)
#b=np.array([1,1,1,1,1,1,1,1])

tab_c="./par_c_tryNOYBG.txt"

c=np.loadtxt(tab_c)

#c = np.array([0,0,0,0,0,0,0,0])

IM=IMC.ImageModel(eps0)

Lpitch = 4.2

ruler=Lpitch*np.array([-1.5,-0.5,0.5,1.5])


####################################################################

#Opens the root file with the data to be processed
data = ROOT.TFile.Open("/home/casarin/Desktop/TCC/gp_23_filtered50000.root", 'read')


#Gets the ntuple from the given file
tree0 = data.Get('tree1')


#Creates the file in which the processed data will be saved afterwards

outfile =  ROOT.TFile("/home/casarin/Desktop/TCC/icpython-23-04-2022/gp_23_filtered50000_pos_2.root",'recreate')


print('Programa rodando')

#Calculates the Pattern Error Function using xguess, yguess and the Yield Data
def fpatxy(X,Y_dat): 
	x=X[0]
	y=X[1]
	f=ft.PatternError(x,y,Y_dat)   #<====PatternErrorOffset??
	return f



def Ytable(xhole,yhole,xoff=0.,yoff=0.,istonorm=1,a=np.ones((4,4))): # theoretical yields calculated using the hole position (xhole,yhole)
	YT=np.zeros((2,4))  
	x= xhole-xoff
	y= yhole-yoff
	print(x,y)        
	Ypm=IM.YpixModel(x,y)*a    
	YT[0][0:4]=IM.YPattern(Ypm,0,istonorm) #lines
	YT[1][0:4]=IM.YPattern(Ypm,1,istonorm) #columns  
	#print('break') 
	print(YT)       
	return YT

def fitpatxy(): # pattern fit of x,y
	# loads model and calib. fit parameters
    IM.eps=np.loadtxt("par_eps_tryNOYBG.txt")
    #IM.Ybg=np.loadtxt("par_Ybg0_2.txt")
    xguess=np.sum(ruler*newExN[0][0:4])*3.
    #print('old xguess=', xguess)
    
    #This part changes the value of xguess and yguess in case they lie outside the detector, as this would lead to miscalculations of the x and y coordinate
    if xguess < -10.:
        xguess = -8.
    if xguess > 10.:
        xguess = 8.
        
    #print('new xguess=', xguess)
    
    yguess=np.sum(ruler*newExN[1][0:4])*3.
    #print('old yguess=', yguess)
    
    if yguess < -10.:
        yguess=-8.
    if yguess > 10.:
        yguess = 8.
    
    #print('new yguess=', yguess)
    
	#print("x,y guess:",xguess,yguess)
    X=np.array([xguess,yguess]) # initial x,y guess
    Lmax=25. #Determines the maximum value for x and y
    bdi=np.array([-Lmax,-Lmax]) # lower limits
    bds=np.array([Lmax,Lmax]) # upper limits	
    
    global result
    #Calculates the x and y coordinates of the interaction point of the alpha particle on the phoswich detector
    result = sp.optimize.least_squares(fpatxy,X,bounds=(bdi,bds),args = ([newExN[0:2]]), ftol=1e-08, xtol=1e-08, gtol=1e-08, loss='soft_l1', tr_solver='lsmr')
    #result = sp.optimize.leastsq(fpatxy,X,bounds=(bdi,bds),args = ([newExN[0:2]]), ftol=1e-08, xtol=1e-08, gtol=1e-08, loss='cauchy', tr_solver='lsmr')
    #print(result)
   
    return result.x

def calib(l1,l2,l3,l4,c1,c2,c3,c4):
    
   
    #Calibrated lines
    calL1 = (l1 + c[0])*b[0]
    calL2 = (l2 + c[1])*b[1]
    calL3 = (l3 + c[2])*b[2]
    calL4 = (l4 + c[3])*b[3]
    sumL = calL1 + calL2 + calL3 + calL4
    
    
    #Calibrated columns
    
    calC1 = (c1 + c[4])*b[4]
    calC2 = (c2 + c[5])*b[5]
    calC3 = (c3 + c[6])*b[6]
    calC4 = (c4 + c[7])*b[7]
    sumC = calC1 + calC2 + calC3 + calC4
    
    
    #Normalized and calibrated lines
    cL1 = calL1/(sumL)
    cL2 = calL2/(sumL)
    cL3 = calL3/(sumL)
    cL4 = calL4/(sumL)
    #print('sumcL =',cL1+cL2+cL3+cL4)
    
    #Normalized and calibrated columns
    cC1 = calC1/(sumC)
    cC2 = calC2/(sumC)
    cC3 = calC3/(sumC)
    cC4 = calC4/(sumC)
   # print('sumcC =',cC1+cC2+cC3+cC4)
    

    global calb
    calb = np.array([cL1,cL2,cL3,cL4,cC1,cC2,cC3,cC4])
   # print(l1,l2,l3,l4,c1,c2,c3,c4)
    #print(calb)

    return calb


def pos0():
    #Creates the tree1 in which the new data will be saved
    tree1  = ROOT.TTree( 'tree1','tree1')
    
    
    nu = tree0.GetEntries()
    
    #Defines the arrays, filled with float type variables, with an initial value of zero
    aL1 = array('f',[0])
    aL2 = array('f',[0])
    aL3 = array('f',[0])
    aL4 = array('f',[0])
    aC1 = array('f',[0])
    aC2 = array('f',[0])
    aC3 = array('f',[0])
    aC4 = array('f',[0])
    aL1c = array('f',[0])
    aL2c = array('f',[0])
    aL3c = array('f',[0])
    aL4c = array('f',[0])
    aC1c = array('f',[0])
    aC2c = array('f',[0])
    aC3c = array('f',[0])
    aC4c = array('f',[0])
    aXp = array('f',[0])
    aYp = array('f',[0])
    aEp = array('f',[0])
    #aCh = array('f',[0])
    #aNCh = array('f',[0])
    anfev = array('f',[0])

    
    
    
    #Creates the Branches, which lay inside the tree1 and divides the data on subcategories. The /F is to indicate we're storing float variables
    bL1 = tree1.Branch('L1',aL1,'Elp[4]/F')
    bL2 = tree1.Branch('L2',aL2,'Elp[5]/F')
    bL3 = tree1.Branch('L3',aL3,'Elp[6]/F')
    bL4 = tree1.Branch('L4',aL4,'Elp[7]/F')
    bC1 = tree1.Branch('C1',aC1,'Elp[0]/F')
    bC2 = tree1.Branch('C2',aC2,'Elp[1]/F')
    bC3 = tree1.Branch('C3',aC3,'Elp[2]/F')
    bC4 = tree1.Branch('C4',aC4,'Elp[3]/F')
    bL1c = tree1.Branch('L1c',aL1c,'L1c/F')
    bL2c = tree1.Branch('L2c',aL2c,'L2c/F')
    bL3c = tree1.Branch('L3c',aL3c,'L3c/F')
    bL4c = tree1.Branch('L4c',aL4c,'L4c/F')
    bC1c = tree1.Branch('C1c',aC1c,'C1c/F')
    bC2c = tree1.Branch('C2c',aC2c,'C2c/F')
    bC3c = tree1.Branch('C3c',aC3c,'C3c/F')
    bC4c = tree1.Branch('C4c',aC4c,'C4c/F')
    bXp = tree1.Branch('Xp',aXp,'Xp/F')
    bYp = tree1.Branch('Yp',aYp,'Yp/F')
    #bCh = tree1.Branch('Ch',aCh,'Ch/F')
    #bNCh = tree1.Branch('NCh',aNCh,'NCh/F')
    bnfev = tree1.Branch('nfev',anfev,'nfev/F')
    bEp = tree1.Branch('Ep',aEp,'Ep/F')
    
    
    #Gets the variables event from event
    for entryNum in range(0, tree0.GetEntries()):    
        tree0.GetEntry(entryNum)
        
        
        print('-------------------------------')
        print(entryNum)
        
        #NCh = getattr(tree0,'NCh')
        NCh = 8
        
        L1 = getattr(tree0, 'L1')
            #print(L1)
        L2 = getattr(tree0, 'L2')
            #print(L2)
        L3 = getattr(tree0, 'L3')
            #print(L3)
        L4 = getattr(tree0, 'L4')
            #print(L4)
        C1 = getattr(tree0, 'C1')
            #print(C1)
        C2 = getattr(tree0, 'C2')
            #print(C2)
        C3 = getattr(tree0, 'C3')
            #print(C3)
        C4 = getattr(tree0, 'C4')
            #print(C4)
        Ep = getattr(tree0, 'Ep')
        
        if NCh>=7 and L1>150 and L2>300 and L3>300 and L4>150 and C1>150 and C2> 100 and C3>100 and C4>150:
            print('passou')
        
        

            calib(L1,L2,L3,L4,C1,C2,C3,C4)
       
            global newExN
        
        #Organizes the calibrated values from each line and column on an array holding two vectors (Lines and Columns)
            #newExN = np.array([[calb[0],calb[1],calb[2],calb[3]],[calb[4],calb[5],calb[6],calb[7]]])
            newExN = np.array([[calb[3],calb[2],calb[1],calb[0]],[calb[4],calb[5],calb[6],calb[7]]])
            #newExN = np.array([[calb[3],calb[2],calb[1],calb[0]],[calb[4],calb[5],calb[6],calb[7]]])
        #Associates the values on the right to the pointers on the left
            aL1[0] = L1
            aL2[0] = L2
            aL3[0] = L3
            aL4[0] = L4
            aC1[0] = C1
            aC2[0] = C2
            aC3[0] = C3
            aC4[0] = C4
            aL1c[0] = newExN[0][3]
            aL2c[0] = newExN[0][2]
            aL3c[0] = newExN[0][1]
            aL4c[0] = newExN[0][0]
            aC1c[0] = newExN[1][0]
            aC2c[0] = newExN[1][1]
            aC3c[0] = newExN[1][2]
            aC4c[0] = newExN[1][3]
            aEp[0] = Ep
            #aNCh[0] = NCh
        #aCh[0] = Ch
        
       
        #Executes the optimize least squared routine and stores the useful information on given pointers
            xypat = fitpatxy()
        
            aYp[0] = result.x[0]
            aXp[0] = result.x[1]
            anfev[0] = result.nfev

       
        #Fills the tree1 with the data of the given event
            tree1.Fill()
       
    
    
    tree1.SetDirectory(0)
    tree1.Fill()
    data.Close()
    print('Acabou')
    
    #Writes all data from every event on the tree1
    tree1.Write()
    
    #Writes everything on the exit file and closes it
    outfile.Write()
    outfile.Close()
    
    #Congratulations! It's the end of the program
