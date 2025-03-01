from __future__ import division 
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
from array import array
from scipy import linalg
import pandas as pd
import openpyxl



#eps0=1.
#Ybg0=1.E-4


tab_eps0 = "./par_eps_Mean_inc_invl_2reflex_switch.txt"
eps0 = np.loadtxt(tab_eps0)

#tab_Ybg0 = "./par_Ybg0_2.txt"
#Ybg0 = np.loadtxt(tab_Ybg0)

#Import of calibration parameters b calculated by fitb() from Fits.py
tab_b="./par_b_Mean_inc_invl_2reflex_switch.txt"

b=np.loadtxt(tab_b)
#b=np.array([1,1,1,1,1,1,1,1])

tab_c="./par_c_Mean_inc_invl_2reflex_switch.txt"

c=np.loadtxt(tab_c)

#c = np.array([0,0,0,0,0,0,0,0])

IM=IMC.ImageModel(eps0)

Lpitch = 4.2

ruler=Lpitch*np.array([-1.5,-0.5,0.5,1.5])


####################################################################

#Opens the root file with the data to be processed
data_address = 'C:/Users/Usuario/Desktop/Incidence-Position-On-Phoswich-Detector/4_2.xlsx'
data = pd.read_excel(data_address)
print(data.columns)


#Calculates the Pattern Error Function using xguess, yguess and the Yield Data
def fpatxy(X,Y_dat): 
	x=X[0]
	y=X[1]
	f=ft.PatternErrorOffset(x,y,-3.5,-2.1,Y_dat)   #<====PatternErrorOffset??
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
    IM.eps=np.loadtxt("par_eps_Mean_inc_invl_2reflex_switch.txt")
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
    calL1 = (l1)*b[0]+c[0]
    calL2 = (l2)*b[1]+c[1]
    calL3 = (l3)*b[2]+c[2]
    calL4 = (l4)*b[3]+c[3]
    sumL = calL1 + calL2 + calL3 + calL4
    
    
    #Calibrated columns
    
    calC1 = (c1)*b[4]+c[4]
    calC2 = (c2)*b[5]+c[5]
    calC3 = (c3)*b[6]+c[6]
    calC4 = (c4)*b[7]+c[7]
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

    #Create empty dataframe for data, calibrated data and positions
    output_columns = ['L1','L2','L3','L4','C1','C2','C3','C4','cL1','cL2','cL3','cL4','cC1','cC2','cC3','cC4','X','Y']
    response = pd.DataFrame(columns=output_columns)


    #Gets the variables event from event
    for entryNum in range(0, len(data)):    
        current_data = data.iloc[entryNum]
        
        
        print('-------------------------------')
        print(entryNum)
        L1, L2, L3, L4, C1, C2, C3, C4 = current_data[['l1', 'l2', 'l3', 'l4', 'c1','c2', 'c3', 'c4']].values.T.tolist()

        #if L1>150 and L2>300 and L3>300 and L4>150 and C1>150 and C2> 100 and C3>100 and C4>150:
        if L1>0 and L2>0 and L3>0 and L4>0 and C1>0 and C2>0 and C3>0 and C4>0:
            calib(L1,L2,L3,L4,C1,C2,C3,C4)  
            global newExN

            #Organizes the calibrated values from each line and column on an array holding two vectors (Lines and Columns)
            newExN = np.array([[calb[3],calb[2],calb[1],calb[0]],[calb[4],calb[5],calb[6],calb[7]]])

            #Executes the optimize least squared routine to calculate position
            xypat = fitpatxy()
            new_row = {'L1':L1,'L2':L2,'L3':L3,'L4':L4,'C1':C1,'C2':C2,'C3':C3,'C4':C4,'cL1':newExN[0][3],'cL2':newExN[0][2],'cL3':newExN[0][1],'cL4':newExN[0][0],'cC1':newExN[1][0],'cC2':newExN[1][1],'cC3':newExN[1][2],'cC4':newExN[1][3],'X':result.x[1],'Y':result.x[0]}
            response = response.append(new_row, ignore_index=True)
    print("It is done.")
    response.to_excel('output_4_2_Mean_inc_invl_2reflex_switch.xlsx', index=False)