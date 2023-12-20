import numpy as np
import scipy.optimize as spopt
import ImageModelClass as IMC
import ExpPat as EX
import ExpReCal as ERC
import scipy.sparse as sp
import matplotlib.pyplot as plt



#Initial values for the optimization parameters
b=np.ones((2,4))
c=np.zeros((2,4))
eps0=0.5
IM=IMC.ImageModel(eps0)

YpixModel=IM.YpixModel
YPattern=IM.YPattern


#List of experimental data point positions
xtab=EX.xtab
ytab=EX.ytab
Nx=xtab.size
Ny=ytab.size


#Starting points for the optimization process
X = np.array([])
X = np.append(X, eps0)  
X = np.append(X, np.ravel(b))
X = np.append(X, np.ravel(c))




# Define or set ExDataTab and RelErrorDat
ExDataTab=EX.ExDataTab
RelErrorDat=EX.RelErrorDat
ExN=ERC.normdata(ERC.recaldata(ExDataTab,b,c))

Y_dat=ExDataTab
Y_sig=RelErrorDat

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

def errorfunc(xoff,yoff,dTab,YTab): #dTab is the data table after calib.+ renorm. for example
	return YTab-dTab

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

def estimate_jacobian_sparsity(func, x, *args):
    # Calculate the number of residuals (size of the function's output)
    n = len(func(x, *args))

    # Calculate sparsity pattern
    sparsity_pattern = sp.lil_matrix((n, len(x)))

    for i in range(len(x)):
        # Perturb one element at a time
        x_perturbed = np.array(x, dtype=float)
        x_perturbed[i] += 1e-8  # Small value for perturbation

        # Calculate the change in function outputs
        residuals_diff = func(x_perturbed, *args) - func(x, *args)

        # Find non-zero elements
        non_zero_elements = np.nonzero(residuals_diff)[0]

        # Update the sparsity pattern matrix
        for row in non_zero_elements:
            sparsity_pattern[row, i] = 1

    return sparsity_pattern


sparsity_pattern = estimate_jacobian_sparsity(fbc, X, Y_dat, Y_sig)
dense_sparsity = sparsity_pattern.toarray()

file_path = "sparsity_pattern_matrix.txt"

# Save the sparsity pattern matrix to a text file
np.savetxt(file_path, dense_sparsity, fmt='%d')

print(dense_sparsity)


plt.imshow(dense_sparsity, cmap='binary', aspect='auto')
plt.xlabel('Variables')
plt.ylabel('Residuals')
plt.title('Jacobian Sparsity Pattern')
plt.show()