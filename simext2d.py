###############################################################################
# Read simextNNNN.dat and simphiavgextNNNN.dat to calculate
# Trphi in the mean fluid frame defined by the phi averaged velocity
# routines which work on 2d data are in this file
###############################################################################

import numpy as np
from matplotlib import pyplot as plt

class simext2d():
    ####################################
    # sorting routines so that
    # quantities will be on an indexable
    # grid q[i,j]
    # quantities will be accessible as simext2d.q_grid
    # for quantity names q
    ####################################
    def sort_scalar2grid(self,quantity):
        #indice range
        imin = int(np.amin(self.i))
        jmin = int(np.amin(self.j))

        imax = int(np.amax(self.i))
        jmax = int(np.amax(self.j))

        grid = np.zeros((imax-imin+1, jmax-jmin+1))
        for i,j,q in zip(self.i, self.j, getattr(self, quantity)):
            grid[int(i)-imin,int(j)-jmin] = q
        setattr(self, quantity+'_grid', grid)
        return

    # works for 4 dimensional vector
    def sort_vector2grid(self,quantity):
        imin = int(np.amin(self.i))
        jmin = int(np.amin(self.j))

        imax = int(np.amax(self.i))
        jmax = int(np.amax(self.j))

        grid = np.zeros((imax-imin+1, jmax-jmin+1,4))
        for i,j,q in zip(self.i, self.j, getattr(self, quantity)):
            grid[int(i)-imin,int(j)-jmin,:] = q
        setattr(self, quantity+'_grid', grid)
        return

    #works for 4x4 D tensor
    def sort_tensor2grid(self,quantity):
        imin = int(np.amin(self.i))
        jmin = int(np.amin(self.j))

        imax = int(np.amax(self.i))
        jmax = int(np.amax(self.j))

        grid = np.zeros((imax-imin+1, jmax-jmin+1,4,4))
        for i,j,q in zip(self.i, self.j, getattr(self, quantity)):
            grid[int(i)-imin,int(j)-jmin,:,:] = q
        setattr(self, quantity+'_grid', grid)
        return

    ##############################
    # initialize the data structure
    # reads required columns
    ##############################
    def __init__(self, filename, radiation = True):
        # read the simfile
        data = np.loadtxt(filename)
        # macro to select columns as they are labeled
        # in fileop.c
        column = lambda i: data[:,i-1]

        #indices
        self.i = column(1)
        self.j = column(2)

        #coordinates
        self.r = column(4)
        self.th = column(5)

        #density and internal energy
        self.rho = column(7)
        self.uint = column(8)

        #contravariant (upper) componenet of B
        self.bcon2 = column(16)
        self.ucon0 = column(9)
        self.ucon3 = column(12)

        #contravariant (upper) componenets of four velocity
        #saved into an array of vectors
        self.ucon = np.zeros((len(self.i),4))
        for i in range(4):
            self.ucon[:,i] = column(9+i)


       



# reads a grid r[i,j], theta[i,j] and computes the metric
# components at each gridpoint
# for 2D data, the metric array will be an NixNjx4x4 dimensional array
def calc_schwarzschild_metric_tensor_on_grid(r_grid,th_grid):
    r = r_grid
    th = th_grid
    ir, ith  = r.shape
    metric_tensor = np.zeros((ir,ith,4,4))
    metric_tensor[:,:,0,0] = -(1.-2./r)
    metric_tensor[:,:,1,1] = 1./(1.-2./r)
    metric_tensor[:,:,2,2] = r*r
    metric_tensor[:,:,3,3] = r*r*np.sin(th)*np.sin(th)
    return metric_tensor

## compute the lorentz boots matrix from frames.c
# this was just copied from Koral and translated into
# python, the inputs are the metric_tensor_array
# the inverse metric tensor array, and an array of
# the contravariant componenets of the fluid four velocity
def calc_Lorentz_lab2ff_on_grid(metric_tensor, inverse_metric_tensor, ucon_grid):

    #four velocity
    ucon = ucon_grid #shortcut

    # four velocity with lower index
    ucov = np.einsum('ijlm,ijm->ijl', metric_tensor, ucon) ## lower index

    # four velocity of the stationary observer (lab frame)
    # wcov = (sqrt(-1/g^tt),0,0,0)
    alpha = np.sqrt(-1./inverse_metric_tensor[:,:,0,0])
    wcov = np.zeros(ucon_grid.shape) ## four velocity shape
    wcov[:,:,0] = -alpha
    #contravariant components
    wcon = np.einsum('ijlm,ijm->ijl',inverse_metric_tensor, wcov)  ## raise index

    #Om is the difference of two outer products
    Om = (np.einsum('ijl,ijm->ijlm', ucon, wcov) -
            np.einsum('ijl,ijm->ijlm', wcon, ucov) )

    #lorentz factor -w^mu u_mu between lab and fluid frame
    gamma = -np.einsum('ijl,ijl->ij',wcon,ucov)

    #Om sum is omega contracted with itself?
    Omsum = np.einsum('ijlm,ijmn->ijln', Om, Om)

    # genergal expression for lorentz boost
    Lorentz = np.zeros(metric_tensor.shape)
    Lorentz += np.diag([1,1,1,1]) #kronecker delta in koral
    Lorentz += 1./(1.+gamma[:,:,np.newaxis, np.newaxis])*Omsum
    Lorentz += Om

    return Lorentz


# applies the lorentz boost array
# to a contravariant vector
# everything is construct to work with
# data on the grid
# copied from koral/frames.c
def boost_vector_lab2ff_on_grid(vector, inverse_metric_tensor, Lorentz):
    boosted_vector = np.einsum('ijlm,ijm->ijl', Lorentz, vector)
    ## make orthnormal?
    alpha = np.sqrt(-1./inverse_metric_tensor[:,:,0,0])
    return boosted_vector*alpha[:,:,np.newaxis]

# same as boos_vector_lab2ff_on_grid but for tensors
# should work for double contra variant tensor (both indices up)
# copied from koral/frames.c
def boost_tensor_lab2ff_on_grid(tensor, inverse_metric_tensor, Lorentz):

    boosted = np.einsum("abik,abjl,abkl->abij", Lorentz, Lorentz, tensor)
    ## make orthnormal?
    alpha = np.sqrt(-1./inverse_metric_tensor[:,:,0,0])
    boosted[:,:,:,0] *= alpha[:,:,np.newaxis]
    boosted[:,:,0,:] *= alpha[:,:,np.newaxis]
    return boosted



