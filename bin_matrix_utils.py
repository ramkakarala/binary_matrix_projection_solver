#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 14:41:08 2020

A collection of utilities for binary matrices.

@author: ramkakarala
"""

import numpy as np

# a Ferrers matrix is binary with all ones in upper left corner,
# and with decreasing sum along rows
def ferrers(X,Y):
    F = np.zeros((X.size,np.amax(X)),dtype=int)
    for k in range(X.size):
        F[k,0:X[k]]=1

    return F

#According to Gale-Ryser theorem, there a matrix 
#with sum of columns "csum" and sum of rows "rsum"
#iff rsum* dominates csum.
# this routine checks dominance of csum star over rsum
# where star is the sum of rows the Ferrers matrix
# with sum of columns csum
def check_X_star_dominates_Y(X,Y):
    # sort them
    X_s = -np.sort(-X) # descending order
    Y_s = -np.sort(-Y)
    F = ferrers(X_s,Y_s)
    
    X_star = np.sum(F,axis=0)
    dominates = False
    if X_star.size > Y_s.size:
        # pad Y_S
        Y_s = np.pad(Y_s,(0,X_star.size-Y_s.size),'constant',constant_values=(0))
    elif Y_s.size > X_star.size:
        # pad X_star
        X_star = np.pad(X_star,(0,Y_s.size-X_star.size),'constant',constant_values=(0))
      
    if np.all(np.cumsum(X_star)>=np.cumsum(Y_s)):
        dominates=True
    return dominates

# same as check_X_star_dominates_Y, but for 2D matrices
def check_X2_star_dominates_Y2(X2,Y2):
    
    #only deal with cubes for now
    assert X2.shape == Y2.shape
    
    dominates = True
    for r in range(X2.shape[0]):
        if not check_X_star_dominates_Y(X2[r,:],Y2[r,:]):
            dominates = False
            break
        
    return dominates

# check all 3 projections to see if a cube exists
def check_cube_projections(axis0sum,axis1sum,axis2sum):
    
    passes = False
    if check_X2_star_dominates_Y2(X2=axis2sum,Y2=axis1sum) and \
        check_X2_star_dominates_Y2(X2=axis0sum,Y2=axis2sum.T) and \
            check_X2_star_dominates_Y2(X2=axis0sum.T,Y2=axis1sum.T):
           passes = True
    
    return passes

#unit test for check_X_dominates_Y
def test_check_X_star_dominates_Y():
    A = np.random.randint(0,2,(5,10))
    sum_of_rows = np.sum(A,axis=0)
    sum_of_cols = np.sum(A,axis=1) 
    cols_vs_rows = check_X_star_dominates_Y(X=sum_of_cols,Y=sum_of_rows)
    rows_vs_cols = check_X_star_dominates_Y(X=sum_of_rows,Y=sum_of_cols)
    if cols_vs_rows and rows_vs_cols:
        print("Dominates test for actual binary matrix passed!")
    else:
        print("Dominates test for actual binary matrix FAILED!")
    # create a case where it will fail
    fake_cols_vs_rows = check_X_star_dominates_Y(X=sum_of_cols,Y=sum_of_rows+100)
    if fake_cols_vs_rows==False:
        print("Dominates test for fake binary matrix passed!")
    else:
        print("Dominates test for fake binary matrix FAILED!")
        
#unit test for dominates condition, but for matrices
def test_check_X2_star_dominates_Y2():
    A = np.random.randint(0,2,(2,2,2))
    # assume A is indexed by y,x,z, in left handed coordinates
    zsum = np.sum(A,axis=0) # sum along z, as function of y and x
    ysum = np.sum(A,axis=1) # sum along y, as function of x and z
    xsum = np.sum(A,axis=2) # sum along x, as function oy y and z

    if check_cube_projections(axis0sum=zsum,axis1sum=ysum,axis2sum=xsum):
        print("Dominates test for actual binary cube passed!")
    else:
        print("Dominates test for actual binary cube FAILED!")
       
     
if __name__=='__main__':
    test_check_X_star_dominates_Y()
    test_check_X2_star_dominates_Y2()
