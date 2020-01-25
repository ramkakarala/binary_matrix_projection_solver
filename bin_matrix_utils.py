#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 14:41:08 2020

A collection of utilities for binary matrices.

@author: ramkakarala
"""

import numpy as np

#According to Gale-Ryser theorem, there a matrix 
#with sum of columns "csum" and sum of rows "rsum"
#iff rsum* dominates csum.
# this routine checks dominance of csum star over rsum
# where star is the sum of rows the Ferrers matrix
# with sum of columns csum
def check_X_star_dominates_Y(X,Y):
    # make them the same length
    X_s = -np.sort(-X) # descending order
    Y_s = -np.sort(-Y)
    F = np.zeros((X.size,Y.size),dtype=int)
    for k in range(X_s.size):
        F[k,0:X_s[k]] = 1

    X_star = np.sum(F,axis=0)
    dominates = False
    if np.all(np.cumsum(X_star)>=np.cumsum(Y_s)):
        dominates=True
    return dominates

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
     
if __name__=='__main__':
    test_check_X_star_dominates_Y()
