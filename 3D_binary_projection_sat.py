# Copyright 2010-2018 Google LLC
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""OR-tools solution to the 3D binary matrix projection problem."""
from __future__ import print_function
import time
import sys
from ortools.sat.python import cp_model
import numpy as np
import bin_matrix_utils as bu
from matplotlib import pyplot as plt


class BinaryProjectionSolutionPrinter(cp_model.CpSolverSolutionCallback):
    """Print solutions."""

    def __init__(self, matrix, x_sum, y_sum, z_sum):
        cp_model.CpSolverSolutionCallback.__init__(self)
        self.__x_sum = x_sum
        self.__y_sum = y_sum
        self.__z_sum = z_sum
        self.__matrix = matrix
        self.__N = x_sum.shape[0] # its a cube
        self.__solution_count = 0
        self.__start_time = time.time()

    def solution_count(self):
        return self.__solution_count

    def on_solution_callback(self):
        current_time = time.time()
        print('Solution %i, time = %f s' % (self.__solution_count,
                                            current_time - self.__start_time))
        self.__solution_count += 1

        x_sum = np.zeros((N,N),dtype=int)
        y_sum = np.zeros((N,N),dtype=int)
        z_sum = np.zeros((N,N),dtype=int)
        for x in range(N):
            print('\n')
            for y in range(N):
                for z in range(N):
                    if self.Value(self.__matrix[x,y,z]) == 1:
                        print('1', end=' ')
                    else:
                        print('0', end=' ')
                print()
         
        print('----\n')
        # compute projections
        for x in range(N):
             for y in range(N):
                for z in range(N):
                    if self.Value(self.__matrix[x,y,z]) == 1:
                        z_sum[x,y] += 1

        for x in range(N):
             for z in range(N):
                for y in range(N):
                    if self.Value(self.__matrix[x,y,z]) == 1:
                        y_sum[x,z] += 1
              
        for y in range(N):
             for z in range(N):
                for x in range(N):
                    if self.Value(self.__matrix[x,y,z]) == 1:
                        x_sum[y,z] += 1


        print("x sum error = {}".format(np.linalg.norm(x_sum-self.__x_sum)))
        print("y sum error = {}".format(np.linalg.norm(y_sum-self.__y_sum)))
        print("z sum error = {}".format(np.linalg.norm(z_sum-self.__z_sum)))
  
        # just the first solution is enough 
        if self.__solution_count > 0:
            self.StopSearch()

def main(x_sum,y_sum,z_sum):
  
    N = x_sum.shape[0] # assume a cube
     
    # Creates the solver.
    model = cp_model.CpModel()
    
    # Creates the variables of binary cube
    matrix = {}
    for x in range(N):
        for y in range(N):
            for z in range(N):
                matrix[x,y,z] = model.NewIntVar(0,1,'x[%i,%i,%i]' % (x,y,z))
    
    # Creates the constraints.

    #add up to x sums
    for y in range(N):
        for z in range(N):
            model.Add(sum(matrix[x,y,z] for x in range(N))==x_sum[y,z])
   
    # add up to y sums
    for x in range(N):
        for z in range(N):
            model.Add(sum(matrix[x,y,z] for y in range(N))==y_sum[x,z])

    # add up to z sums
    for x in range(N):
        for y in range(N):
            model.Add(sum(matrix[x,y,z] for z in range(N))==z_sum[x,y])

    ### solve model
    solver = cp_model.CpSolver()
    solution_printer = BinaryProjectionSolutionPrinter(matrix,x_sum,y_sum,z_sum)
    status = solver.SearchForAllSolutions(model, solution_printer)

    print()
    print('Statistics')
    print('  - conflicts       : %i' % solver.NumConflicts())
    print('  - branches        : %i' % solver.NumBranches())
    print('  - wall time       : %f s' % solver.WallTime())
    print('  - solutions found : %i' % solution_printer.solution_count())


if __name__ == '__main__':
    if len(sys.argv) > 1:
        N = int(sys.argv[1])
    else:
        # default cube size. 
        # For N=2, over determined
        # For N=3, exactly determined
        # for N>=4, under determined
        N = 6
    
    # create a random binary cube
    matrix = np.random.random((N,N,N)) < np.random.random(None)
    
    # compute projections
    x_sum = np.sum(matrix,0)
    y_sum = np.sum(matrix,1)
    z_sum = np.sum(matrix,2)
    
    # show them as images
    plt.figure()
    plt.subplot(131)
    plt.imshow(x_sum)
    plt.title('x sum')
    plt.subplot(132)
    plt.imshow(y_sum)
    plt.title('y sum')
    plt.subplot(133)
    plt.imshow(z_sum)
    plt.title('z sum')
 
  
    # we know that any real cube must pass the Gale-Ryser domination
    # conditions, test it here
    if bu.check_X2_star_dominates_Y2(X2=y_sum,Y2=z_sum) and \
        bu.check_X2_star_dominates_Y2(X2=x_sum,Y2=z_sum.T) and \
          bu.check_X2_star_dominates_Y2(X2=x_sum.T,Y2=y_sum.T):
              
              # solve!
              print("A solution must exist, solving...")
              main(x_sum,y_sum,z_sum)
    else:
        print("Dominance test failed: No solution exists!")
