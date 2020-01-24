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
"""OR-tools solution to the binary matrix projection problem."""
from __future__ import print_function
import time
import sys
from ortools.sat.python import cp_model
import numpy as np


class BinaryProjectionSolutionPrinter(cp_model.CpSolverSolutionCallback):
    """Print solutions."""

    def __init__(self, matrix, row_sum, col_sum):
        cp_model.CpSolverSolutionCallback.__init__(self)
        self.__row_sum = row_sum
        self.__col_sum = col_sum
        self.__matrix = matrix
        self.__solution_count = 0
        self.__start_time = time.time()

    def solution_count(self):
        return self.__solution_count

    def on_solution_callback(self):
        current_time = time.time()
        print('Solution %i, time = %f s' % (self.__solution_count,
                                            current_time - self.__start_time))
        self.__solution_count += 1

        matrix_size = int(np.sqrt(len(self.__matrix)))
        row_sum = np.zeros(matrix_size,dtype=int)
        col_sum = np.zeros(matrix_size,dtype=int)
        for i in range(matrix_size):
            for j in range(matrix_size):
                if self.Value(self.__matrix[i,j]) == 1:
                    # There is a queen in column j, row i.
                    print('1', end=' ')
                    col_sum[i] += 1
                    row_sum[j] += 1
                else:
                    print('0', end=' ')
            print()

        print("row sum - target row sum error =")
        print(np.linalg.norm(row_sum-self.__row_sum))
        print("col sum - target col sum error =")
        print(np.linalg.norm(col_sum-self.__col_sum))
    
         
        if self.__solution_count > 0:
            self.StopSearch()

def main(row_sum,col_sum):

    print("{} {}".format(row_sum.size,col_sum.size))
    
    # Creates the solver.
    model = cp_model.CpModel()
    # Creates the variables of binary matrix
    matrix = {}
    for i in range(col_sum.size):
        for j in range(row_sum.size):
            matrix[i,j] = model.NewIntVar(0,1,'x[%i,%i]' % (i,j))
    
    # Creates the constraints.

    # add up to col sums
    for i in range(col_sum.size):
        model.Add(sum(matrix[i,j] for j in range(col_sum.size))==col_sum[i])
      
     # add up to row sums
    for i in range(row_sum.size):
        model.Add(sum(matrix[j,i] for j in range(row_sum.size))==row_sum[i])
 

    ### solve model
    solver = cp_model.CpSolver()
    solution_printer = BinaryProjectionSolutionPrinter(matrix,row_sum,col_sum)
    status = solver.SearchForAllSolutions(model, solution_printer)

    print()
    print('Statistics')
    print('  - conflicts       : %i' % solver.NumConflicts())
    print('  - branches        : %i' % solver.NumBranches())
    print('  - wall time       : %f s' % solver.WallTime())
    print('  - solutions found : %i' % solution_printer.solution_count())


# set up the matrix
matrix_size = 32

if __name__ == '__main__':
    if len(sys.argv) > 1:
        matrix_size = int(sys.argv[1])
    board = np.random.random((matrix_size,matrix_size)) < np.random.random(1)[0]
    # compute projection
    rsum = np.sum(board,0)
    csum = np.sum(board,1)
    main(rsum,csum)
