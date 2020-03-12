# SimplexMethod
A python scripts for solving linear programming problem by simplex method.

    # for Primal problem , i.e.
    #         max cy - constance
    # subject to  Ay <= b
    #             yi >= 0
    
Example command:

import simplexMethod as sim
import numpy as np

#Testing Problem
cost = np.array([1,1,1])
mat = np.array([[5,0,6],[1,6,4],[4,4,8]])
res = np.array([1,1,1])
c = 0
TestingProblem = sim.linProblem(cost,c,mat,res)

#display the problem
TestingProblem.display()
print("")

#show the final simplex tableau
TestingProblem.solve(0)
TestingProblem.displayResult()

#show the calculation step by step, using *.solve(1)
TestingProblem2 = sim.linProblem(cost,c,mat,res)
TestingProblem2.solve(1)
TestingProblem2.displayResult()

