import numpy as np
class linProblem(object):
    # for Primal problem , i.e.
    #         max cy - constance
    # subject to  Ay <= b
    #             yi >= 0
    
    def __init__ (self,target,cons,A,b):
        
        #Problem var
        
        self.target = target.copy()
        self.cons = cons
        self.A = A.copy()
        self.b = b.copy()
        self.dimension = A.shape
        
        n = self.dimension[0] #Row number
        m = self.dimension[1] #Colunm number , same with unknown var yi
        
        #Solution var
        
        self.basicSol = np.zeros(m).astype('object')
        self.basicInd = np.arange(1,m+1).astype('object')
        self.indepSol = np.zeros(n).astype('object')
        self.indepInd = np.arange(1001,n+1001).astype('object')
        self.maxValue = 0
        
        #setup simplex tableau 
        
        self.tableau = np.zeros([n+1, m+1])
        self.tableau[0:n,0:m] = self.A.copy()
        self.tableau[n,0:m] = self.target.copy()
        self.tableau[0:n,m] = self.b.copy()
        self.tableau[n,m] = -self.cons
        self.tableau_r = self.tableau.copy()
        #final tableau as fraction form
        self.tableau_r_f = self.tableau.copy()
        
        #Calculation stat
        self.iterationCount = 0
        
        
    def display(self):      
        
        #display the problem message
        
        print("Objective function coefficint:")
        if (self.cons == 0):
            print(self.target)
        else:
            print("{target} - {cons}".format(target = self.target, cons = self.cons))       
        print("Subjection:")
        print(self.A,"<=")
        print(self.b)
   
            
    def dpMat_as_fraction(self, mat_, flag = "NoReturn"):
        #display matrix as fraction form
        tempMat = mat_.copy()
        tempMat = tempMat.astype('object')
        from fractions import Fraction
        from sympy import Rational
        n = mat_.shape[0]
        m = mat_.shape[1]
        for indexi in range(n):
            for indexj in range(m):
                tempMat[indexi,indexj] = Fraction.from_float(tempMat[indexi,indexj]).limit_denominator(1000)
                tempMat[indexi,indexj] = Rational(tempMat[indexi,indexj].numerator,tempMat[indexi,indexj].denominator)
        print(tempMat)
        
        if (flag == "Return"):
            return tempMat
        
    def displayResult(self):
        
        n = self.dimension[0]
        m = self.dimension[1]
        
        print("The maximum value of target function is:",self.maxValue)
        print("In fractional form :",-self.tableau_r_f[n,m])
        print("The optimal solution is:",self.basicSol)
        print("The optimal solution to dual problem is:",self.indepSol)
        print("---")
        print("Iteration Count: ",self.iterationCount)
        print("---")

        return

    def solve(self, flag):
        
        #solve the problem by simplex method
        n = self.dimension[0]
        m = self.dimension[1]
        if ((self.tableau_r[n,0:m]<=np.zeros([m])).all()):
            
            print("The Final simplex tableau:")
            self.tableau_r_f = self.dpMat_as_fraction(self.tableau_r,"Return")
            print("")

            #deal with index and solution value
            self.maxValue = -self.tableau_r[n,m]
            for i in range(m):
                if (self.basicInd[i]>1000):
                    self.indepSol[self.basicInd[i]-1001] = -self.tableau_r_f[n,i]
            for i in range(n):
                if (self.indepInd[i]<1000):
                    self.basicSol[self.indepInd[i]-1] = self.tableau_r_f[i,m]            
            return
        
        else:
            
            self._iteration(flag)
            self.solve(flag)    
            
    def _iteration(self, flag):
        #sub iteration function of self.solve
        tableau_ = self.tableau_r.copy()
        tableau_new = tableau_.copy()
        
        n = self.dimension[0]
        m = self.dimension[1]
        i = 0
        j = 0
        #Counting the iteration times
        self.iterationCount = self.iterationCount + 1
        
        if (flag == 1):
            
            print("Step ",self.iterationCount,":")
            self.dpMat_as_fraction(tableau_)
        
        for j in range(m):            
            if (tableau_[n,j] > 0):
                temp =(tableau_[0:n,m]+1e-10)/(tableau_[0:n,j]+1e-10)
                i = np.where(temp >= 0, temp, np.inf).argmin()
                break
              
        d = tableau_[i,j]
        tableau_new[i,:] = tableau_new[i,:]/d
        tableau_new[:,j] =-tableau_new[:,j]/d
        tableau_new[i,j] = 1/d
        
        for i_ in range(n+1):
            for j_ in range(m+1):
                if( i_ != i and j_ != j):
                    tableau_new[i_,j_] = (tableau_[i_,j_]*d - tableau_[i,j_]*tableau_[i_,j] )/d
                    
        #exchange the var index
        
        self.basicInd[j],self.indepInd[i] = self.indepInd[i],self.basicInd[j]
        
        #print the central point
                    
        if (flag == 1):
            print("i,j(math):",i+1,",",j+1)
            print("")
            
        #update the matrix
        
        self.tableau_r = tableau_new
        

  
                
        
        
        
        