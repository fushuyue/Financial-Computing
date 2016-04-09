import numpy as np
import numpy.random as npr
import scipy as sc
from scipy.linalg import eig
class Markov(object):
    def __init__(self,rho,transitionMartix):
        self.rho = rho
        self.transitionMatrix = transitionMartix
        self.numberOfStates = len(rho)

    def print_all(self):
        print('initial condition =  ')
        print(self.rho)
        print('transition matrix = ')
        print(self.transitionMatrix)
        # print('length =%d' % self.numberOfStates)
        print('\n')

    def state(self,n):
        nstate = np.asarray(self.rho*(self.transitionMatrix**n))
        nstate = nstate[0]
        #print('nstate transition matrix= ')
        #print(nstate)
        p=npr.random(1)
        #print(p)
        if p < nstate[0]:
            return 1
        else:
            for i in range(1,self.numberOfStates):
                if p <= sum(nstate[0:i+1]):
                    return i+1

    def eign(self):
        w,v=eig(self.transitionMatrix,left=True,right=False)
        w=w.real
        vector=v[:,np.abs(w-1)<0.000001]
        rho=vector/vector.sum()
        rho=rho.T
        print('the invariant distribution =')
        print(rho)


initialCondition = list(input('Please key in initial condition like "0.2,0.3,0.5" '))
temp=list(input('Please key in transition matrix,use comma to separate like "[0.1,0.2,0.7],[0.3,0.4,0.3],[0.5,0.1,0.4]" '))
transitionMatrix = np.matrix(temp)
# n=input('Please key in how many steps ahead do you want to generate  ')
# count=input('Please key in how many numbers do you want to generate  ')
n=1
'''
initialCondition=[0.2,0.3,0.5]
transitionMatrix = np.matrix([[0.1,0.2,0.7],[0.3,0.4,0.3],[0.5,0.1,0.4]])
n=10
'''

x = Markov(initialCondition,transitionMatrix)
x.print_all()

print('the simulated X is %d ' %x.state(n))

#for i in range(count):
#   print(x.state(n))

print('\n')

# the invariant distribution
x.eign()

