import numpy as np
import random as random
import time
import sys

def compute_N():
 Sum = 0
 for site in S:
  for bound_motors in cs[site]:
   Sum += ms[site]*(Nb - bound_motors)
 return Sum

def initialize_motor_distributions(M,L):
 global u  
 m = [0 for i in range(L)]
 temp = M
 if M < L:
  for i in range(L):
   if temp == 0:
    break
   temp -=1
   m[i] = 1
 if M >= L:
  m = [int(M/L) for i in range(L)]
  temp -= L*int(M/L)
  if temp > 0:
   for i in range(L):
    if temp == 0:
     breakx
    temp -=1
    m[i] += 1
 b = [ind for ind,val in enumerate(m) if val != 0 and ind > 0 and ind < L-1]
 u = sum(m[1:L-1])
 global m0
 m0 = 1 if m[0] > 0 else 0
 return m,b


# add the following lines.
