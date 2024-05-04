import numpy as np
import random as random
import sys

def initialize_motor_distributions(M,L):
 m = np.zeros(L,int)
 m[:] = M//L
 m[0:M%L] += 1
 return m
def initialize_cargo_distributions(L):
 c = np.zeros([L,Nb+1],int)
 return c
def compute_N():
 Sum = 0
 for site in cargo_list:
  Sum += ms[site]*(np.sum(cs[site])*Nb - np.sum(temp_list[site]))
 return Sum  
def add_cargo():
 global m0,Left,B,C,temp_list
 motors_to_bind = 1 if np.minimum(ms[0],Nb) == 1 else random.randint(1,np.minimum(ms[0],Nb))
 ms[0] -= motors_to_bind
 Left -= motors_to_bind
 B += motors_to_bind
 C += 1
 m0 = 1 if ms[0] > 0 else 0
 cs[0,motors_to_bind] += 1
 temp_list[0,motors_to_bind] += motors_to_bind
 if 0 not in cargo_list:
  cargo_list.append(0)
def diffuse_unbound_bulk_motor(P):
 global U
 temp = P
 factor = 2*p*ratio
 for site in range(1,L-1):
  if R < temp + factor*ms[site]:
   ms[site] -= 1  
   if bool(random.getrandbits(1)):
    # Moving a motor to the left.
    ms[site-1] += 1
    if site == 1:
     global m0,Left
     m0 = 1
     Left += 1
     U -= 1
    return
   # Else, moving a motor to the right.
   ms[site+1] += 1
   if site == L-2:
    global Right
    Right += 1
    U -= 1
   return
  temp += factor*ms[site]  
def walk(P):
 global temp_list
 temp = P
 for site in range(L):
  if R < temp + ratio*np.sum(cs[site,1:Nb]):
   indx = random.choice(np.nonzero(cs[site,1:Nb])[0] + 1)
   cs[site,indx] -= 1
   temp_list[site,indx] -= indx
   if np.sum(cs[site]) == 0:
    cargo_list.remove(site)    
   if site < L-1:
    cs[site+1,indx] += 1
    temp_list[site+1,indx] += indx
    if (np.sum(cs[site+1]) > 0) and (site + 1 not in cargo_list):
     cargo_list.append(site+1)
    return 'step'
   if site == L-1:
    global Right,B,C
    ms[L-1] += indx
    Right += indx
    B -= indx
    C -= 1
    return 'success'  
  temp += ratio*np.sum(cs[site,1:Nb]) 
def unbind(P):
 global B,temp_list
 B -= 1
 temp = P
 factor = k_off*ratio
 for site in range(L):
  if R < temp + factor*np.sum(temp_list[site]):   
   for ind,val in enumerate(temp_list[site]):
    if R < temp + factor*val:
     ms[site] += 1
     cs[site,ind] -= 1
     cs[site,ind-1] += 1
     temp_list[site,ind] -= ind
     temp_list[site,ind-1] += ind - 1
     if site == 0:
      global m0,Left
      Left += 1
      m0 = 1
      if ind == 1:
       return(site) 
      return(-1)
     if site == L-1:
      global Right
      Right += 1
      if ind == 1:
       return(site) 
      return(-1)        
     global U
     U += 1
     if ind == 1:
      return(site) 
     return(-1)
    temp += factor*val    
  temp += factor*np.sum(temp_list[site])
def bind(P):
 global B,temp_list
 B += 1
 new_P = P
 factor = k_on*ratio
 for site in range(L):
  factor_prime = factor*ms[site]
  temp = factor_prime*(np.sum(cs[site])*Nb - np.sum(temp_list[site]))
  if R < new_P  + temp:
   for ind,val in enumerate(temp_list[site]):
    temp = factor_prime*(cs[site,ind]*Nb - val)
    if R < new_P + temp:
     ms[site] -= 1        
     cs[site,ind] -= 1
     cs[site,ind+1] += 1
     temp_list[site,ind] -= ind
     temp_list[site,ind+1] += (ind + 1)
     #***************
     if site == 0:
      global m0,Left
      Left -= 1
      m0 = 1 if ms[0] > 0 else 0
      return
     #***************
     if site == L-1:
      global Right
      Right -= 1
      return
     global U
     U -= 1
     return
    new_P += temp
  new_P  += temp    
def diffuse_left_site():
 global Left, m0, U
 Left -= 1  
 U += 1
 ms[0] -= 1
 ms[1] += 1
 m0 = 1 if ms[0] > 0 else 0
def diffuse_right_site():
 global Right,U
 ms[L-1] -= 1
 ms[L-2] += 1
 Right -= 1
 U += 1 
#*************************
# Main body of the code. *
#*************************
list_of_alphas=[0.1,1,10,100,1000]
Motors = int(sys.argv[1])
alpha_indx = int(sys.argv[2])
alpha = list_of_alphas[alpha_indx]

Nb = 3
L = 100

bind_list = list(range(Nb+1))

p = 654
k_walk = 1.
k_on = 0.0125 
k_off = 0.00687

ms = initialize_motor_distributions(Motors,L)
cs = initialize_cargo_distributions(L)

ms_file = "ms_m_"+str(Motors)+"_alpha_"+str(alpha)+".txt"
bound_zero_file = "zero_m_"+str(Motors)+"_alpha_"+str(alpha)+".txt"
bound_one_file = "one_m_"+str(Motors)+"_alpha_"+str(alpha)+".txt"
bound_two_file = "two_m_"+str(Motors)+"_alpha_"+str(alpha)+".txt"
bound_three_file = "three_m_"+str(Motors)+"_alpha_"+str(alpha)+".txt"

falls = "falls_m_"+str(Motors)+"_alpha_"+str(alpha)+".txt"
success = "success_m_"+str(Motors)+"_alpha_"+str(alpha)+".txt"

file1 = open(ms_file, "w")
file2 = open(bound_zero_file, "w")
file3 = open(bound_one_file, "w")
file4 = open(bound_two_file, "w")
file5 = open(bound_three_file, "w")
file6 = open(falls, "w")
file7 = open(success, "w")

measurements = 10000000
measurement_interval = 100000
final_time = measurements*measurement_interval
next_measurement = 0
step = 0

m0 = 1 if ms[0] > 0 else 0
U = sum(ms[1:L-1])
Left = ms[0]
Right = ms[L-1]
temp_list = cs*bind_list
B = np.sum(temp_list)
C = np.sum(cs[:,1:Nb])
cargo_list = []

while step <= final_time:
 #**********************************
 # Computing necessary quantities. *
 #**********************************
 if step == next_measurement:
  next_measurement += measurement_interval
  #*********************************************
  file1.write(' '.join(str(m) for m in ms))
  file1.write('\n')
  file1.flush()
  #*********************************************
  file2.write(' '.join(str(c) for c in cs[:,0]))
  file2.write('\n')
  file2.flush()
  #*********************************************
  file3.write(' '.join(str(c) for c in cs[:,1]))
  file3.write('\n')
  file3.flush()
  #*********************************************
  file4.write(' '.join(str(c) for c in cs[:,2]))
  file4.write('\n')
  file4.flush()
  #*********************************************
  file5.write(' '.join(str(c) for c in cs[:,3]))
  file5.write('\n')
  file5.flush()
  
 step += 1
 N = compute_N()
 
 ratio = 1./(alpha*m0 + 2*p*U + p*(Left + Right) + k_off*B + k_on*N + C)
 R = random.random() 
 P = alpha*m0*ratio
 if R < P:
  add_cargo()
  continue
 if R < P + 2*p*U*ratio:
  diffuse_unbound_bulk_motor(P)
  continue
 P += 2*p*U*ratio 
 if R < P + p*ms[0]*ratio:
  diffuse_left_site()
  continue
 P += p*ms[0]*ratio 
 if R < P + p*ms[L-1]*ratio:
  diffuse_right_site()
  continue
 P += p*ms[L-1]*ratio 
 if R < P + k_off*B*ratio:
  unbind_return = unbind(P)
  if unbind_return > -1:
   file6.write(str(unbind_return)+' '+str(step-1)+'\n')
   file6.flush() 
   continue
  continue   
 P += k_off*B*ratio
 if R < P + k_on*N*ratio:
  bind(P)
  continue
 P += k_on*N*ratio
 if R < P + C*ratio:
  if walk(P) == 'success':
   file7.write(str(step-1)+'\n')
   file7.flush()
   continue
  continue 

file1.close()
file2.close()
file3.close()
file4.close()
file5.close()
file6.close()
file7.close()
