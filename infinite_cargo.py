import numpy as np
import random as random
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
     break
    temp -=1
    m[i] += 1
 b = [ind for ind,val in enumerate(m) if val != 0 and ind > 0 and ind < L-1]
 u = sum(m[1:L-1])
 global m0
 m0 = 1 if m[0] > 0 else 0
 return m,b
def initialize_cargo_distributions(L):
 global b,num_of_cargo
 c = [[] for i in range(L)]
 C = [ind for ind,val in enumerate(c) if len(val) > 0]
 for ind in C:
  num_of_cargo += len(c[ind])
 b = len(C)
 return c,C
def initialize_S_distributions():
 if len(ms) >= len(C):
  S = [ind for ind in C if ms[ind] > 0 and sum(cs[ind]) < len(cs[ind])*Nb]
 S = [ind for ind,val in enumerate(ms) if val > 0 and sum(cs[ind]) < len(cs[ind])*Nb]
 return S 
def add_cargo():
 motors_to_bind = 1 if np.minimum(ms[0],Nb) == 1 else random.randint(1,np.minimum(ms[0],Nb))
 ms[0] -= motors_to_bind
 cs[0].append(motors_to_bind)
 
 global b,num_of_cargo,m0
 num_of_cargo += 1 
 b += motors_to_bind  
    
 if 0 not in C:
  C.append(0)
  if ms[0] > 0 and motors_to_bind < Nb:
   S.append(0)
   return
  if ms[0] == 0:
   m0 = 0    
  return
 if ms[0] > 0 and motors_to_bind < Nb and 0 not in S:
  S.append(0)
  return
 if ms[0] == 0: 
  m0 = 0
  if 0 in S:
   S.remove(0)
  return
 if sum(cs[0]) == len(cs[0])*Nb and 0 in S:
  S.remove(0)
 return
def move_bulk_motor_left(site):  
 l_site = site - 1
 ms[site] -= 1
 ms[l_site] += 1
 if ms[site] == 0:
  B.remove(site)
  if site in S:
   S.remove(site)
 if sum(cs[l_site]) < len(cs[l_site])*Nb and l_site not in S:
   S.append(l_site)  
 if site > 1:
  if l_site not in B:
   B.append(l_site)
  return 
 global u,m0
 u -= 1
 m0 = 1    
def move_bulk_motor_right(site): 
 r_site = site + 1
 ms[site] -= 1
 ms[r_site] += 1
 if ms[site] == 0:
  B.remove(site)
  if site in S:    
   S.remove(site)
 if sum(cs[r_site]) < len(cs[r_site])*Nb and r_site not in S:
  S.append(r_site)
 if r_site < L-1:
  if r_site not in B:
   B.append(r_site)
  return
 global u
 u -= 1
def move_left_boundary_motor(): 
 ms[0] -= 1
 ms[1] += 1
 global u
 u += 1 
 if ms[0] == 0:
  global m0
  m0 = 0 
  if 0 in S:
   S.remove(0)
 if 1 in B:
  if sum(cs[1]) < len(cs[1])*Nb and 1 not in S:
   S.append(1)
  return
 B.append(1) 
def move_right_boundary_motor():  
 site = L-1
 l_site = L-2 
 ms[site] -= 1
 ms[l_site] += 1
 global u
 u += 1
 if ms[site] == 0 and site in S:  
  S.remove(site) 
 if l_site in B:
  if sum(cs[l_site]) < len(cs[l_site])*Nb and l_site not in S:
   S.append(l_site)
  return  
 B.append(l_site) 
def unbind_motor(site,indx):   
 global b,num_of_cargo
 b -= 1
 ms[site] += 1
 cs[site][indx] -= 1
    
 if site > 0:
  if site < L-1:
   global u
   u += 1
   if site not in B:
    B.append(site)
  if cs[site][indx] == 0:
   num_of_cargo -= 1 
   cs[site].pop(indx);
   if len(cs[site]) == 0:
    C.remove(site)
    if site in S:
     S.remove(site)   
   return site if step >= relax_time else -1
  if site not in S:
   S.append(site)
  return -1

 global m0
 m0 = 1
 if cs[0][indx] == 0:
  num_of_cargo -= 1 
  if len(cs[0]) == 0:
   C.remove(0)
   if 0 in S:
    S.remove(0) 
  return 0 if step >= relax_time else -1
 if 0 not in S:
  S.append(0)
 return -1
def bind_motor(site,indx):
 global b
 ms[site] -= 1
 cs[site][indx] += 1
 b += 1
 if site > 0 and site < L-1:
  global u
  u -= 1
  if ms[site] == 0:
   B.remove(site)
   S.remove(site) 
   return
  if sum(cs[site]) == len(cs[site])*Nb:
   S.remove(site)
  return
 if ms[site] == 0 or sum(cs[site]) == len(cs[site])*Nb:
  S.remove(site)
 if site == 0 and ms[0] == 0: 
  global m0
  m0 = 0
 return
def move_a_cargo(site):
 indx = random.choice(range(len(cs[site])))  # choosing the cargo within the site to update.
 if len(cs[site]) == 1:
  C.remove(site)
  if site in S:
   S.remove(site)
 if site < L-1:
  r_site = site + 1
  if cs[site][indx] < Nb and ms[r_site] > 0 and r_site not in S:
   S.append(r_site)   
  cs[r_site].append(cs[site][indx])
  cs[site].pop(indx)
  if r_site not in C: 
   C.append(r_site)   
  if site in S and sum(cs[site]) == len(cs[site])*Nb:
   S.remove(site)
  return
 global b,num_of_cargo
 num_of_cargo -= 1 
 b -= cs[site][indx]
 ms[site] += cs[site][indx]
 cs[site].pop(indx)
 if site in S and sum(cs[site]) == len(cs[site])*Nb:
  S.remove(site)
 return 0 if step >= relax_time else -1
def diffuse_unbound_bulk_motor():
 const1 = R*D/(2*p)
 const2 = P*D/(2*p)
 for indx in B:
  temp = ms[indx]
  if const1 < const2 + temp:
   if bool(random.getrandbits(1)):
    move_bulk_motor_left(indx)
    return
   move_bulk_motor_right(indx)
   return
  const2 += temp
def find_motor_to_unbind():
 const1 = R*D/k_off
 const2 = P*D/k_off
 for site in C:
  temp = sum(cs[site])
  if const1 < const2 + temp:
   for indx,val in enumerate(cs[site]):
    if const1 < const2 + val:
     pass_indx = unbind_motor(site,indx)
     return pass_indx
    const2 += val
  const2 += temp
def find_motor_to_bind():
 const1 = R*D/k_on
 const2 = P*D/k_on
 for site in S:
  temp = ms[site]*(len(cs[site])*Nb - sum(cs[site]))
  if const1 < const2 + temp:
   for indx,val in enumerate(cs[site]):
    if const1 < const2 + val: 
     bind_motor(site,indx)
     return
    const2 += val
  const2 += temp    
def find_cargo_to_step():
 const1 = R*D/k_on
 const2 = P*D/k_on
 for site in C:
  temp = len(cs[site])
  if const1 < const2 + temp:
   pass_val = move_a_cargo(site)
   return pass_val
  const2 += temp
#**************************************************  
# This is where the main part of the code starts. *
#**************************************************
Motors = int(sys.argv[1])
alpha = int(sys.argv[2])

L = 100
Nb = 3

p = 654
k_walk = 1.
k_on = 0.0125 
k_off = 0.00687

u = 0            # number of unbound motors in the bulk
b = 0            # number of bound motors
N = 0            # the sum of u_x(N_b - b_x)
m0 = 0           #indicates whether the site is occupied by unbound motors
num_of_cargo = 0 # number of cargo currently on MT

ms,B = initialize_motor_distributions(Motors,L)  # Initialized unbound motor dist. and occupied site dist.
cs,C = initialize_cargo_distributions(L)         # Initialized cargo dist., occupied cargo dist, and empty adjacent site dist.
S = initialize_S_distributions()                 # Initialized S dist.

fails_file = "fall_m_"+str(Motors)+"_alpha_"+str(alpha)+".txt"
success_file = "success_m_"+str(Motors)+"_alpha_"+str(alpha)+".txt"
ms_file = "ms_m_"+str(Motors)+"_alpha_"+str(alpha)+".txt"
cs_file = "cs_m_"+str(Motors)+"_alpha_"+str(alpha)+".txt"

file1 = open(fails_file, "w")
file2 = open(success_file, "w")
file3 = open(ms_file, "w")
file4 = open(cs_file, "w")

relax_time = 0
next_meansurement = relax_time
measurements = 100
measurement_interval = 100000
final_time = measurements*measurement_interval + relax_time

step = 0
while step <= final_time:
 step += 1
 if step >= next_meansurement:
  next_meansurement += measurement_interval
  file3.write(' '.join(str(m) for m in ms))
  file3.write('\n')
  file3.flush()
  file4.write(' '.join(str(len(c)) for c in cs))
  file4.write('\n')
  file4.flush()

 N = compute_N()
 D = alpha*m0 + 2.*p*u + p*ms[0] + p*ms[L-1] + k_off*b + k_on*N + k_walk*num_of_cargo
 P = alpha*m0/D  

 R = random.random()
 if R < P:
  # A cargo will step onto the MT.
  add_cargo()
  continue
 if R < P + 2.*p*u/D:
  # An unbound bulk-motor will diffuse.
  diffuse_unbound_bulk_motor()
  continue
 P += 2.*p*u/D
 if R < P + p*ms[0]/D:
  # An unbound motor from the left boundary will diffuse.
  move_left_boundary_motor()
  continue
 P += p*ms[0]/D
 if R < P + p*ms[L-1]/D:
  # An unbound motor from the right boundary will diffuse. 
  move_right_boundary_motor()
  continue
 P += p*ms[L-1]/D
 if R < P + k_off*b/D:
  # A bound motor will unbind.
  val = find_motor_to_unbind()
  if val == -1:
   continue
  file1.write(str(val)+' '+str(step))
  file1.write('\n') 
  file1.flush()    
  continue
 P += k_off*b/D
 if R < P + k_on*N/D:
  # An unbound motor will attach to a cargo.
  find_motor_to_bind()
  continue
 # A cargo will step forward if none of the above options were selected.
 val = find_cargo_to_step()
 if val == 0:
  file2.write(str(step))
  file2.write('\n')
  file2.flush()
  
file1.close()
file2.close()
file3.close()
file4.close()