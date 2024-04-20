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
     breakx
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
 global b,num_of_cargo,m0
 motors_to_bind = 1 if np.minimum(ms[0],Nb) == 1 else random.randint(1,np.minimum(ms[0],Nb))
 ms[0] -= motors_to_bind
 cs[0].append(motors_to_bind)
 num_of_cargo += 1            # Incrementing the count for motors on the MT.
 b += motors_to_bind          # Updating the number of motors "b" bound to cargos along the MT.
 #***************************************************************************************************
 # If C was empty. There are a few possibilities for what to update, aside from adding 0 to C:      *
 # (1) We add S if there are unbound motors and  motors_to_bind < Nb.                               *
 # (2) We update m0 = 0 if there are zero unbound motors.                                           *
 # (3) Otherwise, only add C because motors_to_bind = Nb and potentially update m0 depending on the *
 # number of unbound motors at site = 0.                                                            *
 #***************************************************************************************************
 if 0 not in C:
  C.append(0)
  if ms[0] > 0 and motors_to_bind < Nb:
   S.append(0)
   return
  if ms[0] == 0:
   m0 = 0    
  return
 #*********************************************************************
 # If cs[0] is not empty then we may potentially need to update S.    *
 # There are a few possibilities for what to update:                  *
 # (1) We add S if there are unbound motors and  motors_to_bind < Nb. *
 # (2) Otherwise if ms[0] == 0 and 0 in S, remove 0 from S.           *
 #*********************************************************************
 if ms[0] > 0 and motors_to_bind < Nb and 0 not in S:
  S.append(0)
  return
 if ms[0] == 0: 
  m0 = 0
  if 0 in S:
   S.remove(0)
 return
def move_bulk_motor_left(site):  
 l_site = site - 1
 ms[site] -= 1
 ms[l_site] += 1
 #**********************************************************************************************
 # (1)  If ms[site] == 0 then we have to remove "site" from B and potentially S.               *
 # (2) If l_site ia not in S and there are binding sites available at l_site, add l_site to S. *
 # (3)  If l_site = 0 then update u -= 1 and set m0 = 1.                                       *
 #**********************************************************************************************
 if ms[site] == 0:
  B.remove(site)
  if site in S:
   S.remove(site)
 if l_site not in S and sum(cs[l_site]) < len(cs[l_site])*Nb:
  S.append(l_site)
 #*****************************************************************
 # The conditions above should be applied to any site and l_site. *
 # Now for updates for particular sites.                          *
 #*****************************************************************
 if l_site > 0:
  if l_site not in B:
   B.append(l_site)
 else:
  global u,m0
  u -= 1
  m0 = 1
def move_bulk_motor_right(site): 
 r_site = site + 1
 ms[site] -= 1
 ms[r_site] += 1
 #***********************************************************************************************
 # (1)  If ms[site] == 0 then we have to remove "site" from B and potentially S.                *
 # (2)  If r_site ia not in S and there are binding sites available at r_site, add r_site to S. *
 # (3)  If r_site < L-1 then update u -= 1.                                                     *
 #***********************************************************************************************
 if ms[site] == 0:
  B.remove(site)
  if site in S:    
   S.remove(site)
 if r_site not in S and sum(cs[r_site]) < len(cs[r_site])*Nb:
  S.append(r_site)
 #*****************************************************************
 # The conditions above should be applied to any site and r_site. *
 # Now for updates for particular sites.                          *
 #*****************************************************************
 if r_site < L-1:
  if r_site not in B:
   B.append(r_site)
 else:
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
  if 1 not in S and sum(cs[1]) < len(cs[1])*Nb:
   S.append(1)
 else:
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
  if l_site not in S and sum(cs[l_site]) < len(cs[l_site])*Nb:
   S.append(l_site)
 else:
  B.append(l_site) 
def unbind_motor(site,indx):   
 global b
 b -= 1
 ms[site] += 1
 cs[site][indx] -= 1
 len_cs = len(cs[site])
 if site > 0:
  if site < L-1:
   global u
   u += 1
   if site not in B:
    B.append(site)
 else:
  global m0
  m0 = 1
 #*************************************
 # The following applies to any site. *
 #*************************************
 if cs[site][indx] == 0:
  global num_of_cargo    
  num_of_cargo -= 1
  cs[site].pop(indx);
  if (len_cs - 1) == 0:
   C.remove(site)
   if site in S:
    S.remove(site)
  elif site in S and sum(cs[site]) == (len_cs - 1)*Nb:
   S.remove(site)
  return site if step >= relax_time else -1
 elif site not in S:
  S.append(site)
 return -1
def bind_motor(site,indx):
 global b
 b += 1
 ms[site] -= 1
 cs[site][indx] += 1
 #*******************************************************************************
 # We know site is already in S. Now we have to decide whether to remove or not *
 # to remove site from S. If there are no more binding sites on cs[site][indx]  *
 # or if there are no more unbound motors at site, then we remove site from S.  *
 #*******************************************************************************
 if ms[site] == 0 or sum(cs[site]) == len(cs[site])*Nb:
  S.remove(site)
 if site > 0 and site < L-1:
  global u
  u -= 1
  if ms[site] == 0:
   B.remove(site)  
  return 
 elif site == 0 and ms[site] == 0:
  global m0
  m0 = 0
def move_a_cargo(site):
 len_cs = len(cs[site])
 indx = random.choice(range(len_cs))
 if site < L-1:
  #******************************
  # Updating the adjacent site. *
  #******************************
  r_site = site + 1
  if r_site not in S and cs[site][indx] < Nb and ms[r_site] > 0:
   S.append(r_site)   
  if r_site not in C: 
   C.append(r_site)
  cs[r_site].append(cs[site][indx])
  #*****************
  # Updating site. *
  #*****************
  cs[site].pop(indx)
  if (len_cs-1) == 0:
   C.remove(site)
   if site in S:
    S.remove(site)
   return -1
  #************************************************************************************
  # The cargo that stepped may have been the only cargo with binding sites available. *
  # For this reason we are performing the following update.                           *
  #************************************************************************************
  if site in S and sum(cs[site]) == (len_cs-1)*Nb:
   S.remove(site) 
  return -1
 else:
  global b,num_of_cargo
  num_of_cargo -= 1 
  b -= cs[site][indx]
  ms[site] += cs[site][indx]
  cs[site].pop(indx)
  if (len_cs-1) == 0:
   C.remove(site)
   if site in S:
    S.remove(site)
   return 0 if step >= relax_time else -1
  #************************************************************************************
  # The cargo that stepped may have been the only cargo with binding sites available. *
  # For this reason we are performing the following update.                           *
  #************************************************************************************
  if site in S and sum(cs[site]) == (len_cs-1)*Nb:
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
    temp = ms[site]*(Nb - val)
    if const1 < const2 + temp: 
     bind_motor(site,indx)
     return
    const2 += temp
  const2 += temp  
def find_cargo_to_step():
 const1 = R*D
 const2 = P*D
 for site in C:
  temp = len(cs[site])
  if const1 < const2 + temp:
   pass_val = move_a_cargo(site)
   return pass_val
  const2 += temp
#**************************************************  
# This is where the main part of the code starts. *
#**************************************************
list_of_alphas=[0.1,1,10,100,1000]

Motors = int(sys.argv[1])
alpha_indx = int(sys.argv[2])
alpha = list_of_alphas[alpha_indx]

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
bs_file = "bs_m_"+str(Motors)+"_alpha_"+str(alpha)+".txt"

file1 = open(fails_file, "w")
file2 = open(success_file, "w")
file3 = open(ms_file, "w")
file4 = open(cs_file, "w")
file5 = open(bs_file, "w")

relax_time = 0
next_meansurement = relax_time
measurements = 100000
measurement_interval = 100000
final_time = measurements*measurement_interval + relax_time

step = 0
while step <= final_time:
 step += 1
 if step >= next_meansurement:
  next_meansurement += measurement_interval
  #*********************************************
  file3.write(' '.join(str(m) for m in ms))
  file3.write('\n')
  file3.flush()
  #*********************************************
  file4.write(' '.join(str(len(c)) for c in cs))
  file4.write('\n')
  file4.flush()
  #*********************************************
  file5.write(' '.join(str(sum(c)) for c in cs))
  file5.write('\n')
  file5.flush()
 #*************************************************************************************
 N = compute_N()
 D = alpha*m0 + 2.*p*u + p*ms[0] + p*ms[L-1] + k_off*b + k_on*N + num_of_cargo
 R = random.random()
 #*************************************************************************************
 P = alpha*m0/D  
 if R < P:
  add_cargo()
  continue
 if R < P + 2.*p*u/D:
  diffuse_unbound_bulk_motor()
  continue
 P += 2.*p*u/D
 if R < P + p*ms[0]/D:
  move_left_boundary_motor()
  continue
 P += p*ms[0]/D
 if R < P + p*ms[L-1]/D:
  move_right_boundary_motor()
  continue
 P += p*ms[L-1]/D
 if R < P + k_off*b/D:
  val = find_motor_to_unbind()
  if val == -1:
   continue
  file1.write(str(val)+' '+str(step))
  file1.write('\n') 
  file1.flush()    
  continue
 P += k_off*b/D
 if R < P + k_on*N/D:
  find_motor_to_bind()
  continue
 val = find_cargo_to_step()
 if val == 0:
  file2.write(str(step))
  file2.write('\n')
  file2.flush()
  
file1.close()
file2.close()
file3.close()
file4.close()
file5.close()
