#!/usr/bin/env python
# coding: utf-8

# # Work for MA30287 CW
# ## Question 1

# In[2]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root
from scipy.optimize import fsolve
import pandas as pd


# In[2]:


# Solve for the ice line
A = 210 # outgoing radiation
B = 1.9 # outgoing radiation
k = 1.6*B # transport parameter
s = lambda y: 1 - 0.482*(3*y**2 - 1)/2 # solar weighting
aw = 0.32 # water albedo
ai = 0.62 # ice albedo
Tc = -10.0 # critical temperature for ice formation
Q0 = 342.0 # solar constant (1370 W/m^2 divided by 4)

Qmin = ((B+k)*(Tc + A/B))/((1-aw)*(s(1)+k/B)) 
Qmax = ((B+k)*(Tc + A/B))/((1-ai)*(s(0)+k/B))
print("Minimal Q for ice-free state = ", Qmin)
print("Maximal Q for ice-covered state = ", Qmax)

abar = lambda ys: ai + (aw - ai)*ys*(1 - 0.241*(ys**2 - 1))
Qfunc = lambda ys: (Tc + A/B)*(B+k)/(s(ys)*(1 - (ai+aw)/2) + k/B*(1 - abar(ys)))
Tbar = lambda ys, Q: (Q*(1 - abar(ys)) - A)/B 
Tbari = lambda Q: (Q*(1 - ai)- A)/B
Tbarw = lambda Q: (Q*(1 - aw)- A)/B
Qfunc = lambda ys: (Tc + A/B)*(B+k)/(s(ys)*(1 - (ai+aw)/2) + k/B*(1 - abar(ys)))


# Creating the ice line
ys = np.linspace(0, 1, 100);

# Seeing how the ice line evolves when the solar constant changes
Qs = Qfunc(ys)

ys = np.linspace(0, 1, 100);
Qs = Qfunc(ys);
plt.plot(Qs, ys, 'k')
plt.plot([Q0, Q0], [0, 1], '--')
plt.plot([Qmin, 550], [1, 1])
plt.plot([250, Qmax], [0, 0])
plt.xlabel('Q');
plt.ylabel('ys');
plt.grid(1)

# We can obtain the ice-line positions via Newton's. 
# At the value of Q0 = 342, to obtain the ice line 
# solve the equation 342 = Qfunc(ys) for the ys. 
fwd = lambda ys: Q0 - Qfunc(ys)
sol = root(fwd, 0.2)
print(sol.message)
ys_low = sol.x[0]
print("lower ice line = ", ys_low)
sol = root(fwd, 0.9)
print(sol.message)
ys_high = sol.x[0]
print("higher ice line = ", ys_high)
plt.plot([Q0, Q0], [ys_low, ys_high], 'o')


# In[3]:


# Plot mean temperature vs. Q
plt.plot(Qs, Tbar(ys, Qs), 'k')
Qa = np.linspace(250, Qmax, 5)
plt.plot(Qa, Tbari(Qa), 'g')
Qb = np.linspace(Qmin, 550, 5)
plt.plot(Qb, Tbarw(Qb), 'orange')
plt.plot([Q0, Q0], [-60, 100], '--')
plt.xlabel('Q')
plt.ylabel('Mean temperature')
plt.grid(1)


# In[4]:


minQ=min(Qs)
print(minQ)

fwd = lambda ys: minQ - Qfunc(ys)
sol = root(fwd, 0.6)
sol = sol.x

A = 210 # outgoing radiation
B = 1.9 # outgoing radiation
k = 1.6*B # transport parameter
s = lambda y: 1 - 0.482*(3*y**2 - 1)/2 # solar weighting
aw = 0.32 # water albedo
ai = 0.62 # ice albedo
Tc = -10.0 # critical temperature for ice formation
Q0 = 342.0 # solar constant (1370 W/m^2 divided by 4)

Qmin = ((B+k)*(Tc + A/B))/((1-aw)*(s(1)+k/B)) 
Qmax = ((B+k)*(Tc + A/B))/((1-ai)*(s(0)+k/B))
print("Minimal Q for ice-free state = ", Qmin)
print("Maximal Q for ice-covered state = ", Qmax)

abar = lambda ys: ai + (aw - ai)*ys*(1 - 0.241*(ys**2 - 1))
Qfunc = lambda ys: (Tc + A/B)*(B+k)/(s(ys)*(1 - (ai+aw)/2) + k/B*(1 - abar(ys)))

# Creating the ice line
ys = np.linspace(0, 1, 100);

# Seeing how the ice line evolves when the solar constant changes
Qs = Qfunc(ys)

unstableys=np.linspace(0,sol,100)
Qs1 = Qfunc(unstableys)
# Plotting the position of the ice line for certain solar constants
plt.plot(Qs1,unstableys,color='black',linestyle='dotted',label='Partially-iced Earth')
plt.xlim([300,450])

stableys=np.linspace(sol,1,100)
Qs2 = Qfunc(stableys)
# Plotting the position of the ice line for certain solar constants
plt.plot(Qs2,stableys,color='black', label='Partially-iced Earth')
plt.xlim([300,450])

# Adding the current solar constant
plt.axvline(Q0, color='orange', linestyle='--', label='Current Solar Constant, $342Wm^{-2}$')

# Adding the constraints for when the Earth is either totally frozen or iceless
plt.axhline(0, xmax=(Qmax-300)/(450-300), color='darkblue', label='Frozen Earth')
plt.axhline(1, xmin=(Qmin-300)/(450-300), color='red', label='Iceless Earth')
#plt.title('Where the ice line is located on the Earth \n depending on how the solar constant Q varies', y=1.03, size=11)
plt.xlabel('Solar Constant, Q')
plt.ylabel('The Ice-line, $y_s$')
plt.legend(bbox_to_anchor =(0.5,-0.6), loc='lower center')

plt.savefig('ys_vs_Q.pdf', bbox_inches='tight')


# In[5]:


# Plot mean temperature vs. Q
plt.plot(Qs1, Tbar(unstableys,Qs1), color='black', linestyle='dotted', label='Partially ice-covered Earth')
plt.plot(Qs2, Tbar(stableys,Qs2), color='black', label='Partially ice-covered Earth')
Qa = np.linspace(250, Qmax, 5)
plt.plot(Qa, Tbari(Qa), 'darkblue', label='Snowball Earth')
Qb = np.linspace(Qmin, 550, 5)
plt.plot(Qb, Tbarw(Qb), 'red', label='Iceless Earth')
plt.plot([Q0, Q0], [-60, 100], color='darkorange',linestyle='--', label='Current Solar Constant')
plt.xlabel('Solar Constant, Q')
plt.ylabel('Mean temperature, $\overline{T}^*$')
plt.legend(loc='upper left')
plt.grid(1)
plt.savefig('Tstar_vs_Q.pdf', bbox_inches='tight')


# In[6]:


# Define parameters
B = 1.9 
Tc = -10
ai = 0.62
aw = 0.32
k = 1.6*B
abar = lambda ys: ai + (aw - ai)*ys*(1 - 0.241*(ys**2 - 1))
s = lambda ys: 1 - 0.482*(3*ys**2 - 1)/2
Qfunc = lambda A: (B+k)*(Tc+A/B)/((1-ai)*(s(0)+k/B))

varyA = np.linspace(180,220,100)
Qvals = Qfunc(varyA)


# In[7]:


minQ = np.zeros(len(varyA))
ys = np.linspace(0,1,100)
for i in range(len(varyA)):
    A=varyA[i]
    Qfunc = lambda ys: (Tc + A/B)*(B+k)/(s(ys)*(1 - (ai+aw)/2) + k/B*(1 - abar(ys)))
    minQ[i]=min(Qfunc(ys))


# In[10]:


plt.plot(varyA, 342-minQ, label = "Earth's Current State \u2192 Frozen Earth")
plt.plot(varyA, Qvals-minQ, label = 'Frozen Earth \u2192 Iceless Earth')
plt.axvline(202, color='black',linestyle = '--', label = 'Current A value')
plt.ylabel('Change in Solar Constant Q \n required for Frozen/Iceless Earth')
plt.xlabel('A')
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.135))
plt.grid(1)
plt.savefig('A_vs_CriticalQ.pdf', bbox_inches='tight')


# ## Question 2

# In[11]:


# Defining parameters 
Q = 342
A = 202
Tc = -10
k = 100
B = 1.9
ai = 0.62
aw = 0.32
ys = 0.939

# Defining the latitude coordiantes
y = np.linspace(0, 1, 100)

# Defining abar and s
abar = lambda y: ai + (aw - ai)*y*(1 - 0.241*(y**2 - 1))
s = lambda y: 1 - 0.482*(3*y**2 - 1)/2

def afun(y, ys): 
    # albedo function; if T > Tc, set a = aw, elseif T < Tc, a = ai
    a = 0*y
    for i, yy in enumerate(y):
        if yy < ys:
            aa = aw
        elif yy > ys:
            aa = ai
        else:
            aa = (ai+aw)/2
        a[i] = aa
    return a

Temp = lambda y: (Q/(B+k))*(s(y)*(1-afun(y,ys))+k*(1-abar(ys))/B)-A/B
esttemp = lambda y: Q*(s(y)*(1-afun(y,ys)))/B-A/B
# equation for when k is v fat
esttemp2 = lambda y: Q*(1-abar(ys))/B-A/B
est = lambda y: Q*(1-abar(ys))/B-A/B+(1/k)*(Q*(s(y)*(1-afun(y,ys)))-Q*(1-abar(ys)))
#esttemp3 = lambda y: Q*(1-abar(ys))/B-A/B + k*(Q*(s(y)*(1-afun(y,ys)))-Q*(1-abar(ys)))
plt.plot(y,Temp(y), label='$T^*(y)$')
plt.plot(y,est(y), label = '$T_0^*(y)+\epsilon T_1^*(y)$')
plt.legend(loc = 'upper right')
plt.xlabel('Latitude, y')
plt.ylabel(u' Temperature (℃)')
plt.grid(1)
plt.savefig('Temp_vs_y_bigk.pdf', bbox_inches='tight')


# In[12]:


# Defining parameters 
Q = 342
A = 202
Tc = -10
k = 0.1
B = 1.9
ai = 0.62
aw = 0.32
ys = 0.939

# Defining the latitude coordiantes
y = np.linspace(0, 1, 100)

# Defining abar and s
abar = lambda y: ai + (aw - ai)*y*(1 - 0.241*(y**2 - 1))
s = lambda y: 1 - 0.482*(3*y**2 - 1)/2

def afun(y, ys): 
    # albedo function; if T > Tc, set a = aw, elseif T < Tc, a = ai
    a = 0*y
    for i, yy in enumerate(y):
        if yy < ys:
            aa = aw
        elif yy > ys:
            aa = ai
        else:
            aa = (ai+aw)/2
        a[i] = aa
    return a

Temp = lambda y: (Q/(B+k))*(s(y)*(1-afun(y,ys))+k*(1-abar(ys))/B)-A/B
esttemp = lambda y: Q*(s(y)*(1-afun(y,ys)))/B-A/B

# Plotting
plt.plot(y,Temp(y), label='$T^*(y)$')
plt.plot(y,esttemp(y), label='$T_0^*(y)$')
plt.legend(loc = 'upper right')
plt.xlabel('Latitude, y')
plt.ylabel(u' Temperature (℃)')
plt.grid(1)
plt.savefig('Temp_vs_y.pdf', bbox_inches='tight')


# ## Question 5

# In[3]:


# Reading in the Average Albedo Temperatures
df = pd.read_csv('AveAlbedo.csv')

# Converting data into useable format
x = df.to_numpy()

# Northern Hemisphere Data
NH = x[0:35]
NH = NH.T
NHLat = NH[0]/90
NHLat = np.flip(NHLat)
NHAvgAlb = np.flip(NH[1])

# Southern Hemisphere data
SH = x[35:len(x)] 
SH = SH.T
SHLat = (-1)*SH[0]/90
SHAvgAlb = SH[1]
print(SHAvgAlb)

x = x.T
x[0]
plt.plot(x[0]/90,x[1])
plt.xlabel('Latitude, y')
plt.ylabel('Average Albedo')
plt.grid(1)
plt.savefig('AvgAlb_vs_y.pdf', bbox_inches='tight')


# In[5]:


s = lambda y: 1 - 0.482*(3*y**2 - 1)/2
vec = []
# Approximating abar using trapezium rule
for i in range(len(NHLat)-1):
    vec.append((NHLat[i+1]-NHLat[i])*(1/2)*(NHAvgAlb[i]*s(NHLat[i])+NHAvgAlb[i+1]*s(NHLat[i+1])))
    
z = sum(vec)

Q = 342
B = 1.9
k = 1.6*B
A = 202
Temp = lambda y: (Q/(B+k))*(s(y)*(1-NHAvgAlb)+k*(1-z)/B)-A/B
Temp(NHLat)
plt.plot(NHLat, Temp(NHLat))


# In[6]:


vec = []
# Approximating abar using riemann sum
for i in range(len(SHLat)-1):
    vec.append((SHLat[i+1]-SHLat[i])*(1/2)*(SHAvgAlb[i]*s(SHLat[i])+SHAvgAlb[i+1]*s(SHLat[i+1])))
    
w = sum(vec)
print(z)
print(w)

Q = 342
B = 1.9
k = 1.6*B
A = 202

meantemp = (Q*(1-(w+z))-A)/B
print(meantemp)

Temp = lambda y: (Q/(B+k))*(s(y)*(1-SHAvgAlb)+k*(1-(w+z))/B)-A/B
Temp(SHLat)
plt.plot(SHLat, Temp(SHLat))


# In[7]:


Temp = lambda y: (Q/(B+k))*(s(y)*(1-NHAvgAlb)+k*(1-(w+z))/B)-A/B
plt.plot(NHLat, Temp(NHLat), label = 'Northern Hemisphere')
Temp = lambda y: (Q/(B+k))*(s(y)*(1-SHAvgAlb)+k*(1-(w+z))/B)-A/B
plt.plot(SHLat, Temp(SHLat), label = 'Southern Hemisphere')
plt.xlabel('Latitude, y')
plt.ylabel(u' Temperature (℃)')
plt.legend(loc='upper right')
plt.grid(1)
plt.savefig('Temp_vs_ERBE.pdf', bbox_inches='tight')


# In[8]:


Temp = lambda y: (Q/(B+2.88*B))*(s(y)*(1-NHAvgAlb)+2.88*B*(1-(w+z))/B)-A/B
plt.plot(NHLat, Temp(NHLat), label = 'Northern Hemisphere')
Temp = lambda y: (Q/(B+3.79*B))*(s(y)*(1-SHAvgAlb)+3.79*B*(1-(w+z))/B)-A/B
plt.plot(SHLat, Temp(SHLat), label = 'Southern Hemisphere')
plt.xlabel('Latitude, y')
plt.ylabel(u' Temperature (℃)')
plt.legend(loc='upper right')
plt.grid(1)


# In[11]:


def afun(y, ys): 
    # albedo function; if T > Tc, set a = aw, elseif T < Tc, a = ai
    a = 0*y
    for i, yy in enumerate(y):
        if yy < ys:
            aa = aw
        elif yy > ys:
            aa = ai
        else:
            aa = (ai+aw)/2
        a[i] = aa
    return a

ys=0.939
aw=0.32
ai=0.62
Q = 342
B = 1.9
yspace=np.linspace(0,1,100)
abar = lambda y: ai + (aw - ai)*y*(1 - 0.241*(y**2 - 1))
k = 1.6*(B1+B2*Ac)
A1 = 257
A2 = -91
B1 = 1.63
B2 = -0.11
Ac = 0.5
Temp = lambda y: (Q*s(y)*(1-afun(y,ys))-(A1+A2*Ac)+14.9)/(B1+B2*Ac+k)
plt.plot(yspace, Temp(yspace), label = 'Northern Hemisphere')
A1 = 262
A2 = -81
B1 = 1.64
B2 = -0.09
k = 1.6*(B1+B2*Ac)
plt.plot(yspace, Temp(yspace), label = 'Southern Hemisphere')
plt.legend()


# In[10]:


A1 = 257
A2 = -91
B1 = 1.63
B2 = -0.11
Ac = 0.5
k = 1.6*(B1+B2*Ac)
Temp = lambda y: (Q/((B1+B2*Ac)+k))*(s(y)*(1-NHAvgAlb)+k*(1-(w+z))/(B1+B2*Ac))-(A1+A2*Ac)/(B1+B2*Ac)
plt.plot(NHLat, Temp(NHLat), label = 'Northern Hemisphere')
A1 = 262
A2 = -81
B1 = 1.64
B2 = -0.09
k = 1.6*(B1+B2*Ac)
Temp = lambda y: (Q/((B1+B2*Ac)+k))*(s(y)*(1-SHAvgAlb)+k*(1-(w+z))/(B1+B2*Ac))-(A1+A2*Ac)/(B1+B2*Ac)
plt.plot(SHLat, Temp(SHLat), label = 'Southern Hemisphere')
plt.xlabel('Latitude, y')
plt.ylabel(u' Temperature (℃)')
plt.yticks(np.arange(-60,40, step=10))
plt.legend(loc='upper right')
plt.grid(1)
plt.savefig('Temp_vs_ERBE2.pdf', bbox_inches='tight')


# In[ ]:




