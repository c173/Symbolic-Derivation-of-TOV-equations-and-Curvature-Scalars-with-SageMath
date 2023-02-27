#!/usr/bin/env python
# coding: utf-8

# In[1]:


Man = Manifold(4, 'M') #Manifold
X.<t,r,th,ph> = Man.chart(r't r:(0,+oo) th:(0,pi):\theta ph:(0,2*pi):\phi')

rho = function('rho', latex_name='\\rho')
P = function('p')
nu = function('nu')
m = function('m')

g = Man.lorentzian_metric('g')
g[0,0] = -exp(2*nu(r))
g[1,1] = 1/(1-2*m(r)/r)
g[2,2] = r^2
g[3,3] = (r*sin(th))^2

print('Metric:')
show(g.display())



Ricci_scalar = g.ricci_scalar() #Ricci Scalar.
Ricci = g.ricci() #Ricci tensor.
upRic = Ricci.up(g) #Index-uppered Ricci tensor.
gg= g.up(g) #Index-uppered metric.

ET = Ricci - 1/2*g.ricci_scalar() * g
ET.set_name('G')
ETup=ET.up(g)
ETmix=ETup['^ij']*g['_jk']
ETmix.set_name('G')

riem = g.riemann() #Riemann tensor and its lowered-uppered index ver.
riemup = riem.up(g)
riemdown = riem.down(g)

weyl = g.weyl() #Weyl tensor and its lowered-uppered index ver.
weylup = weyl.up(g)
weyldown= weyl.down(g)

Weyl_scalar = weyldown['_ijkl']*weylup['^ijkl']
Contracted_ricci = Ricci['_ij']*upRic['^ij']
Kretschmann_scalar = riemdown['_ijkl']*riemup['^ijkl']

show('EFEs:')
show(ET.display_comp())

show('EFEs with Mixed Components:')
show(ETmix.display_comp())


# In[2]:


#Defining the 4-velocity
u = Man.vector_field('u')
u[0] = exp(-nu(r))

print('The 4-velocity:')
show(u.display())

u_form = u.down(g) #index lowered 4-velocity.


#Stress-Energy Tensor
print('Stress-Energy tensor:')
T = (rho(r)+p(r))* (u_form * u_form) + p(r) * g
T.set_name('T')
show(T.display_comp())
Tmix = g.inverse()['^ij']*T['_jk']


#defining the E=G-kT=0 equation to solve the Einstein's equations.
print('E=G-kT=0 equation:')
E = ET - 8*pi*T
E.set_name('E')
show(E.display_comp())


# In[3]:


#solving the m', nu' ve nu''
EE0_solv = solve(E[0,0].expr()==0, diff(m(r),r))
EE_0 = EE0_solv[0]
show(EE_0.full_simplify())

EE1_solv = solve(E[1,1].expr()==0, diff(nu(r),r))
EE_1 = EE1_solv[0]
show(EE_1.full_simplify())

EE2_solv = solve(E[2,2].expr()==0, diff(diff(nu(r),r),r))
EE_2 = EE2_solv[0]
show(EE_2.full_simplify())


#I redefined derivative expressions to replace them in Curvature scalars.
diffm=EE_0.rhs()
diffnu=EE_1.rhs()
diffdiffnu=(8*pi*r^3*p(r) - (r^3 - 2*r^2*m(r))*((4*pi*(r**3)*p(r)+m(r))/((r**2)-2*r*m(r)))^2 + r*(4*pi*(r**2)*rho(r)) + (r^2*(4*pi*(r**2)*rho(r)) - r^2 + r*m(r))*((4*pi*(r**3)*p(r)+m(r))/((r**2)-2*r*m(r))) - m(r))/(r^3 - 2*r^2*m(r))

#diffdiffnu=EE_2.rhs()
#((4*pi*(r**3)*p(r)+m(r))/((r**2)-2*r*m(r)))
#(4*pi*(r**2)*rho(r))


# In[4]:


#d(nu)/dr d^2(nu)/dr^2 ve d(m)/dr'yi is replaced in terms of p and rho in the curvature scalars. -This is done by hand-

print('Ricci scalar in terms of Energy-Momentum tensor')
naked_Ricci_scalar=-2*((r^2 - 2*r*m(r))*((diffnu))^2 - (r*(diffm) - 2*r + 3*m(r))*(diffnu) + (r^2 - 2*r*m(r))*(diffdiffnu) - 2*(diffm))/r^2
show(naked_Ricci_scalar.full_simplify())

print('Fully Contracted Ricci Tensor in terms of Energy-Momentum tensor')
naked_Contracted_ricci=2*((r^6 - 4*r^5*m(r) + 4*r^4*m(r)^2)*((diffnu))^4 + 3*r^2*((diffm))^2 + 2*(r^5 - 3*r^4*m(r) + 2*r^3*m(r)^2 - (r^5 - 2*r^4*m(r))*(diffm))*((diffnu))^3 - 2*r*m(r)*(diffm) + (r^4*((diffm))^2 + 3*r^4 - 8*r^3*m(r) + 5*r^2*m(r)^2 - 2*(2*r^4 - 3*r^3*m(r))*(diffm))*((diffnu))^2 + (r^6 - 4*r^5*m(r) + 4*r^4*m(r)^2)*((diffdiffnu))^2 + 3*m(r)^2 + 2*(r^3*((diffm))^2 - r^3*(diffm) - r^2*m(r) + 3*r*m(r)^2)*(diffnu) + 2*(r^3*m(r) - 2*r^2*m(r)^2 + (r^6 - 4*r^5*m(r) + 4*r^4*m(r)^2)*((diffnu))^2 - (r^4 - 2*r^3*m(r))*(diffm) + (r^5 - 3*r^4*m(r) + 2*r^3*m(r)^2 - (r^5 - 2*r^4*m(r))*(diffm))*(diffnu))*(diffdiffnu))/r^6
show(naked_Contracted_ricci.full_simplify())

print('Kretschmann scalar in terms of Energy-Momentum tensor')
naked_Kretschmann_scalar=4*((r^6 - 4*r^5*m(r) + 4*r^4*m(r)^2)*((diffnu))^4 + 2*r^2*((diffm))^2 + 2*(r^4*m(r) - 2*r^3*m(r)^2 - (r^5 - 2*r^4*m(r))*(diffm))*((diffnu))^3 - 4*r*m(r)*(diffm) + (r^4*((diffm))^2 - 2*r^3*m(r)*(diffm) + 2*r^4 - 8*r^3*m(r) + 9*r^2*m(r)^2)*((diffnu))^2 + (r^6 - 4*r^5*m(r) + 4*r^4*m(r)^2)*((diffdiffnu))^2 + 6*m(r)^2 + 2*((r^6 - 4*r^5*m(r) + 4*r^4*m(r)^2)*((diffnu))^2 + (r^4*m(r) - 2*r^3*m(r)^2 - (r^5 - 2*r^4*m(r))*(diffm))*(diffnu))*(diffdiffnu))/r^6
show(naked_Kretschmann_scalar.full_simplify().expand())

print('Weyl scalar in terms of Energy-Momentum tensor')
naked_Weyl_scalar=4/3*((r^6 - 4*r^5*m(r) + 4*r^4*m(r)^2)*((diffnu))^4 + r^2*((diffm))^2 - 2*(r^5 - 5*r^4*m(r) + 6*r^3*m(r)^2 + (r^5 - 2*r^4*m(r))*(diffm))*((diffnu))^3 - 6*r*m(r)*(diffm) + (r^4*((diffm))^2 + r^4 - 12*r^3*m(r) + 21*r^2*m(r)^2 + 2*(2*r^4 - 5*r^3*m(r))*(diffm))*((diffnu))^2 + (r^6 - 4*r^5*m(r) + 4*r^4*m(r)^2)*((diffdiffnu))^2 + 9*m(r)^2 - 2*(r^3*((diffm))^2 - 3*r^2*m(r) + 9*r*m(r)^2 + (r^3 - 6*r^2*m(r))*(diffm))*(diffnu) - 2*(3*r^3*m(r) - 6*r^2*m(r)^2 - (r^6 - 4*r^5*m(r) + 4*r^4*m(r)^2)*((diffnu))^2 - (r^4 - 2*r^3*m(r))*(diffm) + (r^5 - 5*r^4*m(r) + 6*r^3*m(r)^2 + (r^5 - 2*r^4*m(r))*(diffm))*(diffnu))*(diffdiffnu))/r^6
show(naked_Weyl_scalar.full_simplify().expand())


# In[5]:


##TOV Equation

Tmix=g.inverse()['^ij']*T['_jk']
nabla=g.connection() #Levi-Civita connection

dTmix=nabla(Tmix) #Covariance derivative form.
divT = dTmix['^i_{ji}'] #Trace is taken to find the Divergence.
show(divT.display())

#Equated the divergence of T to zero to find the necessary equations for TOV equations.
EEp_sol = solve(divT[1].expr()==0, diff(p(r),r))
EEp = EEp_sol[0]
show(EEp.full_simplify())


# In[6]:


#TOV equations using a simple polytropic EoS as an example.
k=var('k')
p_eos(r)=k*rho(r)^2
print('for a basic EoS')
show(p,'=',p_eos(r))


#I substituted p_eos and nu to the TOV equations I solved.
EE_1_rho = EE_1.substitute_function(p, p_eos)

EEp_rho = EEp.substitute_function(p, p_eos)
EEp_rho = (EEp_rho / (2*k*rho(r))).simplify_full()
EEp_rho = EEp_rho.subs({diff(nu(r),r): EE_1_rho.rhs()}).simplify_full()


print('TOV equations')
#TOV equations
for i in [EE_0, EE_1_rho, EEp_rho]:
    show(i)


# In[42]:


test2=2*(6*r^5 + (r^6 - 4*r^5*m(r) + 4*r^4*m(r)^2)*((diffnu))^4 + 6*r^4 + 3*r^2*((diffm))^2 + 2*(r^4*m(r) - 2*r^3*m(r)^2 - (r^5 - 2*r^4*m(r))*(diffm))*((diffnu))^3 + 3*(16*r^4 + 3)*m(r)^2 - (4*r^8 + 4*r^7 - r^4*((diffm))^2 - 5*r^4 + (16*r^6 - 13*r^2)*m(r)^2 - 8*(2*r^7 + r^6 - 2*r^3)*m(r) - 2*(r^4 - 3*r^3*m(r))*(diffm))*((diffnu))^2 + (r^6 - 4*r^5*m(r) + 4*r^4*m(r)^2)*((diffdiffnu))^2 - 12*(2*r^5 + 3*r^4)*m(r) - 4*(r^6 + r^5 - (2*r^5 + r)*m(r))*(diffm) + 2*(2*r^7 + 2*r^6 - r^3*((diffm))^2 + 12*(r^5 + r)*m(r)^2 - (10*r^6 + 6*r^5 + 5*r^2)*m(r) + (2*r^7 + 2*r^6 - (4*r^6 + r^2)*m(r))*(diffm))*(diffnu) - 2*(2*r^8 + 2*r^7 + 4*(2*r^6 + r^2)*m(r)^2 - (r^6 - 4*r^5*m(r) + 4*r^4*m(r)^2)*((diffnu))^2 - 2*(4*r^7 + 2*r^6 + r^3)*m(r) - (r^4 - 2*r^3*m(r))*(diffm) - (r^4*m(r) - 2*r^3*m(r)^2 - (r^5 - 2*r^4*m(r))*(diffm))*(diffnu))*(diffdiffnu))/r^6


# In[47]:


show(Tmix.display_comp())


# In[103]:


k=var('k')

a=g[1,1]
b=a^(-1)
test=2*Tmex[0,0]*Tmex[1,1]-2*Tmex[0,0]*Tmex[2,2]-2*Tmex[0,0]*Tmex[3,3]-2*Tmex[2,2]*Tmex[1,1]-2*Tmex[2,2]*Tmex[3,3]+3*(Tmex[0,0])^2+3*(Tmex[1,1])^2+2*(Tmex[2,2])^2+2*(Tmex[3,3])^2+8*(Tmex[0,0]+Tmex[1,1]-(1/2)*(Tmex[2,2]+Tmex[3,3]))*(1-b)*r^2+12*b^2*(1-a)^2
test.display()
damn=16*(12*pi^2*r^2*p(r)^2 + 12*pi^2*r^2*rho(r)^2 + 3*m(r)^2 - 8*(pi*r^3*m(r) - pi^2*r^2*p(r))*rho(r))/r^2


# In[62]:


ett=2*ETmix[0,0]*ETmix[1,1]-2*ETmix[0,0]*ETmix[2,2]-2*ETmix[0,0]*ETmix[3,3]-2*ETmix[2,2]*ETmix[1,1]-2*ETmix[2,2]*ETmix[3,3]+3*(ETmix[0,0])^2+3*(ETmix[1,1])^2+2*(ETmix[2,2])^2+2*(ETmix[3,3])^2+8*(ETmix[0,0]+ETmix[1,1]-(1/2)*(ETmix[2,2]+ETmix[3,3]))*(1-b)*r^2+12*b^2*(1-a)^2
ett.display()


# In[104]:


esd=2*((r^6 - 4*r^5*m(r) + 4*r^4*m(r)^2)*((diffnu))^4 + 3*r^2*((diffm))^2 + 2*(r^4*m(r) - 2*r^3*m(r)^2 - (r^5 - 2*r^4*m(r))*(diffm))*((diffnu))^3 - 4*(2*r^5 - r)*m(r)*(diffm) + (r^4*((diffm))^2 + 5*r^4 + (16*r^6 + 13*r^2)*m(r)^2 - 8*(r^7 + 2*r^3)*m(r) + 2*(r^4 - 3*r^3*m(r))*(diffm))*((diffnu))^2 + (r^6 - 4*r^5*m(r) + 4*r^4*m(r)^2)*((diffdiffnu))^2 + 9*m(r)^2 - 2*(r^3*((diffm))^2 + 12*(r^5 - r)*m(r)^2 - (4*r^6 - r^2)*m(r)*(diffm) - (4*r^6 - 5*r^2)*m(r))*(diffnu) + 2*(4*(2*r^6 - r^2)*m(r)^2 + (r^6 - 4*r^5*m(r) + 4*r^4*m(r)^2)*((diffnu))^2 - 2*(2*r^7 - r^3)*m(r) + (r^4 - 2*r^3*m(r))*(diffm) + (r^4*m(r) - 2*r^3*m(r)^2 - (r^5 - 2*r^4*m(r))*(diffm))*(diffnu))*(diffdiffnu))/r^6
show(esd.full_simplify().expand())
show(damn.full_simplify().expand())


# In[ ]:




