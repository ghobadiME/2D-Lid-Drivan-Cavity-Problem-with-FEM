#2D lid driven cavity FDM
import numpy as np
import matplotlib.pyplot as plt
Re=float(input('please enter Reynolds number='))
if (Re>99):
   dx=0.01
   dy=0.01
elif (9<Re<99) :
   dx=0.02
   dy=0.02  
else :
  dx=0.05
  dy=0.05
dt = min((dx**2)*Re*0.25,4.0/Re)
tmax=10.0
U_up=1.0
nu = 1.0/Re
nx = int(1.0/dx)
ny=int(1.0/dy)
nt=int(tmax/dt)
xmax=1.0
ymax=1.0
#initial condition
omega=np.zeros([nx, ny])
omega_1=np.zeros([nx, ny])
psi=np.zeros([nx, ny])
psi_1=np.zeros([nx, ny])
u=np.zeros([nx,ny])
v=np.zeros([nx,ny])
P=np.ones([nx,ny])
#boundary condition
psi[0:nx,ny-1]=0.5*dy
psi[0:nx, 0] = 0.0
psi[0,0:ny]=0.0
psi[nx-1,0:ny]=0.0
u[0:nx,0]=0.0
u[0:nx,ny-1]=1.0
u[0,0:ny]=0.0
u[nx-1,0:ny]=0.0
v[0:nx,0]=0.0
v[0:nx,ny-1]=0.0
v[0,0:ny]=0.0
v[nx-1,0:ny]=0.0
# solve 
for n in range (nt):
  for i in range(nx-2):
   for j in range(ny-2):
      psi[i+1,j+1]=0.25*(psi[i+2,j+1]+psi[i,j+1]+psi[i+1,j+2]+psi[i+1,j]+(dx**2)*omega[i+1,j+1])
      omega[0:nx,0] = (2.0/dy**2)*(psi[0:nx,0]-psi[0:nx,1])
      omega[0:nx,ny-1] = (2.0/(dy**2))*(psi[0:nx,ny-1]-psi[0:nx,ny-2])+2.0*U_up/dy
      omega[0, 0:ny] = (2.0/dx**2)*(psi[0,0:ny]-psi[1,0:ny])
      omega[nx-1,0:ny] = (2.0/dx**2)*(psi[nx-1,0:ny]-psi[nx-2,0:ny]) 
      omega_1[i+1,j+1]=  (nu* ((omega[i+2,j+1] - 2.0* omega[i+1,j+1] + omega[i,j+1]) / dx**2)+nu*((omega[i+1,j+2]-2.0*omega[i+1,j+1]+omega[i+1,j])/dy**2)+((psi[i+2,j+1]-psi[i,j+1])/(2.0*dx))*((omega[i+1,j+2]-omega[i+1,j])/(2.0*dy))-((psi[i+1,j+2]-psi[i+1,j])/(2.0*dy))*((omega[i+2,j+1]-omega[i,j+1])/(2.0*dx)))   
      omega[i+1,j+1]=omega[i+1,j+1]+dt*omega_1[i+1,j+1] 
for i in range(nx-2):
 for j in range(ny-2):
      u[i+1,j+1]=(psi[i+1,j+1]-psi[i+1,j])/dy
      v[i+1,j+1]=(psi[i,j+1]-psi[i+1,j+1])/dx 
for n in range (nt):  
 for o in range(nx-2):
  for p in range(ny-2):
      P[o+1,p+1]=0.25*(P[o+2,p+1]+P[o,p+1]+P[o+1,p+2]+P[o+1,p])-(dx/16.0)*((2.0/dt)*(u[o+2,p+1]-u[o,p+1]+v[o+1,p+2]-v[o+1,p])-(2.0/dx)*(u[o+1,p+2]-u[o+1,p])*(v[o+2,p+1]-v[o,p+1])-(1.0/dx)*((u[o+2,p+1]-u[o,p+1])**2.0)-(1.0/dy)*((v[o+1,p+2]-v[o+1,p])**2.0))
#plot results  
x = np.linspace(0.0, xmax,nx)
y = np.linspace(0.0, ymax,ny)
X,Y= np.meshgrid(x,y,indexing='ij')
z1=omega[0:nx,0:ny]
z2=psi[0:nx,0:ny]
z3=u[0:nx,0:ny]
z4=v[0:nx,0:ny]
z5=P[0:nx,0:ny]
fig,ax1=plt.subplots(1,1)
cp1 = ax1.contourf(X,Y,z1)
fig.colorbar(cp1)
fig,ax2=plt.subplots(1,1)
cp2 = ax2.contourf(X,Y,z2)
fig.colorbar(cp2)
fig,ax3=plt.subplots(1,1)
cp3 = ax3.contourf(X,Y,z3)
fig.colorbar(cp3)
fig,ax4=plt.subplots(1,1)
cp4 = ax4.contourf(X,Y,z4)
fig.colorbar(cp4)
fig,ax5=plt.subplots(1,1)
cp5 = ax5.contourf(X,Y,z5)
fig.colorbar(cp5)
ax1.set_title('omega, Reynolds number='+str(Re))
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax2.set_title('psi, Reynolds number='+str(Re))
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax3.set_title('u, Reynolds number='+str(Re))
ax3.set_xlabel('x')
ax3.set_ylabel('y')
ax4.set_title('v, Reynolds number='+str(Re))
ax4.set_xlabel('x')
ax4.set_ylabel('y')
ax5.set_title('P, Reynolds number='+str(Re))
ax5.set_xlabel('x')
ax5.set_ylabel('y')
plt.show()
