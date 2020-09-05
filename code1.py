import numpy as np
import matplotlib.pyplot as plt
from numpy import sin, cos, sqrt
from scipy.integrate import quad
from mpl_toolkits.mplot3d import axes3d
from scipy.integrate import dblquad

pi = 3.14159
mu = 4*pi*10**-7
M0 = 1.21*1000/mu # Magnetization, equal to residual magnetic flux divided by mu
b = 6.35/1000 # radius of magnet
L = 12.7/1000

# quad returns a tuple
# First value is the estimated integral, second one is upper bound on error

def By(y, z):
    return 1.21*1000*b/(4*pi)*(dblquad(lambda phi, zprime: sin(phi)*(z-zprime)/((b**2+y**2+(z-zprime)**2-2*y*b*sin(phi))**1.5),
               -L, 0, lambda phi: 0, lambda phi: 2*pi)[0])
    
y_b = np.arange(0, 3.1, 0.1)
ByVector = [By(y_*b, b*1)/(1.21*1000*b/(4*pi)) for y_ in y_b]
plt.plot(y_b, ByVector)



# Visualize 3d mesh plot of function
pts1, pts2 = 30, 30

m = np.linspace(-b, 4*b, pts1)
n = np.linspace(0, 0.3*b, pts2)

M, N = np.meshgrid(m, n)
O = np.zeros((pts1, pts2))
for i in range(pts1):
    for j in range(pts2):
        O[i][j] = By(M[i][j], N[i][j])
        
fig = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(M, N, O, 50, cmap='binary')



#--------- Do not run these two codes together-------------- #

# Visualize the contour plot of function
plt.axis("equal")
plt.contour(M, N, O)






