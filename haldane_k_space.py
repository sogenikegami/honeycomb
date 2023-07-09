import math
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

a = 1 # a:lattice constant
sigmax = np.array([[0,1],[1,0]])
sigmay = np.array([[0,-1j],[1j,0]])
sigmaz = np.array([[1,0],[0,-1]])

a1 = np.array([a,0])
a2 = np.array([-a/2,math.sqrt(3)*a/2])
a3 = np.array([-a/2,-math.sqrt(3)*a/2])
A = [a1,a2,a3]
b1 = np.array([0,math.sqrt(3)*a])
b2 = np.array([-3*a/2,-math.sqrt(3)*a/2])
b3 = np.array([3*a/2,-math.sqrt(3)*a/2])
B = [b1,b2,b3]

t1 = 1
t2 = 1
M = 0.5

path = "/Users/sogenikegami/Documents/UT4S/non-crystal/honeycomb/image/"

def h0(kx,ky,t1):
    H = np.zeros((2,2))*0j
    for i in range(3):
        H += t1*math.cos(kx*A[i][0] + ky*A[i][1])*sigmax
        H -= t1*math.sin(kx*A[i][0] + ky*A[i][1])*sigmay
    return H

def h1(kx,ky,M,t2):
    H = np.zeros((2,2))*0j
    temp = 0
    for i in range(3):
        temp += math.sin(kx*B[i][0]+ky*B[i][1])
    H += sigmaz*(M + 2*t2*temp)
    return H

def hamiltonian(kx,ky,t1,t2,M):
    H = h0(kx,ky,t1) + h1(kx,ky,M,t2)
    return H

def bandplot(t1,t2,M,imgname):
    n = 256
    x = np.linspace(-math.pi/a,math.pi/a,n)
    y = np.linspace(-math.pi/a,math.pi/a,n)
    X,Y = np.meshgrid(x,y)
    z1 = np.zeros((n,n))
    z2 = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            kx = x[i]
            ky = y[j]
            eigval,eigvec = np.linalg.eigh(hamiltonian(kx,ky,t1,t2,M))
            energy = [eigval[0].real,eigval[1].real]
            z1[i][j] += max(energy)
            z2[i][j] += min(energy)


    fig = plt.figure(figsize = (8, 8))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(X, Y, z1, cmap = "plasma_r")
    ax.plot_surface(X, Y, z2, cmap = "plasma_r")
    fig.show()
    plt.savefig(path + imgname)
    



def main():
    imgname = "test"
    bandplot(t1,t2,M,imgname)
    print(A[0][0])
    return 0
if __name__ == "__main__":
    main()