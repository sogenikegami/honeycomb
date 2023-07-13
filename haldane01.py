import matplotlib.pyplot as plt
import math
import numpy as np
import honeycomb_lattice
import sys
sys.path.append('..')
from crosshair import crosshair_marker01
from local_chern_marker import local_chern_marker01


a = 1#lattice constant
n = 15 #systemsize x axis
m = 30 #systemsize y axis

#boundary condition 
xPBC = False
yPBC = True

pointdata = honeycomb_lattice.vertex_create(a,n,m,xPBC,yPBC)
#hamiltonian parameter
M = 1 #haldane parameter site energy
t1 = 1.0 #haldane nn hopping parameter
t2 = 1/3 + 0j #haldae nnn hopping parameter
phi = math.pi/3 #local flux mod pi
filepath = "/Users/sogenikegami/Documents/UT4S/non-crystal/honeycomb/image/"

#fermi_energy = 0.1
Rx = 5 #crosshair position
Ry = 8 #crosshair position




def hsite(pointdata,M):
    vertexdata = pointdata[0]
    #plaquettedata = pointdata[1]
    N = len(vertexdata)
    H = np.zeros((N,N))*0j
    for i in range(N):
        if i % 2 == 0:
            H[i][i] += (-1)*M
        else:
            H[i][i] += M
    return H

def hhop1(pointdata,t1):
    vertexdata = pointdata[0]
    #plaquettedata = pointdata[1]
    N = len(vertexdata)
    H = np.zeros((N,N))*0j
    for i in range(N):
        neighbor = vertexdata[i]["neighbor"]
        for j in range(len(neighbor)):
            H[neighbor[j]][i] += t1
    return H

def hhop2(pointdata,t2,phi):
    vertexdata = pointdata[0]
    plaquettedata = pointdata[1]
    N = len(vertexdata)
    H = np.zeros((N,N))*0j
    for i in range(N):
        nnn = vertexdata[i]["secondneighbor"]
        for j in range(len(nnn)):
            eiphi = math.cos(phi) + 1j * math.sin(phi) * nnn[j][1]
            H[nnn[j][0]][i] += t2 * eiphi
    return H

def hamiltonian(pointdata,M,t1,t2,phi):
    H = hsite(pointdata,M) + hhop1(pointdata,t1) + hhop2(pointdata,t2,phi)
    return H

def fourier_y(vertexdata,H,ky):
    N = len(vertexdata)

    h = np.zeros((2*(n+1),2*(n+1)))*0j
    for i in range(2*(n+1)):
        for j in range(2*(n+1)):
            for k in range(m+1):
                for l in range(m+1):
                    yk = vertexdata[i+2*(n+1)*k]["pos"][1]
                    yl = vertexdata[j+2*(n+1)*l]["pos"][1]
                    phase = math.cos(ky*yl-ky*yk) + 1j*math.sin(ky*yl-ky*yk)
                    h[i][j] += phase/2/math.pi * H[i+2*(n+1)*k][j+l*2*(n+1)]
    eigval,eigvec = np.linalg.eig(h)
    energy = []
    for i in range(len(eigval)):
        energy.append(eigval[i].real)
    energy.sort()
    """
    fig=plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.scatter(ky,energy)
    plt.savefig(filepath + imgname)
    fig.show()
    """
    return energy

def y_band_plot(vertexdata,H,imgname):
    L = 2*math.pi/math.sqrt(3)/a/(m+1)
    K = np.linspace(0,L*(m+1),m+1)
    E = np.zeros((2*n+2,m+1))
    fig=plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i in range(m+1):
        energy = fourier_y(vertexdata,H,i*L,)
        for j in range(len(energy)):
            E[j][i] += energy[j]
    for i in range(2*n+2):
        plt.plot(K,E[i],color='black')
    plt.savefig(filepath + imgname)
    fig.show()
    return 0


def eigenenergy(H,filepath,imgname):
    eigval,eigvec = np.linalg.eigh(H)
    num = list(range(0,len(eigval)))
    eigenergy = []
    for i in range(len(eigval)):
        eigenergy.append(eigval[i].real)
    eigenergy.sort()
    fermi_energy = (eigenergy[int(len(H[0])/2)]+eigenergy[int(len(H[0])/2)-1])/2
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.scatter(num,eigenergy,s=10)
    plt.hlines(fermi_energy,0,len(eigval),linestyle="dashed",color='red')
    fig.show()
    plt.savefig(filepath + imgname)
    
    return fermi_energy

def main():
    H = hamiltonian(pointdata,M,t1,t2,phi)
    imgname = "haldane_eigval" + str(n) + "times" + str(m)
    fermi_energy = eigenenergy(H,filepath,imgname)


    imgname_crosshair = "crosshair" + str(n) + "times" + str(m)
    #crosshair_marker01.plot(H,pointdata[0],fermi_energy,Rx,Ry,filepath,imgname_crosshair)
    imgname_local = "local_C_from_ch" + str(n) + "times" + str(m) 
    #crosshair_marker01.local_marker(H,pointdata[0],fermi_energy,filepath,imgname_local)

    imgname_localC = "local_C" + str(n) + "times" + str(m) 
    #local_chern_marker01.plot(H,pointdata[0],fermi_energy,filepath,imgname_localC)

    imgname_yfourier = "ky_band" + str(n) + "times" + str(m) 
    xnumber =1
    y_band_plot(pointdata[0],H,imgname_yfourier)
    imgname_yfourier = "ky_band" + str(n) + "times" + str(m) +"(02)"
    return 0

if __name__ == "__main__":
    main()