import matplotlib.pyplot as plt
import math
import numpy as np
import honeycomb_lattice
import sys
sys.path.append('..')
from crosshair import crosshair_marker01
from local_chern_marker import local_chern_marker01


a = 1#lattice constant
n = 10 #systemsize x axis
m = 10 #systemsize y axis
pointdata = honeycomb_lattice.vertex_create(a,n,m)
#hamiltonian parameter
M = 1.1 #haldane parameter site energy
t1 = 1.0 #haldane nn hopping parameter
t2 = 0.3 + 0j #haldae nnn hopping parameter
phi = math.pi*0.17 #local flux mod pi
filepath = "/Users/sogenikegami/Documents/UT4S/non-crystal/honeycomb/image/"

#fermi_energy = 0.1
Rx = 4 #crosshair position
Ry = 4 #crosshair position

def hsite(pointdata,M):
    vertexdata = pointdata[0]
    plaquettedata = pointdata[1]
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
    plaquettedata = pointdata[1]
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
    H = hsite(pointdata,M) - hhop1(pointdata,t1) - hhop2(pointdata,t2,phi)
    return H


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
    imgname = "haldane_eigval" + str(n) + "times" + str(m)+"(b)"
    fermi_energy = eigenenergy(H,filepath,imgname)


    imgname_crosshair = "haldane_crosshair" + str(n) + "times" + str(m)
    #crosshair_marker01.plot(H,pointdata[0],fermi_energy,Rx,Ry,filepath,imgname_crosshair)
    imgname_local = "local_C_from_ch" + str(n) + "times" + str(m) +"(b)"
    crosshair_marker01.local_marker(H,pointdata[0],fermi_energy,filepath,imgname_local)

    imgname_localC = "local_C" + str(n) + "times" + str(m) +"(b)"
    local_chern_marker01.plot(H,pointdata[0],fermi_energy,filepath,imgname_localC)

if __name__ == "__main__":
    main()