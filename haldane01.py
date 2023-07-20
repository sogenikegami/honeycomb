import matplotlib.pyplot as plt
import math
import numpy as np
import honeycomb_lattice
import sys
sys.path.append('..')
from crosshair import crosshair_marker01
from local_chern_marker import local_chern_marker01


a = 1#lattice constant
n = 9 #systemsize x axis
m = 9 #systemsize y axis

#boundary condition 
xPBC =True
yPBC =False

pointdata = honeycomb_lattice.vertex_create(a,n,m,xPBC,yPBC)

#hamiltonian parameter
M = 1.0 #haldane parameter site energy
t1 = 1.0 #haldane nn hopping parameter
t2 = 1/3 #haldae nnn hopping parameter
phi = math.pi/3 #local flux mod pi
filepath = "/Users/sogenikegami/Documents/UT4S/non-crystal/honeycomb/image/9times9xPBC/"

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
            eiphi = math.cos(phi) + 1j * math.sin(phi) * nnn[j][1] * (1-2*(i%2))
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
                    h[i][j] += phase/a/math.sqrt(3)/m * H[i+2*(n+1)*k][j+l*2*(n+1)]
    eigval,eigvec = np.linalg.eig(h)
    energy = []
    for i in range(len(eigval)):
        energy.append(eigval[i].real)
    energy.sort()
    return energy

def y_band_plot(vertexdata,H,imgname,fermi_energy):
    L = 2*math.pi/math.sqrt(3)/a/(m+1)
    K = np.linspace(-L*(m+1)/2,L*(m+1)/2,m+1)
    E = np.zeros((2*n+2,m+1))
    fig=plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i in range(m+1):
        energy = fourier_y(vertexdata,H,(i-(m+1)/2)*L,)
        for j in range(len(energy)):
            E[j][i] += energy[j]
    for i in range(2*n+2):
        plt.plot(K,E[i],color='black')
    plt.hlines(fermi_energy,-L*(m+1)/2,L*(m+1)/2,linestyle="dashed",color='red')
    plt.savefig(filepath + imgname)
    fig.show()
    return E

def fourier_x(vertexdata,H,kx):
    N = len(vertexdata)

    h = np.zeros((4*(m+1),4*(m+1)))*0j
    for i in range(4*(m+1)):
        for j in range(4*(m+1)):
            for k in range(int((n+1)/2)):
                for l in range(int((n+1)/2)):
                    xk = vertexdata[2*(n+1)*int(i/4)+4*k+(i%4)]["pos"][0]
                    xl = vertexdata[2*(n+1)*int(j/4)+4*l+(j%4)]["pos"][0]
                    phase = math.cos(kx*xl-kx*xk) + 1j*math.sin(kx*xl-kx*xk)
                    h[i][j] += phase/2/math.pi * H[2*(n+1)*int(i/4)+4*k+(i%4)][2*(n+1)*int(j/4)+4*l+(j%4)]
    eigval,eigvec = np.linalg.eigh(h)
    energy = []
    for i in range(len(eigval)):
        energy.append(eigval[i].real)
    energy.sort()
    return energy

def x_band_plot(vertexdata,H,imgname):
    L = 2*math.pi/3/a/(n+1)*2
    K = np.linspace(-L*(n+1)/4,L*(n+1)/4,int((n+1)/2))
    E = np.zeros((4*(m+1),int((n+1)/2)))
    fig=plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i in range(int((n+1)/2)):
        energy = fourier_x(vertexdata,H,(i-(n+1)/4)*L)
        for j in range(len(energy)):
            E[j][i] += energy[j]
    for i in range(4*(m+1)):
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

def plotmap(vertexdata,list1,imgname):
    fig = plt.figure(figsize=(8.66+2,10))
    ax = fig.add_subplot(1,1,1)
    X= []
    Y = []
    for i in range(len(vertexdata)):
        X.append(vertexdata[i]["pos"][0])
        Y.append(vertexdata[i]["pos"][1])
        for j in range(len(vertexdata[i]["neighbor"])):
            neinum = vertexdata[i]["neighbor"][j]
            neipos = vertexdata[neinum]["pos"]
            nowpos = vertexdata[i]["pos"]
            plt.plot([nowpos[0],neipos[0]],[nowpos[1],neipos[1]],color = 'black',linewidth = 1)
    plt.vlines(6,0,math.sqrt(3)*a*m,linestyle="dashed",color='green')
    markermax = max(np.abs(list1))
    mappable = ax.scatter(X,Y,c=list1,cmap='bwr',s=800*np.abs(list1)/markermax,alpha=1,linewidth=0.5,edgecolors='black')
    fig.colorbar(mappable,ax=ax)
    fig.show()
    plt.savefig(filepath + imgname)

def diff(local_c,local_c_from_ch,imgname):
    difflist = []
    for i in range(len(local_c)):
        difflist.append(local_c[i]-local_c_from_ch[i])
    plotmap(pointdata[0],difflist,imgname)
    return difflist

def diff1d_plot(vertexdata,nx,local_c,local_c_from_ch,imgname):
    Y = []
    X = []
    difflist = []
    lc_ch = []
    lc = []
    difference = []
    for i in range(len(local_c)):
        difflist.append(local_c[i]-local_c_from_ch[i])
    for i in range(m):
        Y.append(vertexdata[nx+2*(n+1)*i]["pos"][1])
        X.append(vertexdata[nx+2*(n+1)*i]["pos"][0])
        lc_ch.append(local_c_from_ch[nx+2*(n+1)*i])
        lc.append(local_c[nx+2*(n+1)*i])
        difference.append(difflist[nx+2*(n+1)*i])
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.plot(Y,difference,color = 'black',linewidth = 1,label = 'difference')
    plt.plot(Y,lc,color = 'red',linewidth = 1,label = 'local_C')
    plt.plot(Y,lc_ch,color = 'blue',linewidth = 1,label = 'local_C_from_crosshair')
    plt.hlines(-1,0,Y[-1],linestyle="dashed",color='green')
    plt.legend()
    plt.savefig(filepath + imgname)


def main():
    H = hamiltonian(pointdata,M,t1,t2,phi)
    imgname = "haldane_eigval" + str(n) + "times" + str(m) +"yPBC"
    fermi_energy = eigenenergy(H,filepath,imgname)

    #crosshair
    imgname_crosshair = "crosshair" + str(n) + "times" + str(m) +"yPBC"
    crosshair_marker01.plot(H,pointdata[0],fermi_energy,Rx,Ry,filepath,imgname_crosshair)

    #local chern from crosshair
    imgname_local = "local_C_from_ch" + str(n) + "times" + str(m) +"yPBC"
    local_c_from_ch = crosshair_marker01.local_marker(H,pointdata[0],fermi_energy,filepath,imgname_local)

    #local chern
    imgname_localC = "local_C" + str(n) + "times" + str(m) +"yPBC"
    local_c = local_chern_marker01.plot(H,pointdata[0],fermi_energy,filepath,imgname_localC)

    #yband
    imgname_yfourier = "ky_band" + str(n) + "times" + str(m) +"yPBC"
    #Epbc = y_band_plot(pointdata[0],H,imgname_yfourier,fermi_energy)

    #xband
    imgname4xfourier = "kx_band" + str(n) + "times" + str(m) +"xPBC"
    x_band_plot(pointdata[0],H,imgname4xfourier)

    #difference btw local_c and local_c from crosshair
    imgname4diff = 'diff_localC'+ str(n) + "times" + str(m) +"yPBC"
    difflist = diff(local_c,local_c_from_ch,imgname4diff)   

    #difference 1D cut plot
    imgname4diff_1D = 'diff_1D'+ str(n) + "times" + str(m) +"yPBC"
    #diff1d_plot(pointdata[0],8,local_c,local_c_from_ch,imgname4diff_1D)
    return 0

if __name__ == "__main__":
    main()