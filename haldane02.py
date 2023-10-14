import matplotlib.pyplot as plt
import matplotlib.colors as colors
import math
import cmath
import numpy as np
import pandas as pd
import honeycomb_tetra_lattice
import sys
import os
from scipy.linalg import logm
sys.path.append('..')
from crosshair import crosshair_marker01
from local_chern_marker import local_chern_marker01
from bott import bott_index01
from mpl_toolkits.mplot3d import Axes3D
import honeycomb_tetra_lattice
from scipy.interpolate import griddata

filepath = '/Users/sogenikegami/Documents/UT4S/non-crystal/honeycomb/image/clover_honeycomb/clover/hopping_weaken/alpha01/beta01/PBC/32_32/'

#lattice parameter
a=1
n=32
m=32
xPBC=True
yPBC=False
nnn = True
#nnn=True: next nearest neighbor hopping via defect point exists.
#nnn=False: nnn hopping via defect point does NOT exist.
#defectratio=0.2
mode=0
#mode0 : defect sites are chosen from the defect points in bishamon-kikko lattice
#mode1 : defect sites are randomly chosen.
hopping_weaken = True
nnn_weaken = True
alpha = 0.1
beta = 0.1

#Hamiltonian/energy parameter
#fe_control=1/2*(1-defectratio)+2/5*defectratio
fe_control = 0.666
M = 0.5
t1= 1.0
t2= 1/3
phi = math.pi /3
#crosshair position
Rx=18
Ry=18
#see the value on 1D
x_dir=True
y_dir=False
rx_num=30 #30
ry_num=20 #18

integral_factor = math.sqrt(3)*3/4*a*a

def main():
    #parameter_var()
    defectlist = [0]
    func_defect(defectlist)
    return 0

def parameter_var():
    flag = True
    defectratio = 1
    chern = [[],[],[],[],[]]
    delta_devided_by_t2 =[]
    phi_list = []
    for Delta1 in range(49):
        Delta = (-6 + Delta1*0.25)*t2
        delta_devided_by_t2.append(Delta/t2)
        for phase1 in range(41):
            phase = -1 + phase1 *0.05
            if Delta == -1: phi_list.append(phase)
            M = Delta
            phi = phase*math.pi
            M_str=''
            phi_str=''
            for i in range(len(str(M))):
                if str(M)[i] == '.': M_str+='_'
                else: M_str+=str(M)[i]
            for i in range(len(str(phase))):
                if str(phase)[i] == '.': phi_str+='_'
                else: phi_str+=str(phase)[i]


            vertexdata = honeycomb_tetra_lattice.vertex_create(a,n,m,xPBC,yPBC,mode,defectratio,nnn)
            H = Hamiltonian(vertexdata,M,t1,t2,phi,hopping_weaken,alpha,nnn_weaken,beta)
            imgname1 = "eigval"
            energy = eigenenergy(H,filepath,imgname1)
            eigval = energy[0]
            eigvec = energy[1]
            fermi_energy = energy[2]

            bott_index=0
            #bott_index = bott_index01.B(vertexdata,H,eigval,eigvec,fermi_energy)
            #print('bott index:' + str(bott_index))

            #crosshair
            imgname_crosshair = "crosshair" 
            #crosshairlist = crosshair_marker01.crosshair(H,eigval,eigvec,vertexdata,fermi_energy,Rx,Ry)
            #plotmap(vertexdata,crosshairlist,crosshairlist,'bwr',imgname_crosshair,800,filepath)

            #local chern
            imgname_localC = "phi" + phi_str + "_M_" + M_str
            #local_c = local_chern_marker01.local_chern_marker(H,eigval,eigvec,vertexdata,fermi_energy)
            filepath4localc2d = '/Users/sogenikegami/Documents/UT4S/non-crystal/honeycomb/image/clover_honeycomb/clover/defect1/PBC/phase_diag/local_chern_2d/'
            #plotmap(vertexdata,local_c,local_c,'bwr',imgname_localC,800,filepath4localc2d)

            filepath4localc1d = '/Users/sogenikegami/Documents/UT4S/non-crystal/honeycomb/image/clover_honeycomb/clover/defect1/PBC/phase_diag/local_chern_1d/'
            imgname41d = "phi" + phi_str + "_M_" + M_str
            #one_d_plot(vertexdata,local_c,filepath4localc1d,imgname41d)

            #wavefunc_plot(H,eigval,eigvec,vertexdata,filepath)
            nonzero_list=[]
            anomaly_list=[]
            h_axis=[]
            bott_indexes=[]
            chern_num=[]
            #nonzero_list,anomaly_list,h_axis,bott_indexes = bott_energy(vertexdata,H,eigval,eigvec,5,filepath,defectratio)
            #nonzero_list,anomaly_list = bott_add(vertexdata,H,eigval,eigvec,10)
            #print(bott_extra(vertexdata,H,eigval,eigvec,10))

            #bishamon-band
            #x_band_plot(vertexdata,H,"kx_band" )
            filepath4band = '/Users/sogenikegami/Documents/UT4S/non-crystal/honeycomb/image/clover_honeycomb/clover/defect1/PBC/phase_diag/band/'
            #xy_band_plot(a,"phi" + phi_str + "_M_" + M_str,filepath4band,M,phi,band_extension)
            #band_plot_3d(vertexdata,H,M,phi)

            #honeycomb-band
            #haldane_band("honeycomb_band",20)

            #Chern-number
            filepath4berrycurv = '/Users/sogenikegami/Documents/UT4S/non-crystal/honeycomb/image/clover_honeycomb/clover/defect1/PBC/phase_diag/'
            chern_num = chern_number4(50,filepath4berrycurv,M,phi)
            for i in range(5):
                chern[i].append(chern_num[i])
            #print(haldane_chern2(10))
            #print(chern_num)
            if flag:
                filepath4winfo = '/Users/sogenikegami/Documents/UT4S/non-crystal/honeycomb/image/clover_honeycomb/clover/nnnoff/defect10/PBC/phase_diag/'
                writeinfo(bott_index,H,nonzero_list,anomaly_list,h_axis,bott_indexes,chern_num,filepath4winfo,defectratio)
                flag = False
    with open('/Users/sogenikegami/Documents/UT4S/non-crystal/honeycomb/image/clover_honeycomb/clover/nnnoff/defect10/PBC/phase_diag/output3.txt' , 'w') as f:
        f.write('Delta_axis:' + str(delta_devided_by_t2) + '\n')
        f.write('phi_devided_by_t2: ' + str(phi_list) + '\n')
        for i in range(5):
            f.write('#' +str(i)+':' + str(chern[i])+'\n')
    
    return 0

def func_defect(defectlist):
    
    for i in range(len(defectlist)):
        defectratio = defectlist[i]/10
        #filepath = '/Users/sogenikegami/Documents/UT4S/non-crystal/honeycomb/image/clover_honeycomb/clover/defect' + str(defectlist[i]) + '/xPBC/m-12/32_32_3/'
        #os.makedirs(filepath)

        vertexdata = honeycomb_tetra_lattice.vertex_create(a,n,m,xPBC,yPBC,mode,defectratio,nnn)
        H = Hamiltonian(vertexdata,M,t1,t2,phi,hopping_weaken,alpha,nnn_weaken,beta)
        
        np.set_printoptions(threshold=5000)
        #print(H)
        imgname1 = "eigval"
        energy = eigenenergy(H,filepath,imgname1)
        eigval = energy[0]
        eigvec = energy[1]
        fermi_energy = energy[2]

        bott_index = bott_index01.B(vertexdata,H,eigval,eigvec,fermi_energy)
        print('bott index:' + str(bott_index))

        #crosshair
        imgname_crosshair = "crosshair" 
        #crosshairlist = crosshair_marker01.crosshair(H,eigval,eigvec,vertexdata,fermi_energy,Rx,Ry)
        #plotmap(vertexdata,crosshairlist,crosshairlist,'bwr',imgname_crosshair,800,filepath)

        #local chern
        imgname_localC = "local_C"
        local_c = local_chern_marker01.local_chern_marker(H,eigval,eigvec,vertexdata,fermi_energy)
        plotmap(vertexdata,local_c,local_c,'bwr',imgname_localC,800,filepath)

        one_d_plot(vertexdata,local_c,filepath,'local_c_xdir')

        wavefunc_plot(H,eigval,eigvec,vertexdata,filepath)
        nonzero_list=[]
        anomaly_list=[]
        h_axis=[]
        bott_indexes=[]
        chern_num=[]
        #nonzero_list,anomaly_list,h_axis,bott_indexes = bott_energy(vertexdata,H,eigval,eigvec,5,filepath,defectratio)
        #nonzero_list,anomaly_list = bott_add(vertexdata,H,eigval,eigvec,10)
        #print(bott_extra(vertexdata,H,eigval,eigvec,10))

        #bishamon-band
        #x_band_plot(vertexdata,H,"kx_band" )
        #xy_band_plot(a,'2D-band',filepath,M,phi,True)
        #band_plot_3d(vertexdata,H)
        #chern_num = chern_number4(50,filepath,M,phi)

        #honeycomb-band
        #haldane_band("honeycomb_band",20)
        #chern_num = haldane_chern2(50)

        #hopping_weaken
        hopping_weaken_band(True,'band')
        chern_num = hopping_weaken_chern(70,filepath,M,phi,alpha,t1,t2,beta)

        print(chern_num)
        writeinfo(bott_index,H,nonzero_list,anomaly_list,h_axis,bott_indexes,chern_num,filepath,defectratio)
    return 0

def Hamiltonian(vertexdata,M,t1,t2,phi,hopping_weaken,alpha,nnn_weaken,beta):
    H = hsite(vertexdata,M) + hhop1(vertexdata,t1,hopping_weaken,alpha) + hhop2(vertexdata,t2,phi,hopping_weaken,alpha,nnn_weaken,beta)
    return H

def hsite(vertexdata,M):
    N = len(vertexdata)
    H = np.zeros((N,N))*0j
    for i in range(N):
        H[i][i] += M*(vertexdata[i]["parity"]*2-1)
    return H

def hhop1(vertexdata,t1,hopping_weaken,alpha):
    N = len(vertexdata)
    H = np.zeros((N,N))*0j
    candidate = candidate_gen()
    if hopping_weaken:
        for i in range(N):
            neighbor = vertexdata[i]["neighbor"]
            for j in range(len(neighbor)):
                if (i in candidate) or (neighbor[j] in candidate):
                    H[neighbor[j]][i] += alpha*t1
                else:
                    H[neighbor[j]][i] += t1
    else:
        for i in range(N):
            neighbor = vertexdata[i]["neighbor"]
            for j in range(len(neighbor)):
                H[neighbor[j]][i] += t1
    return H

def hhop2(vertexdata,t2,phi,hopping_weaken,alpha,nnn_weaken,beta):
    N = len(vertexdata)
    H = np.zeros((N,N))*0j
    if hopping_weaken:
        candidate = candidate_gen()
        if nnn_weaken:
            for i in range(N):
                neighbor = vertexdata[i]["neighbor"]
                for j in range(len(neighbor)):
                    nnn = vertexdata[neighbor[j]]["neighbor"]
                    for k in range(len(nnn)):
                        nnnnum = nnn[k]
                        if i!=nnnnum:
                            for l in range(len(vertexdata[i]["secondneighbor"])):
                                if nnnnum == vertexdata[i]["secondneighbor"][l][0]:
                                    if (i in candidate) or (nnnnum in candidate):
                                        H[nnnnum][i] += alpha * t2*(math.cos(phi)+math.sin(phi)*1j*vertexdata[i]["secondneighbor"][l][1])
                                    elif neighbor[j] in candidate:
                                        H[nnnnum][i] += beta * t2*(math.cos(phi)+math.sin(phi)*1j*vertexdata[i]["secondneighbor"][l][1])
                                    else:
                                        H[nnnnum][i] += t2*(math.cos(phi)+math.sin(phi)*1j*vertexdata[i]["secondneighbor"][l][1])                        
        else:
            for i in range(N):
                nnn = vertexdata[i]["secondneighbor"]
                for j in range(len(nnn)):
                    if (i in candidate) or (nnn[j][0] in candidate):
                        H[nnn[j][0]][i] += alpha * t2*(math.cos(phi)+math.sin(phi)*1j*nnn[j][1])
                    else:
                        H[nnn[j][0]][i] += t2*(math.cos(phi)+math.sin(phi)*1j*nnn[j][1])
    else:
        if nnn_weaken:
            candidate = candidate_gen()
            for i in range(N):
                neighbor = vertexdata[i]["neighbor"]
                for j in range(len(neighbor)):
                    nnn = vertexdata[neighbor[j]]["neighbor"]
                    for k in range(len(nnn)):
                        nnnnum = nnn[k]
                        if i!=nnnnum:
                            for l in range(len(vertexdata[i]["secondneighbor"])):
                                if nnnnum == vertexdata[i]["secondneighbor"][l][0]:
                                    if neighbor[j] in candidate:
                                        H[nnnnum[0]][i] += beta * t2*(math.cos(phi)+math.sin(phi)*1j*vertexdata[i]["secondneighbor"][l][1])
                                    else:
                                        H[nnnnum[0]][i] += t2*(math.cos(phi)+math.sin(phi)*1j*vertexdata[i]["secondneighbor"][l][1])      
        else:
            for i in range(N):
                nnn = vertexdata[i]["secondneighbor"]
                for j in range(len(nnn)):
                    H[nnn[j][0]][i] += t2*(math.cos(phi)+math.sin(phi)*1j*nnn[j][1])
    return H

def candidate_gen():
    candidate = []
    if xPBC:
        for i in range(int(m/4)):
            for j in range(int(n/2)+1):
                candidate.append(2*n+7+(6*n+12)*i+4*j)
            for j in range(int(n/2)+1):
                candidate.append(4*n+11+(6*n+12)*i+4*j)
        if yPBC:
            for i in range(int(n/2)+1):
                candidate.append(int(2*n+5+(6*n+12)*m/4+2*i))
    else:
        for i in range(int(m/4)):
            for j in range(int(n/2)):
                candidate.append(2*n+4+(6*n+8)*i+4*j)
            for j in range(int(n/2)):
                candidate.append(4*n+7+(6*n+8)*i+4*j)
        if yPBC:
            for i in range(int(n/2)):
                candidate.append(2*n+4+(6*n+8)*int(m/4)+2*i)

    return candidate


def bott_energy_one(vertexdata,H,eigval,eigvec,h_axis):
    eigenergy = []
    bott_index = []        
    for i in range(len(eigval)):
        eigenergy.append(eigval[i].real)
    eigenergy.sort()
    for i in range(len(h_axis)):
        E = eigenergy[int(len(eigenergy)*h_axis[i])]
        bott_index.append(bott_index01.B(vertexdata,H,eigval,eigvec,E))

    return bott_index

def bott_energy(vertexdata,H,eigval,eigvec,iter,filepath,defectratio):
    """
    h_axis=[0]
    for i in range(5):
        for j in range(3): 
            if i==2 and j==1:
                h_axis.append(0.46)
                h_axis.append(0.47)
                h_axis.append(0.48)
                h_axis.append(0.49)
            elif i==2 and j==2:
                h_axis.append(0.51) 
                h_axis.append(0.52)
                h_axis.append(0.53)
                h_axis.append(0.54)
            h_axis.append(0.2*i+0.05+0.05*j)
            if i==2 and j==2:
                h_axis.append(0.56)
                h_axis.append(0.57)
        if i!=4:
            for j in range(11):
                h_axis.append(0.18+0.2*i+0.005*j)
    """
    #h_axis = [0.1,0.2,0.3,0.38,0.4,0.42,0.44,0.46,0.48,0.5,0.52,0.54,0.56,0.58,0.6,0.62,0.7,0.78,0.8,0.82,0.9]
    #h_axis = [0.1,0.2,0.3,0.38,0.39,0.4,0.405,0.41,0.415,0.42,0.425,0.43,0.435,0.44,0.445,0.45,0.455,0.46,0.465,0.47,0.475,0.48,0.485,0.49,0.495,0.5,0.505,0.51,0.515,0.52,0.525,0.53,0.535,0.54,0.545,0.55,0.555,0.56,0.565,0.57,0.575,0.58,0.585,0.59,0.595,0.6,0.61,0.62,0.63,0.64,0.65,0.7,0.78,0.79,0.8,0.81,0.82,0.9]
    h_axis = [0.35,0.36,0.37]
    bott_index=np.zeros(len(h_axis))
    for i in range(iter):
        if i==0:
            bott_one = bott_energy_one(vertexdata,H,eigval,eigvec,h_axis)
            for j in range(len(bott_one)):
                bott_index[j] += bott_one[j]/iter
        else:
            vertexdata = honeycomb_tetra_lattice.vertex_create(a,n,m,xPBC,yPBC,mode,defectratio,nnn)
            H = Hamiltonian(vertexdata,M,t1,t2,hopping_weaken,alpha)
            eigval,eigvec = np.linalg.eigh(H)
            bott_one = bott_energy_one(vertexdata,H,eigval,eigvec,h_axis)
            for j in range(len(bott_one)):
                bott_index[j] += bott_one[j]/iter

    nonzero_list=[]
    anomaly_list = []
    for i in range(len(bott_index)):
        if abs(bott_index[i]) > 0.05:
            nonzero_list.append(h_axis[i])
        if abs(bott_index[i]) > 1.05:
            anomaly_list.append(h_axis[i])

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.plot(h_axis,bott_index,marker='.',markersize=10)
    ax.set_ylim(-1.1,1.1)
    plt.savefig(filepath + 'bott_index_fe')
    return [nonzero_list,anomaly_list,h_axis,bott_index]
    
def bott_extra(vertexdata,H,eigval,eigvec,iter):
    bott = [0,0,0,0]
    for i in range(iter):
        if i==0:
            eigenergy = []     
            for j in range(len(eigval)):
                eigenergy.append(eigval[j].real)
            eigenergy.sort()
            bott[0] += bott_index01.B(vertexdata,H,eigval,eigvec,eigenergy[int(len(eigenergy)*0.53)])/iter
            bott[1] += bott_index01.B(vertexdata,H,eigval,eigvec,eigenergy[int(len(eigenergy)*0.54)])/iter
            bott[2] += bott_index01.B(vertexdata,H,eigval,eigvec,eigenergy[int(len(eigenergy)*0.56)])/iter
            bott[3] += bott_index01.B(vertexdata,H,eigval,eigvec,eigenergy[int(len(eigenergy)*0.57)])/iter

        else:
            vertexdata = honeycomb_tetra_lattice.vertex_create(a,n,m,xPBC,yPBC,mode,defectratio,nnn)
            H = Hamiltonian(vertexdata,M,t1,t2,hopping_weaken,alpha)
            eigval,eigvec = np.linalg.eigh(H)
            eigenergy = []     
            for j in range(len(eigval)):
                eigenergy.append(eigval[j].real)
            eigenergy.sort()
            bott[0] += bott_index01.B(vertexdata,H,eigval,eigvec,eigenergy[int(len(eigenergy)*0.53)])/iter
            bott[1] += bott_index01.B(vertexdata,H,eigval,eigvec,eigenergy[int(len(eigenergy)*0.54)])/iter
            bott[2] += bott_index01.B(vertexdata,H,eigval,eigvec,eigenergy[int(len(eigenergy)*0.56)])/iter
            bott[3] += bott_index01.B(vertexdata,H,eigval,eigvec,eigenergy[int(len(eigenergy)*0.57)])/iter

    return bott


def bott_add(vertexdata,H,eigval,eigvec,iter):
    h_axis=[0]
    for i in range(5):
        for j in range(3): 
            if i==2 and j==1:
                h_axis.append(0.46)
                h_axis.append(0.47)
                h_axis.append(0.48)
                h_axis.append(0.49)
            elif i==2 and j==2:
                h_axis.append(0.51) 
                h_axis.append(0.52)
                h_axis.append(0.53)
                h_axis.append(0.54)
            h_axis.append(0.2*i+0.05+0.05*j)
            if i==2 and j==2:
                h_axis.append(0.56)
                h_axis.append(0.57)
        if i!=4:
            for j in range(11):
                h_axis.append(0.18+0.2*i+0.005*j)
    bott = []
    for i in range(66):
        if i==32: bott.append(0)
        elif i==33: bott.append(0)
        else: bott.append(0)
    
    for i in range(iter):
        if i==0:
            eigenergy = []     
            for j in range(len(eigval)):
                eigenergy.append(eigval[j].real)
            eigenergy.sort()
            E = eigenergy[int(len(eigenergy)*0.46)]
            bott[30] += bott_index01.B(vertexdata,H,eigval,eigvec,E)/iter

        else:
            vertexdata = honeycomb_tetra_lattice.vertex_create(a,n,m,xPBC,yPBC,mode,defectratio,nnn)
            H = Hamiltonian(vertexdata,M,t1,t2,hopping_weaken,alpha)
            eigval,eigvec = np.linalg.eigh(H)
            eigenergy = []     
            for j in range(len(eigval)):
                eigenergy.append(eigval[j].real)
            eigenergy.sort()
            E = eigenergy[int(len(eigenergy)*0.46)]
            bott[30] += bott_index01.B(vertexdata,H,eigval,eigvec,E)/iter

    for i in range(iter):
        if i==0:
            eigenergy = []     
            for j in range(len(eigval)):
                eigenergy.append(eigval[j].real)
            eigenergy.sort()
            E = eigenergy[int(len(eigenergy)*0.47)]
            bott[31] += bott_index01.B(vertexdata,H,eigval,eigvec,E)/iter

        else:
            vertexdata = honeycomb_tetra_lattice.vertex_create(a,n,m,xPBC,yPBC,mode,defectratio,nnn)
            H = Hamiltonian(vertexdata,M,t1,t2,hopping_weaken,alpha)
            eigval,eigvec = np.linalg.eigh(H)
            eigenergy = []     
            for j in range(len(eigval)):
                eigenergy.append(eigval[j].real)
            eigenergy.sort()
            E = eigenergy[int(len(eigenergy)*0.47)]
            bott[31] += bott_index01.B(vertexdata,H,eigval,eigvec,E)/iter
    
    nonzero=[]
    anomaly = []
    for i in range(len(bott)):
        if bott[i] >= 0.1: nonzero.append(h_axis[i])
        if bott[i] > 1.05: anomaly.append(h_axis[i])
    
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.plot(h_axis,bott,marker='.',markersize=10)
    ax.set_ylim(-1.1,1.1)
    plt.savefig(filepath + 'bott_index_fe_2')
    return [nonzero,anomaly]
            


def eigenenergy(H,filepath,imgname):
    eigval,eigvec = np.linalg.eigh(H)
    num = list(range(0,len(eigval)))
    eigenergy = []
    for i in range(len(eigval)):
        eigenergy.append(eigval[i].real)
    eigenergy.sort()
    #fermi_energy = (eigenergy[int(len(H[0])*fe_control)]+eigenergy[int(len(H[0])*fe_control)-1])/2
    fermi_energy = eigenergy[int(len(H[0])*fe_control)]
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.scatter(num,eigenergy,s=10)
    plt.hlines(fermi_energy,0,len(eigval),linestyle="dashed",color='red')
    fig.show()
    plt.savefig(filepath + imgname)
    
    return [eigval,eigvec,fermi_energy]

def wavefunc_plot(H,eigval,eigvec,vertexdata,filepath):
    absval=[]
    phase=[]
    for j in range(len(eigval)):
        absval.append(abs(eigvec[j][int(len(H[0])*fe_control)]))
        phase.append(cmath.phase(eigvec[j][int(len(H[0])*fe_control)]))
    plotmap(vertexdata,phase,absval,'jet','wavefunc',300,filepath)
    return 0

def plotmap(vertexdata,listc,lists,mapcolor,imgname,Size,filepath):
    fig = plt.figure(figsize=(24,6))
    ax = fig.add_subplot(1,1,1)
    if xPBC: shiftx = 3*a*(n/2+1)
    else: shiftx = 0
    if yPBC: shifty = [3/2*a*(m/2+1),3/2*math.sqrt(3)*a*(m/2+1)]
    else: shifty=[0,0]
    X=[]
    Y=[]

    for i in range(len(vertexdata)):
        X.append(vertexdata[i]["pos"][0])
        Y.append(vertexdata[i]["pos"][1])
        for j in range(len(vertexdata[i]["neighbor"])):
            dot = False
            nowpos = vertexdata[i]["pos"]
            neinum = vertexdata[i]["neighbor"][j]
            neipos0= vertexdata[neinum]["pos"][0]
            neipos1= vertexdata[neinum]["pos"][1]
            
            if nowpos[1]-neipos1 > 2*a:
                neipos0 += shifty[0]
                neipos1 += shifty[1]
                dot = True
            elif neipos1-nowpos[1] > 2*a:
                neipos0 -= shifty[0]
                neipos1 -= shifty[1]
                dot = True
            
            if nowpos[0]-neipos0 > 3/2*a*(m/2+1):  
                neipos0+=shiftx
                dot = True
            elif neipos0-nowpos[0] > 3/2*a*(m/2+1): 
                neipos0-=shiftx
                dot = True
            if dot: plt.plot([nowpos[0],neipos0],[nowpos[1],neipos1],color="blue",linestyle="dotted")
            else:  plt.plot([nowpos[0],neipos0],[nowpos[1],neipos1],color="black")

    if "crosshair" in str(imgname):
        ymax = max(vertexdata[-1]["pos"][1],vertexdata[-2]["pos"][1])
        xmax = vertexdata[-1]["pos"][0]
        plt.vlines(Rx,0,ymax,linestyle="dashed",color='pink')
        plt.hlines(Ry,0,xmax,linestyle="dashed",color='pink')
    if y_dir:
        rx=vertexdata[rx_num]["pos"][0]
        ymax = max(vertexdata[-1]["pos"][1],vertexdata[-2]["pos"][1])
        plt.vlines(rx,0,ymax,linestyle="dashed",color='purple')
    if x_dir:
        xmax = vertexdata[-1]["pos"][0]
        ry=ry_num*a*math.sqrt(3)/2
        plt.hlines(ry,0,xmax,linestyle="dashed",color='green')
    markermax = max(np.abs(lists))
    #mappable = ax.scatter(X,Y,c=listc,cmap=mapcolor,s=Size*np.abs(lists)/markermax,alpha=1,linewidth=0.5,edgecolors='black',norm=colors.TwoSlopeNorm(0))
    mappable = ax.scatter(X,Y,c=listc,cmap=mapcolor,s=100,alpha=1,linewidth=0.5,edgecolors='black',norm=colors.TwoSlopeNorm(0))
    fig.colorbar(mappable,ax=ax)
    fig.show()
    plt.savefig(filepath + imgname)

    return 0

def one_d_plot(vertexdata,list,filepath,imgname):
    X = []
    Y = []
    listx=[]
    listy=[]
    rx=vertexdata[rx_num]["pos"][0]
    ry=ry_num*a*math.sqrt(3)/2
    if x_dir:
        for i in range(len(vertexdata)):
            if abs(vertexdata[i]["pos"][1]-ry)<0.0001:
                X.append(vertexdata[i]["pos"][0])
                listx.append(list[i])
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        xmin = X[0]
        xmax = X[-1]
        plt.hlines(1,xmin,xmax,linestyle="dashed",color='green')
        plt.rcParams["font.size"] = 15
        plt.plot(X,listx,color = 'red',linewidth = 1,label = 'local chern marker')
        plt.legend()
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.xlabel("x",fontsize=15)
        plt.savefig(filepath + imgname)
    if y_dir:
        for i in range(len(vertexdata)):
            if abs(vertexdata[i]["pos"][0]-rx)<0.0001:
                Y.append(vertexdata[i]["pos"][1])
                listy.append(list[i])
        fig = plt.figure()
        ax = fig.add_subplot(1,1,1)
        plt.plot(Y,listy,color = 'red',linewidth = 1,label = 'local Chern marker')
        plt.legend()
        plt.xlabel("y")
        plt.savefig(filepath +  'local_chern_ydir')
    return 0

def writeinfo(bott,H,nonzero_list,anomaly_list,h_axis,bott_indexes,chern_num,filepath,defectratio):
    filename = 'parameter_information.txt'
    is_file = os.path.isfile(filepath + filename)
    with open(filepath + filename,'w') as f:
        f.write("bott index B=" + str(bott) +'\n')
        f.write("lattice constant a=" +str(a)+'\n')
        if xPBC:    f.write("boundary condition xdir:periodic"+'\n')
        else:  f.write("boundary condition xdir:open"+'\n')
        if yPBC:    f.write("boundary condition ydir:periodic"+'\n')
        else:  f.write("boundary condition ydir:open"+'\n')
        f.write("system size n=" + str(n)+'\n')
        f.write("system size m=" + str(m)+'\n')
        if nnn: f.write("nnn hopping via defect: ON" +'\n')
        else: f.write("nnn hopping via defect: OFF" +'\n')
        if hopping_weaken: f.write("defect related hopping weaken: ON" +'\n')
        else: f.write("defect related hopping weaken: OFF" +'\n')
        f.write("defect candidate hopping strength coeff: alpha = " + str(alpha) +'\n')
        f.write("nnn hopping strength coedd via defect cand. beta = " + str(beta) +'\n')
        f.write("Crosshair position Rx=" + str(Rx)+'\n')
        f.write("Crosshair position Ry=" + str(Ry)+'\n')
        f.write("integral factor(weight)=" + str(integral_factor)+'\n')
        #f.write("Hamiltonian mode=" + str(mode)+'\n')
        #f.write("Circle radius r=" + str(r)+'\n')
        #f.write("The number of stripe stripe_num=" + str(stripe_num)+'\n')
        f.write("fermi energy selecting: fe=" + str(fe_control)+'\n')
        #f.write("to illustrate the 1D crosshair value for x-direction, y=" + str(ycutline)+'\n')
        #f.write("lattice number to 1D value, y-direction:"+str(ycutline_num)+'\n')
        #f.write("to illustrate the 1D crosshair value for y-direction, x=" + str(xcutline)+'\n')
        #f.write("lattice number to 1D value, x-direction:"+str(xcutline_num)+'\n')
        f.write("defect ratio:" + str(defectratio) + '\n')
        f.write("site potential: M=" +  str(M) + "\n")
        f.write("t1:" + str(t1) + "\n")
        f.write("t2:" + str(t2) + "\n")
        f.write("phi:" + str(phi) + "\n")
        f.write("system size:" + str(len(H[0])) + '\n')
        f.write("Chern number (from the bottom band):" + str(chern_num) + '\n')
        f.write("nonzero-list:" + str(nonzero_list) + '\n')
        f.write("anomaly-list:" + str(anomaly_list) + '\n')
        f.write("h_axis:" + str(h_axis) + '\n')
        f.write("bott_indexes:" + str(bott_indexes) + '\n')

        f.write(str(H))
    return 0

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
    eigval,eigvec = np.linalg.eigh(h)
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
    h = np.zeros((int(5*m/2+4),int(5*m/2+4)))*0j
    for i in range(int(5*m/2+4)):#y1
        for j in range(int(5*m/2+4)):#y2
            for k in range(int(n/2+1)):
                for l in range(int(n/2+1)):
                    if 0<=i%10 and i%10<4: knum = int(i%10)+10*int(n/2+1)*int(i/10)+4*k
                    elif 4<= i%10 and i%10 < 7: knum = int(i%10)-4+2*n+4+10*int(n/2+1)*int(i/10)+3*k
                    else: knum = int(i%10)-7+7*int(n/2+1)+10*int(n/2+1)*int(i/10)+3*k
                    if 0<=j%10 and j%10<4: lnum = int(j%10)+10*int(n/2+1)*int(j/10)+4*l
                    elif 4<= j%10 and j%10 < 7: lnum = int(j%10)-4+2*n+4+10*int(n/2+1)*int(j/10)+3*l
                    else: lnum = int(j%10)-7+7*int(n/2+1)+10*int(n/2+1)*int(j/10)+3*l
                    xk = vertexdata[knum]["pos"][0]
                    xl = vertexdata[lnum]["pos"][0]
                    phase = math.cos(kx*(xk-xl)) + 1j*math.sin(kx*(xk-xl))
                    h[i][j] += phase/2/math.pi * H[knum][lnum]
    eigval,eigvec = np.linalg.eig(h)
    energy = []
    for i in range(len(eigval)):
        energy.append(eigval[i].real)
    energy.sort()
    return energy


def x_band_plot(vertexdata,H,imgname):
    L = 2*math.pi/3/a/(n/2+1)
    K = np.linspace(-L*n/4,L*n/4,int(n/2+1))
    E = np.zeros((int(5*m/2+4),int(n/2+1)))
    fig=plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i in range(int(n/2+1)):
        kx = K[i]
        energy = fourier_x(vertexdata,H,kx)
        for j in range(len(energy)):
            E[j][i] += energy[j]
    for i in range(int(5*m/2+4)):
        plt.plot(K,E[i],color='black')
    plt.savefig(filepath + imgname)
    fig.show()
    return 0

def fourier_2d(vertexdata,H,kx,ky):
    h = np.zeros((5,5))*0j
    for i in range(5):
        for j in range(5):
            for k in range(int(n/2+1)):
                for l in range(int(n/2+1)):
                    for e in range(int(m/2+1)):
                        for f in range(int(m/2+1)):
                            if i==0:
                                if e==0: num1 = int(4*(n/2+1) + 10+(n/2+1)*m/4 + k)
                                else: num1 = int(2*n+4 + 3*(n/2+1)*(e%2) + 10*(n/2+1)*int((e-1)/2) + 3*k)
                            elif i==1:
                                if e==0 and k==0: num1 = int(2*n+2)
                                elif e==0: num1 = int(2 + (k-1)*4)
                                elif e!=0 and e%2==1: num1 = int(2*n+5 + 10*(n/2+1)*int((e-1)/2) + 3*k)
                                else: num1 = int(10*(n/2+1) + 10*(n/2+1)*int((e-1)/2) + 4*k)
                            elif i==2:
                                if e==0 and k==0: num1 = int(2*n+1)
                                elif e==0: num1 = int(3+(k-1)*4)
                                elif e!=0 and e%2==1: num1 = int(2*n+6 + 10*(n/2+1)*int((e-1)/2) + 3*k)
                                else: num1 = int(10*(n/2+1)+1 + 10*(n/2+1)*int((e-1)/2) + 4*k)
                            elif i==3: 
                                if e==0: num1 = int(4*k)
                                elif e%2==1: num1 = int(7*(n/2+1)+1 + 10*(n/2+1)*int((e-1)/2) + 3*k)
                                else: num1 = int(10*(n/2+1)+2 + 10*(n/2+1)*int((e-1)/2) + 4*k)
                            else:
                                if e==0: num1 = int(1+4*k)
                                elif e%2==1: num1 = int(7*(n/2+1)+2 + 10*(n/2+1)*int((e-1)/2) + 3*k)
                                else: num1 = int(10*(n/2+1)+3 + 10*(n/2+1)*int((e-1)/2) + 4*k)
                            if j==0:
                                if f==0: num2 = int(4*(n/2+1) + 10+(n/2+1)*m/4 + l)
                                else: num2 = int(2*n+4 + 3*(n/2+1)*(f%2) + 10*(n/2+1)*int((f-1)/2) + 3*l)
                            elif j==1:
                                if f==0 and l==0: num2 = int(2*n+2)
                                elif f==0: num2 = int(2 + (l-1)*4)
                                elif f!=0 and f%2==1: num2 = int(2*n+5 + 10*(n/2+1)*int((f-1)/2) + 3*l)
                                else: num2 = int(10*(n/2+1) + 10*(n/2+1)*int((f-1)/2) + 4*l)
                            elif j==2:
                                if f==0 and l==0: num2 = int(2*n+1)
                                elif f==0: num2 = int(3+(l-1)*4)
                                elif f!=0 and f%2==1: num2 = int(2*n+6 + 10*(n/2+1)*int((f-1)/2) + 3*l)
                                else: num2 = int(10*(n/2+1)+1 + 10*(n/2+1)*int((f-1)/2) + 4*l)
                            elif j==3: 
                                if f==0: num2 = int(4*l)
                                elif f%2==1: num2 = int(7*(n/2+1)+1 + 10*(n/2+1)*int((f-1)/2) + 3*l)
                                else: num2 = int(10*(n/2+1)+2 + 10*(n/2+1)*int((f-1)/2) + 4*l)
                            else:
                                if f==0: num2 = int(1+4*l)
                                elif f%2==1: num2 = int(7*(n/2+1)+2 + 10*(n/2+1)*int((f-1)/2) + 3*l)
                                else: num2 = int(10*(n/2+1)+3 + 10*(n/2+1)*int((f-1)/2) + 4*l)
                            pos1 = vertexdata[num1]["pos"]
                            pos2 = vertexdata[num2]["pos"]
                            phasex = (math.cos(kx*(pos1[0]-pos2[0])) - 1j * math.sin(kx*(pos1[0]-pos2[0])))
                            phasey = (math.cos(ky*(pos1[1]-pos2[1])) - 1j * math.sin(ky*(pos1[1]-pos2[1])))
                            phase = phasex*phasey
                            h[i][j] += phase/(9*math.sqrt(3)/2*a*a)**2 * H[num1][num2]

    return h

"""
def fourier_2d_2(vertexdata,H,kx,ky):
    h = np.zeros((5,5))*0j
    for i in range(5):
        for j in range(5):
            if i==0: num1 = 2*n+7
            elif i==1: num1 = 2*n+8
            elif i==2: num1 = 2*n+9
            elif i==3: num1 = int(7*(n/2+1)+4)
            else: num1 = int(7*(n/2+1)+5)
            if j==0: num2 = 2*n+7
            elif j==1: num2 = 2*n+8
            elif j==2: num2 = 2*n+9
            elif j==3: num2 = int(7*(n/2+1)+4)
            else: num2 = int(7*(n/2+1)+5)
            pos1 = vertexdata[num1]["pos"]
            pos2 = vertexdata[num2]["pos"]
            for k in range(len(vertexdata)):
                nowpos = vertexdata[k]["pos"]
                h[i][j]+=H[nowpos][nowpos]
                for l in range(len(vertexdata[nowpos]["neighbor"])):
                    neinum = vertexdata[nowpos]["neighbor"][l]
                    phasex = (math.cos(kx*(pos1[0]-pos2[0])) + 1j * math.sin(kx*(pos1[0]-pos2[0])))
                    phasey = (math.cos(ky*(pos1[1]-pos2[1])) + 1j * math.sin(ky*(pos1[1]-pos2[1])))
                    phase = phasex*phasey
                    h[i][j] += phase*H[k][neinum]
"""

def xy_band_plot(a,imgname,filepath,M,phi,band_extension):
    #gamma-K
    klist = []
    x_axis=[]
    N=20
    """
    for i in range(int(2/3*(n/2+1))):
        klist.append([2*math.pi*i/3/a/(n/2+1),0])
        x_axis.append(2*math.pi*i/3/a/(n/2+1))
    if (int(n/2+1))%3!=0:
        klist.append([4*math.pi/9/a,0])
        x_axis.append(4*math.pi/9/a)

    for i in range(int(2/3*(n/2+1)),int((n/2+1)/2),-1):
        x = 2*math.pi/3/a/(n/2+1)*i
        klist.append([x,-math.sqrt(3)*x + 4*math.pi/3/a/math.sqrt(3)])
        x_axis.append(4*math.pi/9/a + 2*(4*math.pi/9/a-x))
    if (n+2)%4!=0:
        klist.append([math.pi/3/a,math.pi/3/a/math.sqrt(3)])
        x_axis.append(4*math.pi/9/a + 2*(4*math.pi/9/a-math.pi/3/a))
    for i in range(int((n+2)/4),-1,-1):
        x = 2*math.pi/3/a/(n/2+1)*i
        klist.append([x,x/math.sqrt(3)])
        x_axis.append(4*math.pi/9/a + 2*(4*math.pi/9/a-math.pi/3/a) + 2/math.sqrt(3)*(math.pi/3/a-x))
    """
    for i in range(int(2*N)):
        x = 4*math.pi/9/a/2/N*i
        klist.append([x,0])
        x_axis.append(x)
    for i in range(N):
        x = 4*math.pi/9/a - (math.pi/9/a/N*i)
        klist.append([x,(4*math.pi/9/a-x)*math.sqrt(3)])
        x_axis.append(4*math.pi/9/a + (4*math.pi/9/a-x)*2)
    N3 = int(math.sqrt(3)*N)
    for i in range(N3):
        x = math.pi/3/a*(1-i/N3)
        klist.append([x,x/math.sqrt(3)])
        x_axis.append(2*math.pi/3/a+ math.pi/3/a/N3*i*2/math.sqrt(3))
    if band_extension:
        for i in range(N3):
            x = math.pi/3/a*i/N3
            klist.append([x,-x/math.sqrt(3)])
            x_axis.append(2*math.pi/a/9*(3+math.sqrt(3)) + x *2/math.sqrt(3))
        for i in range(N+1):
            x = math.pi/3/a + math.pi/9/a*i/N
            klist.append([x,-(4*math.pi/9/a - x)*math.sqrt(3)])
            x_axis.append(2*math.pi/a/9*(3+2*math.sqrt(3)) + math.pi/9/a*i/N*2 )
    else:
        klist.append([0,0])
        x_axis.append(2*math.pi/9/a*(3+math.sqrt(3)))


    E = np.zeros((5,len(klist)))*0j
    for i in range(len(klist)):
        #h = fourier_2d(vertexdata,H,klist[i][0],klist[i][1])
        if nnn:h = fourier_2d_2(klist[i][0],klist[i][1],M,phi)
        else: h = fourier_2d_nnnoff(klist[i][0],klist[i][1],M,phi)
        eigval,eigvec = np.linalg.eig(h)
        energy = []
        for j in range(len(eigval)):
            energy.append(eigval[j].real)
        energy.sort()
        for k in range(5):
            E[k][i] += energy[k]

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i in range(5):
        plt.plot(x_axis,E[i],color='black')
        if i==0: Emin = min(E[i])
        if i==4: Emax = max(E[i])
    plt.axvline(4*math.pi/9/a)
    plt.axvline(6*math.pi/9/a)
    if band_extension:
        plt.axvline(2*math.pi/9/a*(3+math.sqrt(3)))
        plt.axvline(2*math.pi/9/a*(3+2*math.sqrt(3)))
    ax.tick_params(bottom=False)
    if band_extension:
        ax.set_xticks([0,4*math.pi/9/a,6*math.pi/9/a,2*math.pi*(3+math.sqrt(3))/9/a,2*math.pi*(3+2*math.sqrt(3))/9/a,2*math.pi*(4+2*math.sqrt(3))/9/a])
        ax.set_xticklabels(['\N{greek capital letter gamma}', 'K', 'M' , '\N{greek capital letter gamma}','M\'','K'])
    else:
        ax.set_xticks([0,4*math.pi/9/a,6*math.pi/9/a,2*math.pi*(3+math.sqrt(3))/9/a])
        ax.set_xticklabels(['\N{greek capital letter gamma}', 'K', 'M' , '\N{greek capital letter gamma}'])
    plt.savefig(filepath + imgname)
    fig.show()
    return 0


def band_plot_3d(vertexdata,H,M,phi):
    Kx=[]
    Ky=[]
    E=[[],[],[],[],[]]
    for i in range(-int(2/3*(n/2+1)),int(2/3*(n/2+1))+1,1):
        kx = 2*math.pi/3/a/(n/2+1)*i
        for j in range(-int(2/3*(m/2+1)),int(2/3*(m/2+1))+1,1):
            ky = (4*math.pi*j/3/a/(m/2+1)-kx)/math.sqrt(3)
            if abs(ky) < 2*math.sqrt(3)*math.pi/9/a and abs(ky+math.sqrt(3)*kx) < 4*math.pi/3/math.sqrt(3)/a and abs(ky-math.sqrt(3)*kx) < 4*math.pi/3/math.sqrt(3)/a:
                Kx.append(kx)
                Ky.append(ky)
                if nnn:h = fourier_2d_2(kx,ky,M,phi)
                else: fourier_2d_nnnoff(kx,ky,M,phi)
                eigval,eigvec = np.linalg.eig(h)
                energy = []
                for i in range(len(eigval)):
                    energy.append(eigval[i].real)
                energy.sort()
                for k in range(5):
                    E[k].append(energy[k].real)

    np_x = np.array(Kx)
    np_y = np.array(Ky)
    np_E = np.array(E)
    x_new, y_new = np.meshgrid(np.unique(np_x), np.unique(np_y))

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(projection='3d')
    ax.set_xlabel('kx')
    ax.set_ylabel('ky')
    ax.set_zlabel('E')
    for i in range(5):
        z_new = griddata((np_x, np_y), np_E[i], (x_new, y_new))
        ax.plot_surface(x_new,y_new,z_new)
    plt.savefig(filepath + '3D band plot')
    
    return 0

def chern_number(vertexdata,H):
    mux=2*math.pi/3/a/(n/2+1)
    muy=4*math.pi/3/math.sqrt(3)/a/(m/2+1)
    phase = np.zeros((5,4))*0j #+x,+y,-x,-y
    for i in range(5):
        for k in range(int(n/2+1)):
            for e in range(int(m/2+1)):
                if i==0:
                    if e==0: num1 = int(4*(n/2+1) + 10+(n/2+1)*m/4 + k)
                    else: num1 = int(2*n+4 + 3*(n/2+1)*(e%2) + 10*(n/2+1)*int((e-1)/2) + 3*k)
                elif i==1:
                    if e==0 and k==0: num1 = int(2*n+2)
                    elif e==0: num1 = int(2 + (k-1)*4)
                    elif e!=0 and e%2==1: num1 = int(2*n+5 + 10*(n/2+1)*int((e-1)/2) + 3*k)
                    else: num1 = int(10*(n/2+1) + 10*(n/2+1)*int((e-1)/2) + 4*k)
                elif i==2:
                    if e==0 and k==0: num1 = int(2*n+1)
                    elif e==0: num1 = int(3+(k-1)*4)
                    elif e!=0 and e%2==1: num1 = int(2*n+6 + 10*(n/2+1)*int((e-1)/2) + 3*k)
                    else: num1 = int(10*(n/2+1)+1 + 10*(n/2+1)*int((e-1)/2) + 4*k)
                elif i==3: 
                    if e==0: num1 = int(4*k)
                    elif e%2==1: num1 = int(7*(n/2+1)+1 + 10*(n/2+1)*int((e-1)/2) + 3*k)
                    else: num1 = int(10*(n/2+1)+2 + 10*(n/2+1)*int((e-1)/2) + 4*k)
                else:
                    if e==0: num1 = int(1+4*k)
                    elif e%2==1: num1 = int(7*(n/2+1)+2 + 10*(n/2+1)*int((e-1)/2) + 3*k)
                    else: num1 = int(10*(n/2+1)+3 + 10*(n/2+1)*int((e-1)/2) + 4*k)
                pos = vertexdata[num1]["pos"]
                phase[i][0]+=math.cos(mux*pos[0]) + 1j*math.sin(mux*pos[0])
                phase[i][1]+=math.cos(muy*pos[1]) + 1j*math.sin(muy*pos[1])
                phase[i][2]+=math.cos(mux*pos[0]) + 1j*math.sin(mux*pos[0])
                phase[i][3]+=math.cos(muy*pos[1]) + 1j*math.sin(muy*pos[1])
    
    
    chern_number=[0,0,0,0,0]
    for i in range(-int(2/3*(n/2+1)),int(2/3*(n/2+1))+1,1):
        kx = 2*math.pi/3/a/(n/2+1)*i
        for j in range(-int(2/3*(m/2+1)),int(2/3*(m/2+1))+1,1):
            ky = (4*math.pi*j/3/a/(m/2+1)-kx)/math.sqrt(3)
            if abs(ky) < 2*math.sqrt(3)*math.pi/9/a and abs(ky+math.sqrt(3)*kx) < 4*math.pi/3/math.sqrt(3)/a and abs(ky-math.sqrt(3)*kx) < 4*math.pi/3/math.sqrt(3)/a:
                h1 = fourier_2d(vertexdata,H,kx,ky)
                h2 = fourier_2d(vertexdata,H,kx+mux,ky)
                h3 = fourier_2d(vertexdata,H,kx+mux,ky+muy)
                h4 = fourier_2d(vertexdata,H,kx,ky+muy)
                eigval1,eigvec1 = np.linalg.eig(h1)
                eigval2,eigvec2 = np.linalg.eig(h2)
                eigval3,eigvec3 = np.linalg.eig(h3)
                eigval4,eigvec4 = np.linalg.eig(h4)
                energy = [[],[],[],[]]
                eigvec = [eigvec1,eigvec2,eigvec3,eigvec4]
                for k in range(5): 
                    energy[0].append([eigval1[k].real,k])
                    energy[1].append([eigval2[k].real,k])
                    energy[2].append([eigval3[k].real,k])
                    energy[3].append([eigval4[k].real,k])
                for k in range(4): energy[k].sort(key=lambda x:x[0])
                #if i==2 and j==2 : print(energy)
                chern_temp=[1+0j,1+0j,1+0j,1+0j,1+0j]
                for k in range(4):
                    for l in range(5): #energy
                        U=0
                        for e in range(5): #alpha inner prod
                            if k<=1: U += eigvec[k][e][energy[k][l][1]].conj() * eigvec[int((k+1)%4)][e][energy[int((k+1)%4)][l][1]] * phase[e][k] 
                            else:U += eigvec[k][e][energy[k][l][1]] * eigvec[int((k+1)%4)][e][energy[int((k+1)%4)][l][1]].conj() * phase[e][k] 
                        if k<=1: chern_temp[l] *= U/abs(U)
                        else: chern_temp[l] *= abs(U)/U
                for k in range(5):
                    z = 0
                    if abs(chern_temp[k].real) < 0.0000001:
                        z += 0.00000000001
                    else:
                        z += chern_temp[k].real
                    if abs(chern_temp[k].imag) < 0.00000001:
                        z += 0j
                    else:
                        z += chern_temp[k].imag*1j
                    chern_number[k] += cmath.log(z).imag/2/math.pi

    return chern_number

def chern_number2(vertexdata,H,N):
    k1 = [4*math.pi/9/a/N,0]
    k2 = [2*math.pi/9/a/N,2*math.pi/9/a/N/math.sqrt(3)]
    phase = np.zeros((5,4))*0j #+x,+y,-x,-y
    for i in range(5):
        for k in range(int(n/2+1)):
            for e in range(int(m/2+1)):
                if i==0:
                    if e==0: num1 = int(4*(n/2+1) + 10+(n/2+1)*m/4 + k)
                    else: num1 = int(2*n+4 + 3*(n/2+1)*(e%2) + 10*(n/2+1)*int((e-1)/2) + 3*k)
                elif i==1:
                    if e==0 and k==0: num1 = int(2*n+2)
                    elif e==0: num1 = int(2 + (k-1)*4)
                    elif e!=0 and e%2==1: num1 = int(2*n+5 + 10*(n/2+1)*int((e-1)/2) + 3*k)
                    else: num1 = int(10*(n/2+1) + 10*(n/2+1)*int((e-1)/2) + 4*k)
                elif i==2:
                    if e==0 and k==0: num1 = int(2*n+1)
                    elif e==0: num1 = int(3+(k-1)*4)
                    elif e!=0 and e%2==1: num1 = int(2*n+6 + 10*(n/2+1)*int((e-1)/2) + 3*k)
                    else: num1 = int(10*(n/2+1)+1 + 10*(n/2+1)*int((e-1)/2) + 4*k)
                elif i==3: 
                    if e==0: num1 = int(4*k)
                    elif e%2==1: num1 = int(7*(n/2+1)+1 + 10*(n/2+1)*int((e-1)/2) + 3*k)
                    else: num1 = int(10*(n/2+1)+2 + 10*(n/2+1)*int((e-1)/2) + 4*k)
                else:
                    if e==0: num1 = int(1+4*k)
                    elif e%2==1: num1 = int(7*(n/2+1)+2 + 10*(n/2+1)*int((e-1)/2) + 3*k)
                    else: num1 = int(10*(n/2+1)+3 + 10*(n/2+1)*int((e-1)/2) + 4*k)
                pos = vertexdata[num1]["pos"]
                phase[i][0]+=math.cos(k1[0]*pos[0]+k1[1]*pos[1]) + 1j*math.sin(k1[0]*pos[0]+k1[1]*pos[1])
                phase[i][1]+=math.cos(k2[0]*pos[0]+k2[1]*pos[1]) + 1j*math.sin(k2[0]*pos[0]+k2[1]*pos[1])
                phase[i][2]+=math.cos(k1[0]*pos[0]+k1[1]*pos[1]) + 1j*math.sin(k1[0]*pos[0]+k1[1]*pos[1])
                phase[i][3]+=math.cos(k2[0]*pos[0]+k2[1]*pos[1]) + 1j*math.sin(k2[0]*pos[0]+k2[1]*pos[1])

    chern_num=[0,0,0,0,0]
    for c in range(-N,N,1):
        if c < 0:
            dmin=-c-N
            dmax=N
        else:
            dmin=-N
            dmax=N-c
        for d in range(dmin,dmax,1):
            kx = 4*math.pi*d/9/a/N + 2*math.pi*c/9/a/N
            ky = 2*math.pi*c/3/math.sqrt(3)/a/N
            h1 = fourier_2d(vertexdata,H,kx,ky)
            h2 = fourier_2d(vertexdata,H,kx+k1[0],ky+k1[1])
            h3 = fourier_2d(vertexdata,H,kx+k1[0]+k2[0],ky+k1[1]+k2[1])
            h4 = fourier_2d(vertexdata,H,kx+k2[0],ky+k2[1])
            eigval1,eigvec1 = np.linalg.eig(h1)
            eigval2,eigvec2 = np.linalg.eig(h2)
            eigval3,eigvec3 = np.linalg.eig(h3)
            eigval4,eigvec4 = np.linalg.eig(h4)
            energy=[[],[],[],[]]
            eigvec = [eigvec1,eigvec2,eigvec3,eigvec4]
            for i in range(5):
                energy[0].append([eigval1[i].real,i])
                energy[1].append([eigval2[i].real,i])
                energy[2].append([eigval3[i].real,i])
                energy[3].append([eigval4[i].real,i])
            for i in range(4): energy[i].sort(key=lambda x:x[0])

            Uprod=[1+0j,1+0j,1+0j,1+0j,1+0j]
            for e in range(5):#energy
                for i in range(4):#link
                    U=0
                    for j in range(5):#alpha
                        if i<=1: U+=phase[j][i] * eigvec[i][j][energy[i][e][1]].conj() * eigvec[i+1][j][energy[i+1][e][1]]
                        else: U+=phase[j][i] * eigvec[i][j][energy[i][e][1]] * eigvec[int((i+1)%4)][j][energy[int((i+1)%4)][e][1]].conj()
                    if i<=1: Uprod[e] *= U/abs(U)
                    else: Uprod[e] *= abs(U)/U
            for e in range(5):
                chern_num[e] += cmath.log(Uprod[e]).imag/2/math.pi
    return chern_num


def chern_number3(vertexdata,H,N):
    k1 = [4*math.pi/9/a/N,0]
    k2 = [2*math.pi/9/a/N,2*math.pi/9/a/N/math.sqrt(3)]
    phase = np.zeros((5,4))*0j #+x,+y,-x,-y
    for i in range(5):
        for k in range(int(n/2+1)):
            for e in range(int(m/2+1)):
                if i==0:
                    if e==0: num1 = int(4*(n/2+1) + 10+(n/2+1)*m/4 + k)
                    else: num1 = int(2*n+4 + 3*(n/2+1)*(e%2) + 10*(n/2+1)*int((e-1)/2) + 3*k)
                elif i==1:
                    if e==0 and k==0: num1 = int(2*n+2)
                    elif e==0: num1 = int(2 + (k-1)*4)
                    elif e!=0 and e%2==1: num1 = int(2*n+5 + 10*(n/2+1)*int((e-1)/2) + 3*k)
                    else: num1 = int(10*(n/2+1) + 10*(n/2+1)*int((e-1)/2) + 4*k)
                elif i==2:
                    if e==0 and k==0: num1 = int(2*n+1)
                    elif e==0: num1 = int(3+(k-1)*4)
                    elif e!=0 and e%2==1: num1 = int(2*n+6 + 10*(n/2+1)*int((e-1)/2) + 3*k)
                    else: num1 = int(10*(n/2+1)+1 + 10*(n/2+1)*int((e-1)/2) + 4*k)
                elif i==3: 
                    if e==0: num1 = int(4*k)
                    elif e%2==1: num1 = int(7*(n/2+1)+1 + 10*(n/2+1)*int((e-1)/2) + 3*k)
                    else: num1 = int(10*(n/2+1)+2 + 10*(n/2+1)*int((e-1)/2) + 4*k)
                else:
                    if e==0: num1 = int(1+4*k)
                    elif e%2==1: num1 = int(7*(n/2+1)+2 + 10*(n/2+1)*int((e-1)/2) + 3*k)
                    else: num1 = int(10*(n/2+1)+3 + 10*(n/2+1)*int((e-1)/2) + 4*k)
                pos = vertexdata[num1]["pos"]
                phase[i][0]+=math.cos(k1[0]*pos[0]+k1[1]*pos[1]) + 1j*math.sin(k1[0]*pos[0]+k1[1]*pos[1])
                phase[i][1]+=math.cos(k2[0]*pos[0]+k2[1]*pos[1]) + 1j*math.sin(k2[0]*pos[0]+k2[1]*pos[1])
                phase[i][2]+=math.cos(k1[0]*pos[0]+k1[1]*pos[1]) + 1j*math.sin(k1[0]*pos[0]+k1[1]*pos[1])
                phase[i][3]+=math.cos(k2[0]*pos[0]+k2[1]*pos[1]) + 1j*math.sin(k2[0]*pos[0]+k2[1]*pos[1])

    phase_mat1=np.zeros((5,5))*0j
    phase_mat2=np.zeros((5,5))*0j
    """
    for i in range(5):
        for j in range(5):
            for k in range(int(n/2+1)):
                for l in range(int(n/2+1)):
                    for e in range(int(m/2+1)):
                        for f in range(int(m/2+1)):
                            if i==0:
                                if e==0: num1 = int(4*(n/2+1) + 10+(n/2+1)*m/4 + k)
                                else: num1 = int(2*n+4 + 3*(n/2+1)*(e%2) + 10*(n/2+1)*int((e-1)/2) + 3*k)
                            elif i==1:
                                if e==0 and k==0: num1 = int(2*n+2)
                                elif e==0: num1 = int(2 + (k-1)*4)
                                elif e!=0 and e%2==1: num1 = int(2*n+5 + 10*(n/2+1)*int((e-1)/2) + 3*k)
                                else: num1 = int(10*(n/2+1) + 10*(n/2+1)*int((e-1)/2) + 4*k)
                            elif i==2:
                                if e==0 and k==0: num1 = int(2*n+1)
                                elif e==0: num1 = int(3+(k-1)*4)
                                elif e!=0 and e%2==1: num1 = int(2*n+6 + 10*(n/2+1)*int((e-1)/2) + 3*k)
                                else: num1 = int(10*(n/2+1)+1 + 10*(n/2+1)*int((e-1)/2) + 4*k)
                            elif i==3: 
                                if e==0: num1 = int(4*k)
                                elif e%2==1: num1 = int(7*(n/2+1)+1 + 10*(n/2+1)*int((e-1)/2) + 3*k)
                                else: num1 = int(10*(n/2+1)+2 + 10*(n/2+1)*int((e-1)/2) + 4*k)
                            else:
                                if e==0: num1 = int(1+4*k)
                                elif e%2==1: num1 = int(7*(n/2+1)+2 + 10*(n/2+1)*int((e-1)/2) + 3*k)
                                else: num1 = int(10*(n/2+1)+3 + 10*(n/2+1)*int((e-1)/2) + 4*k)
                            if j==0:
                                if f==0: num2 = int(4*(n/2+1) + 10+(n/2+1)*m/4 + l)
                                else: num2 = int(2*n+4 + 3*(n/2+1)*(f%2) + 10*(n/2+1)*int((f-1)/2) + 3*l)
                            elif j==1:
                                if f==0 and l==0: num2 = int(2*n+2)
                                elif f==0: num2 = int(2 + (l-1)*4)
                                elif f!=0 and f%2==1: num2 = int(2*n+5 + 10*(n/2+1)*int((f-1)/2) + 3*l)
                                else: num2 = int(10*(n/2+1) + 10*(n/2+1)*int((f-1)/2) + 4*l)
                            elif j==2:
                                if f==0 and l==0: num2 = int(2*n+1)
                                elif f==0: num2 = int(3+(l-1)*4)
                                elif f!=0 and f%2==1: num2 = int(2*n+6 + 10*(n/2+1)*int((f-1)/2) + 3*l)
                                else: num2 = int(10*(n/2+1)+1 + 10*(n/2+1)*int((f-1)/2) + 4*l)
                            elif j==3: 
                                if f==0: num2 = int(4*l)
                                elif f%2==1: num2 = int(7*(n/2+1)+1 + 10*(n/2+1)*int((f-1)/2) + 3*l)
                                else: num2 = int(10*(n/2+1)+2 + 10*(n/2+1)*int((f-1)/2) + 4*l)
                            else:
                                if f==0: num2 = int(1+4*l)
                                elif f%2==1: num2 = int(7*(n/2+1)+2 + 10*(n/2+1)*int((f-1)/2) + 3*l)
                                else: num2 = int(10*(n/2+1)+3 + 10*(n/2+1)*int((f-1)/2) + 4*l)
                            pos1 = vertexdata[num1]["pos"]
                            pos2 = vertexdata[num2]["pos"]
                            phase_mat1[i][j] += math.cos(k1[0]*(pos1[0]-pos2[0])+k1[1]*(pos1[1]-pos2[1])) - 1j*math.sin(k1[0]*(pos1[0]-pos2[0])+k1[1]*(pos1[1]-pos2[1]))
                            phase_mat2[i][j] += math.cos(k2[0]*(pos1[0]-pos2[0])+k2[1]*(pos1[1]-pos2[1])) - 1j*math.sin(k2[0]*(pos1[0]-pos2[0])+k2[1]*(pos1[1]-pos2[1]))
    phase_mat = [phase_mat1,phase_mat2]
    """
    hlist = []
    #horigin = fourier_2d(vertexdata,H,-2*math.pi/9/a,-2*math.pi/3/math.sqrt(3)/a)
    for c in range(-N,N+1,1):
        if c < 0:
            dmin=-c-N
            dmax=N
        else:
            dmin=-N
            dmax=N-c
        for d in range(dmin,dmax+1,1):
            kx = 4*math.pi*d/9/a/N + 2*math.pi*c/9/a/N
            ky = 2*math.pi*c/3/math.sqrt(3)/a/N
            if d==dmin:
                hlist.append([fourier_2d_2(kx,ky,M,phi)])
            else:
                hlist[c+N].append(fourier_2d_2(kx,ky,M,phi))
            """
            if d == dmin and c==-N:
                hlist.append([horigin])
            elif d==dmin and c<=0:
                temp = np.zeros((5,5))*0j
                for i in range(5):
                    for j in range(5):
                        temp[i][j] += horigin[i][j] * (phase_mat[1][i][j] * phase_mat[0][i][j].conj())**(c+N)
                hlist.append([temp])
            elif d==dmin and c > 0:
                temp = np.zeros((5,5))*0j
                for i in range(5):
                    for j in range(5):
                        temp[i][j] += hlist[N][0][i][j] * (phase_mat[1][i][j])**(c)
                hlist.append([temp])
            else:
                temp = np.zeros((5,5))*0j
                for i in range(5):
                    for j in range(5):
                        temp[i][j] += hlist[c+N][-1][i][j] * (phase_mat[0][i][j])
                hlist[c+N].append(temp)
            """
    """
    print(len(hlist))
    for i in range(len(hlist)):
        print(len(hlist[i]))
    """
    eiglist = []
    for c in range(-N,N+1,1):
        if c < 0:
            dmin=-c-N
            dmax=N
        else:
            dmin=-N
            dmax=N-c
        for d in range(0,dmax-dmin+1,1):
            eigval,eigvec = np.linalg.eig(hlist[c+N][d])
            energy = []
            for i in range(5): energy.append([eigval[i].real,i])
            energy.sort(key=lambda x:x[0])
            if d==0: eiglist.append([[energy,eigvec]])
            else: eiglist[c+N].append([energy,eigvec])
    for i in range(N):
        eiglist[N+i+1].append(eiglist[i+1][1])
    """
    print(len(eiglist))
    for i in range(len(eiglist)):
        print(len(eiglist[i]))
    """
    chern_num = [0,0,0,0,0]
    for c in range(-N,N,1):
        if c < 0:
            dmin=-c-N
            dmax=N
        else:
            dmin=-N
            dmax=N-c
        for d in range(0,dmax-dmin,1):
            uprod=[1+0j,1+0j,1+0j,1+0j,1+0j]
            if c<0:
                clist = [0,0,1,1]
                dlist = [0,1,2,1]
            else:
                clist = [0,0,1,1]
                dlist = [0,1,1,0]
            for link in range(4):
                for Energy in range(5):
                    U = 0
                    for alpha in range(5):
                        enum1 = eiglist[c+N+clist[link]][d+dlist[link]][0][Energy][1]
                        enum2 = eiglist[c+N+clist[int((link+1)%4)]][d+dlist[int((link+1)%4)]][0][Energy][1]
                        if link <=1:
                            U += eiglist[c+N+clist[link]][d+dlist[link]][1][alpha][enum1].conj() * eiglist[c+N+clist[link+1]][d+dlist[link+1]][1][alpha][enum2] * phase[alpha][link]
                        else:
                            U += eiglist[c+N+clist[link]][d+dlist[link]][1][alpha][enum1] * eiglist[c+N+clist[int((link+1)%4)]][d+dlist[int((link+1)%4)]][1][alpha][enum2].conj() * phase[alpha][link]
                    if link <= 1: uprod[Energy] *= U/abs(U)
                    else: uprod[Energy] *= abs(U)/U
            for i in range(5):
                chern_num[i] += cmath.log(uprod[i]).imag/2/math.pi
    
    return chern_num

def chern_number4(N,filepath,M,phi):
    k1 = [4*math.pi/9/a/N,0]
    k2 = [2*math.pi/9/a/N,2*math.sqrt(3)*math.pi/9/a/N]
    Kx=[]
    Ky=[]
    chern_map = [[],[],[],[],[]]
    chern_num = [0,0,0,0,0]
    for c in range(-N,N,1):
        if c < 0:
            dmin=-c-N
            dmax=N
        else:
            dmin=-N
            dmax=N-c
        for d in range(dmin,dmax,1):
            kx = 4*math.pi*d/9/a/N + 2*math.pi*c/9/a/N
            ky = 2*math.pi*c/3/math.sqrt(3)/a/N
            kx += 1/10/N
            ky += 1/10/N
            """
            if (c==0 and d==dmax-1) or (c==-1 and d==dmax-1):
                kx -= 1/2/N
            elif (c==0 and d==dmin):
                kx += 1/2/N
            elif (c==-N and d==dmin):
                kx += 1/2/N
                ky += 1/2/N
            elif (c==N-1 and d==dmin):
                kx += 1/2/N
                ky -= 1/2/N
            elif (c==-N and d==dmax-1):
                kx -= 1/2/N
                ky += 1/2/N
            elif (c==N-1 and d==dmax-1) or (c==N-1 and d==dmax-2):
                kx -= 1/2/N
                ky -= 1/2/N
            elif (c==0 or c==-1) and (d==-1):
                kx -= 1/5/N
            elif (c==0 or c==-1) and (d==0):
                kx += 1/5/N
            """
            
            Kx.append(kx)
            Ky.append(ky)
            if nnn:
                h1 = fourier_2d_2(kx,ky,M,phi)
                h2 = fourier_2d_2(kx+k1[0],ky+k1[1],M,phi)
                h3 = fourier_2d_2(kx+k1[0]+k2[0],ky+k1[1]+k2[1],M,phi)
                h4 = fourier_2d_2(kx+k2[0],ky+k2[1],M,phi)
            else:
                h1 = fourier_2d_nnnoff(kx,ky,M,phi)
                h2 = fourier_2d_nnnoff(kx+k1[0],ky+k1[1],M,phi)
                h3 = fourier_2d_nnnoff(kx+k1[0]+k2[0],ky+k1[1]+k2[1],M,phi)
                h4 = fourier_2d_nnnoff(kx+k2[0],ky+k2[1],M,phi)
            eigval1,eigvec1 = np.linalg.eig(h1)
            eigval2,eigvec2 = np.linalg.eig(h2)
            eigval3,eigvec3 = np.linalg.eig(h3)
            eigval4,eigvec4 = np.linalg.eig(h4)
            energy=[[],[],[],[]]
            eigvec = [eigvec1,eigvec2,eigvec3,eigvec4]
            for i in range(5):
                energy[0].append([eigval1[i].real,i])
                energy[1].append([eigval2[i].real,i])
                energy[2].append([eigval3[i].real,i])
                energy[3].append([eigval4[i].real,i])
            for i in range(4): energy[i].sort(key=lambda x:x[0])

            Uprod=[1+0j,1+0j,1+0j,1+0j,1+0j]
            for e in range(5):#energy
                for i in range(4):#link
                    U=0
                    for j in range(5):#alpha
                        if i<=1: U+= eigvec[i][j][energy[i][e][1]].conj() * eigvec[i+1][j][energy[i+1][e][1]]
                        else: U+= eigvec[i][j][energy[i][e][1]] * eigvec[int((i+1)%4)][j][energy[int((i+1)%4)][e][1]].conj()
                    if i<=1: Uprod[e] *= U/abs(U)
                    else: Uprod[e] *= abs(U)/U
            for e in range(5):
                chern_num[e] += cmath.log(Uprod[e]).imag/2/math.pi
                chern_map[e].append(cmath.log(Uprod[e]).imag/2/math.pi)
    
    fig = plt.figure(figsize=(12,10))
    ax = fig.add_subplot(1,1,1)
    markermax = max(np.abs(chern_map[0]))
    mappable = ax.scatter(Kx,Ky,c=chern_map[0],cmap='bwr',s=20,alpha=1,linewidth=0.5,edgecolors='black',norm=colors.TwoSlopeNorm(0))
    fig.colorbar(mappable,ax=ax)
    fig.show()
    plt.savefig(filepath + 'berry_curv_0')

    fig = plt.figure(figsize=(12,10))
    ax = fig.add_subplot(1,1,1)
    markermax = max(np.abs(chern_map[1]))
    mappable = ax.scatter(Kx,Ky,c=chern_map[1],cmap='bwr',s=20,alpha=1,linewidth=0.5,edgecolors='black',norm=colors.TwoSlopeNorm(0))
    fig.colorbar(mappable,ax=ax)
    fig.show()
    plt.savefig(filepath + 'berry_curv_1')

    fig = plt.figure(figsize=(12,10))
    ax = fig.add_subplot(1,1,1)
    markermax = max(np.abs(chern_map[2]))
    mappable = ax.scatter(Kx,Ky,c=chern_map[2],cmap='bwr',s=20,alpha=1,linewidth=0.5,edgecolors='black',norm=colors.TwoSlopeNorm(0))
    fig.colorbar(mappable,ax=ax)
    fig.show()
    plt.savefig(filepath + 'berry_curv_2')

    fig = plt.figure(figsize=(12,10))
    ax = fig.add_subplot(1,1,1)
    markermax = max(np.abs(chern_map[3]))
    mappable = ax.scatter(Kx,Ky,c=chern_map[3],cmap='bwr',s=20,alpha=1,linewidth=0.5,edgecolors='black',norm=colors.TwoSlopeNorm(0))
    fig.colorbar(mappable,ax=ax)
    fig.show()
    plt.savefig(filepath + 'berry_curv_3')

    fig = plt.figure(figsize=(12,10))
    ax = fig.add_subplot(1,1,1)
    markermax = max(np.abs(chern_map[4]))
    mappable = ax.scatter(Kx,Ky,c=chern_map[4],cmap='bwr',s=20,alpha=1,linewidth=0.5,edgecolors='black',norm=colors.TwoSlopeNorm(0))
    fig.colorbar(mappable,ax=ax)
    fig.show()
    plt.savefig(filepath + 'berry_curv_4')
    

    return chern_num



def honeycomb_fourier(H,vertexdata,kx,ky):
    h = np.zeros((2,2))*0j
    for i in range(2):
        for j in range(2):
            for k in range(len(vertexdata)):
                for l in range(len(vertexdata)):
                    if i==0: kparity = 0
                    else: kparity = 1
                    if j==0: lparity = 0
                    else: lparity = 1
                    if vertexdata[k]["parity"] == kparity:
                        if vertexdata[l]["parity"] == lparity:
                            posk = vertexdata[k]["pos"]
                            posl = vertexdata[l]["pos"]
                            phasex = math.cos(kx*(posk[0]-posl[0])) - 1j*math.sin(kx*(posk[0]-posl[0]))
                            phasey = math.cos(ky*(posk[1]-posl[1])) - 1j*math.sin(ky*(posk[1]-posl[1]))
                            h[i][j] += H[k][l] * phasex * phasey
    return h

def haldane_chern(H,vertexdata,N):
    k1 = [2*math.pi/3/math.sqrt(3)/a/N*math.sqrt(3),2*math.pi/3/math.sqrt(3)/a/N]
    k2 = [-2*math.pi/3/math.sqrt(3)/a/N*math.sqrt(3),2*math.pi/3/math.sqrt(3)/a/N]
    phase = np.zeros((2,4))*0j

    for i in range(len(vertexdata)):
        pos = vertexdata[i]["pos"]
        if vertexdata[i]["parity"] ==0 : temp=0
        else:   temp = 1
        phase[temp][0] += math.cos(k1[0]*pos[0] + k1[1]*pos[1]) + 1j * math.sin(k1[0]*pos[0] + k1[1]*pos[1])
        phase[temp][1] += math.cos(k2[0]*pos[0] + k2[1]*pos[1]) + 1j * math.sin(k2[0]*pos[0] + k2[1]*pos[1])
        phase[temp][2] += math.cos(k1[0]*pos[0] + k1[1]*pos[1]) + 1j * math.sin(k1[0]*pos[0] + k1[1]*pos[1])
        phase[temp][3] += math.cos(k2[0]*pos[0] + k2[1]*pos[1]) + 1j * math.sin(k2[0]*pos[0] + k2[1]*pos[1])

    hmat = []
    for n2 in range(N,-N-1,-1):
        if n2 > 0:
            n1min = -N
            n1max = N-n2
        else:
            n1min = -N-n2
            n1max = N
        for n1 in range(n1min,n1max+1,1):
            kx = 2*math.pi/3/a/N*(n1+n2)
            ky = 2*math.sqrt(3)*math.pi/9/a/N*(n1-n2)
            if n1 == n1min:
                hmat.append([honeycomb_fourier2(kx,ky)])
            else:
                hmat[N-n2].append(honeycomb_fourier2(kx,ky))
    
    eigmat = []
    for n2 in range(N,-N-1,-1):
        if n2 > 0:
            n1min = -N
            n1max = N-n2
        else:
            n1min = -N-n2
            n1max = N
        for n1 in range(n1min,n1max+1,1):
            h = hmat[N-n2][n1-n1min]
            eigval,eigvec = np.linalg.eig(h)
            energy = []
            for i in range(len(eigval)):
                energy.append([eigval[i].real,i])
            energy.sort(key=lambda x:x[0])
            if n1 == n1min:
                eigmat.append([[energy,eigvec]])
            else:
                eigmat[N-n2].append([energy,eigvec])
    

    chern_num = [0,0]
    for n2 in range(N,-N,-1):
        if n2 > 0:
            n1min = -N
            n1max = N-n2
        else:
            n1min = -N-n2
            n1max = N
        for n1 in range(n1min,n1max,1): 
            uprod = [1,1]
            if n2<=0 and n1==n1min:
                n2list = [0,0,1,-N]
                n1list = [0,1,0,N-n2]
            elif n2 > 0:
                n2list = [0,0,1,1]
                n1list = [0,1,1,0]
            else:
                n2list = [0,0,1,1]
                n1list = [0,1,0,-1]
            for link in range(4):
                for ene in range(2):
                    U = 0
                    for alpha in range(2):
                        temp1 = eigmat[N-n2+n2list[link]][n1-n1min+n1list[link]][0][ene][1]
                        temp2 = eigmat[N-n2+n2list[int((link+1)%4)]][n1-n1min+n1list[int((link+1)%4)]][0][ene][1]
                        if link <= 1: U += eigmat[N-n2+n2list[link]][n1-n1min+n1list[link]][1][alpha][temp1].conj() * eigmat[N-n2+n2list[int((link+1)%4)]][n1-n1min+n1list[int((link+1)%4)]][1][alpha][temp2]
                        else: U+=eigmat[N-n2+n2list[link]][n1-n1min+n1list[link]][1][alpha][temp1] * eigmat[N-n2+n2list[int((link+1)%4)]][n1-n1min+n1list[int((link+1)%4)]][1][alpha][temp2].conj()
                    if link <= 1:
                        uprod[ene] *= U/abs(U)
                    else:
                        uprod[ene] *= abs(U)/U
            for i in range(2):
                chern_num[i] += cmath.log(uprod[i]).imag/2/math.pi

    return chern_num

def haldane_chern2(N):
    k1 = [2*math.pi/3/math.sqrt(3)/a/N*math.sqrt(3),2*math.pi/3/math.sqrt(3)/a/N]
    k2 = [-2*math.pi/3/math.sqrt(3)/a/N*math.sqrt(3),2*math.pi/3/math.sqrt(3)/a/N]
    chern_num = [0,0]
    chern_map = [[],[]]
    Kx = []
    Ky = []
    for n2 in range(N,-N,-1):
        if n2 > 0:
            n1min = -N
            n1max = N-n2
        else:
            n1min = -N-n2
            n1max = N
        for n1 in range(n1min,n1max,1):
            kx = 2*math.pi/3/a/N*(n1+n2)
            ky = 2*math.sqrt(3)*math.pi/9/a/N*(n1-n2)
            Kx.append(kx)
            Ky.append(ky)
            h1 = honeycomb_fourier2(kx,ky)
            h2 = honeycomb_fourier2(kx+k1[0],ky+k1[1])
            h3 = honeycomb_fourier2(kx+k1[0]+k2[0],ky+k1[1]+k2[1])
            h4 = honeycomb_fourier2(kx+k2[0],ky+k2[1])
            eigval1,eigvec1 = np.linalg.eig(h1)
            eigval2,eigvec2 = np.linalg.eig(h2)
            eigval3,eigvec3 = np.linalg.eig(h3)
            eigval4,eigvec4 = np.linalg.eig(h4)
            energy=[[],[],[],[]]
            eigvec = [eigvec1,eigvec2,eigvec3,eigvec4]
            for i in range(2):
                energy[0].append([eigval1[i].real,i])
                energy[1].append([eigval2[i].real,i])
                energy[2].append([eigval3[i].real,i])
                energy[3].append([eigval4[i].real,i])
            for i in range(4): energy[i].sort(key=lambda x:x[0])

            Uprod=[1+0j,1+0j]
            for e in range(2):#energy
                for i in range(4):#link
                    U=0
                    for j in range(2):#alpha
                        if i<=1: U+= eigvec[i][j][energy[i][e][1]].conj() * eigvec[i+1][j][energy[i+1][e][1]]
                        else: U+= eigvec[i][j][energy[i][e][1]] * eigvec[int((i+1)%4)][j][energy[int((i+1)%4)][e][1]].conj()
                    if i<=1: Uprod[e] *= U/abs(U)
                    else: Uprod[e] *= abs(U)/U
            for e in range(2):
                chern_num[e] += cmath.log(Uprod[e]).imag/2/math.pi
                chern_map[e].append(cmath.log(Uprod[e]).imag/2/math.pi)

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(1,1,1)
    markermax = max(np.abs(chern_map[0]))
    mappable = ax.scatter(Kx,Ky,c=chern_map[0],cmap='bwr',s=20,alpha=1,linewidth=0.5,edgecolors='black',norm=colors.TwoSlopeNorm(0))
    fig.colorbar(mappable,ax=ax)
    fig.show()
    plt.savefig(filepath + 'berry_curv_0')

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(1,1,1)
    markermax = max(np.abs(chern_map[1]))
    mappable = ax.scatter(Kx,Ky,c=chern_map[1],cmap='bwr',s=20,alpha=1,linewidth=0.5,edgecolors='black',norm=colors.TwoSlopeNorm(0))
    fig.colorbar(mappable,ax=ax)
    fig.show()
    plt.savefig(filepath + 'berry_curv_1')

    return chern_num

def haldane_band(imgname,N):
    klist = []
    x_axis=[]
    M=int(math.sqrt(3)*N)
    for i in range(M):
        klist.append([2*math.pi*i/3/a/M,0])
        x_axis.append(2*math.pi*i/3/a/M)
    for i in range(N):
        klist.append([2*math.pi/3/a , 2*math.pi*i/3/math.sqrt(3)/a/N])
        x_axis.append(2*math.pi/3/a + 2*math.pi*i/3/math.sqrt(3)/a/N)
    M=int(2*N)
    for i in range(M):
        klist.append([2*math.pi/3/a*(M-i)/M , 2*math.pi/3/a*(M-i)/M/math.sqrt(3)])
        x_axis.append(2*math.pi/3/a*(1+1/math.sqrt(3)) + (2*math.pi/3/a/math.sqrt(3) - 2*math.pi/3/a*(M-i)/M/math.sqrt(3))*2)
    klist.append([0,0])
    x_axis.append(2*math.pi/3/a/math.sqrt(3) * (3+math.sqrt(3)))
    

    E = np.zeros((2,len(klist)))*0j
    for i in range(len(klist)):
        #h = fourier_2d(vertexdata,H,klist[i][0],klist[i][1])
        h = honeycomb_fourier2(klist[i][0],klist[i][1])
        eigval,eigvec = np.linalg.eig(h)
        energy = []
        for j in range(len(eigval)):
            energy.append(eigval[j].real)
        energy.sort()
        for k in range(2):
            E[k][i] += energy[k]

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i in range(2):
        plt.plot(x_axis,E[i],color='black')
        if i==0: Emin = min(E[i])
        if i==1: Emax = max(E[i])
    plt.axvline(2*math.pi/3/a)
    plt.axvline(2*math.pi/3/a * (1+1/math.sqrt(3)))
    ax.tick_params(bottom=False)
    ax.set_xticks([0,2*math.pi/3/a,2*math.pi/3/a * (1+1/math.sqrt(3)),2*math.pi/3/a/math.sqrt(3)*(3+math.sqrt(3))])
    ax.set_xticklabels(['\N{greek capital letter gamma}', 'M', 'K' , '\N{greek capital letter gamma}'])
    plt.savefig(filepath + imgname)
    fig.show()
    return 0



def honeycomb_fourier2(kx,ky):
    h = np.zeros((2,2))*0j
    h[0][0] +=  M+t2*(4*math.cos(3*a*kx/2)*math.cos(math.sqrt(3)*a*ky/2-phi) + 2*math.cos(math.sqrt(3)*a*ky+phi))
    h[1][1] +=  -M+ t2*(4*math.cos(3*a*kx/2)*math.cos(math.sqrt(3)*a*ky/2+phi) + 2*math.cos(math.sqrt(3)*a*ky-phi))
    h[0][1] += t1*(math.cos(a*kx)+1j*math.sin(a*kx) + 2*math.cos(math.sqrt(3)*a*ky/2)*(math.cos(a*kx/2)-1j*math.sin(a*kx/2)))
    h[1][0] += h[0][1].conj()
    return h

def fourier_2d_2(kx,ky,M,phi):
    h = np.zeros((5,5))*0j
    h[0][0]+=M
    h[0][1]+=t1*(math.cos(a*kx/2+math.sqrt(3)*a*ky/2) + 1j* math.sin(a*kx/2+math.sqrt(3)*a*ky/2))
    h[0][2]+=t2*cmath.exp(1j*phi) * ( 2*cmath.exp(1j*math.sqrt(3)*a*ky/2)*math.cos(3*a*kx/2) + cmath.exp(-1j*math.sqrt(3)*a*ky) )
    h[0][3]+=t1*cmath.exp(1j*(a*kx/2 - math.sqrt(3)*a*ky/2))
    h[0][4]+=t2*cmath.exp(-1j*phi) * (cmath.exp(1j*math.sqrt(3)*a*ky) + 2*cmath.exp(-1j*math.sqrt(3)*a*ky/2) * math.cos(3*a*kx/2))
    h[1][0]+=h[0][1].conj()
    h[1][1]+=-M
    h[1][2]+=t1*cmath.exp(1j*a*kx)
    h[1][3]+=t2*cmath.exp(-1j*phi) * ( 2*cmath.exp(1j*math.sqrt(3)*a*ky/2)*math.cos(3*a*kx/2) + cmath.exp(-1j*math.sqrt(3)*a*ky) )
    h[1][4]+=h[0][3].conj()
    h[2][0]+=h[0][2].conj()
    h[2][1]+=h[1][2].conj()
    h[2][2]+=M
    h[2][3]+=h[0][1]
    h[2][4]+=h[0][2]
    h[3][0]+=h[0][3].conj()
    h[3][1]+=h[1][3].conj()
    h[3][2]+=h[2][3].conj()
    h[3][3]+=-M
    h[3][4]+=h[1][2]
    h[4][0]+=h[0][4].conj()
    h[4][1]+=h[1][4].conj()
    h[4][2]+=h[2][4].conj()
    h[4][3]+=h[3][4].conj()
    h[4][4]+=M
    return h

def fourier_2d_nnnoff(kx,ky,M,phi):
    h = np.zeros((5,5))*0j
    h[0][0]+=M
    h[0][1]+=t1*(math.cos(a*kx/2+math.sqrt(3)*a*ky/2) + 1j* math.sin(a*kx/2+math.sqrt(3)*a*ky/2))
    h[0][2]+=t2*cmath.exp(1j*phi) * (cmath.exp(1j*(3/2*a*kx+math.sqrt(3)/2*a*ky)) + cmath.exp(-1j*math.sqrt(3)*a*ky) )
    h[0][3]+=t1*cmath.exp(1j*(a*kx/2 - math.sqrt(3)*a*ky/2))
    h[0][4]+=t2*cmath.exp(-1j*phi) * (cmath.exp(1j*math.sqrt(3)*a*ky) + cmath.exp(1j*(3/2*a*kx-math.sqrt(3)/2*a*ky)))
    h[1][0]+=h[0][1].conj()
    h[1][1]+=-M
    h[1][2]+=t1*cmath.exp(1j*a*kx)
    h[1][3]+=t2*cmath.exp(-1j*phi) * ( 2*cmath.exp(1j*math.sqrt(3)*a*ky/2)*math.cos(3*a*kx/2) + cmath.exp(-1j*math.sqrt(3)*a*ky) )
    h[1][4]+=h[0][3].conj()
    h[2][0]+=h[0][2].conj()
    h[2][1]+=h[1][2].conj()
    h[2][2]+=M
    h[2][3]+=h[0][1]
    h[2][4]+=t2*cmath.exp(1j*phi) * 2 * cmath.exp(1j*math.sqrt(3)/2*a*ky) * math.cos(3/2*a*kx)
    h[3][0]+=h[0][3].conj()
    h[3][1]+=h[1][3].conj()
    h[3][2]+=h[2][3].conj()
    h[3][3]+=-M
    h[3][4]+=h[1][2]
    h[4][0]+=h[0][4].conj()
    h[4][1]+=h[1][4].conj()
    h[4][2]+=h[2][4].conj()
    h[4][3]+=h[3][4].conj()
    h[4][4]+=M
    return h

def hopping_weaken_ft(kx,ky,alpha,M,t1,t2,phi,beta):
    h = np.zeros((6,6))*0j
    h[0][0] += -M
    h[0][1] += alpha*t1*cmath.exp(1j*a*kx)
    h[0][2] += alpha*t2*cmath.exp(-1j*phi) * (2*math.cos(3/2*a*kx)*cmath.exp(1j*math.sqrt(3)/2*a*ky) + cmath.exp(-1j*math.sqrt(3)*a*ky))
    h[0][3] += alpha*t1*cmath.exp(1j*(-a/2*kx + math.sqrt(3)/2*a*ky))
    h[0][4] += alpha*t2*cmath.exp(1j*phi) * (cmath.exp(1j*math.sqrt(3)*a*ky) + 2*math.cos(3/2*a*kx) * cmath.exp(-1j*math.sqrt(3)/2*a*ky))
    h[0][5] += alpha*t1*cmath.exp(1j*(-a/2*kx-math.sqrt(3)/2*a*ky))
    h[1][0] += h[0][1].conj()
    h[1][1] += M
    h[1][2] += t1 * cmath.exp(1j*(a/2*kx + math.sqrt(3)/2*a*ky))
    h[1][3] += t2*cmath.exp(1j*phi) * ( cmath.exp(1j*(3/2*a*kx + math.sqrt(3)/2*a*ky)) + cmath.exp(-1j*math.sqrt(3)*a*ky) + beta * cmath.exp(1j*(-3/2*a*kx + math.sqrt(3)/2*a*ky)))
    h[1][4] += t1*cmath.exp(1j*(a*kx/2 - math.sqrt(3)*a*ky/2))
    h[1][5] += t2*cmath.exp(-1j*phi) * (cmath.exp(1j*math.sqrt(3)*a*ky) + cmath.exp(1j*(3/2*a*kx - math.sqrt(3)/2*a*ky)) + beta*cmath.exp(1j*(-3/2*a*kx-math.sqrt(3)/2*a*ky)))
    h[2][0] += h[0][2].conj()
    h[2][1] += h[1][2].conj()
    h[2][2] += -M
    h[2][3] += t1*cmath.exp(1j*a*kx)
    h[2][4] += t2*cmath.exp(-1j*phi) * ( 2*cmath.exp(1j*math.sqrt(3)*a*ky/2)*math.cos(3*a*kx/2) + cmath.exp(-1j*math.sqrt(3)*a*ky) )
    h[2][5] += t1*cmath.exp(1j*(-a/2*kx + math.sqrt(3)/2*a*ky))
    h[3][0] += h[0][3].conj()
    h[3][1] += h[1][3].conj()
    h[3][2] += h[2][3].conj()
    h[3][3] += M
    h[3][4] += t1*cmath.exp(1j*(1/2*a*kx + math.sqrt(3)/2*a*ky))
    h[3][5] += t2*cmath.exp(1j*phi) *(2*math.cos(3/2*a*kx) * cmath.exp(1j*math.sqrt(3)/2*a*ky) + beta* cmath.exp(-1j*math.sqrt(3)*a*ky))
    h[4][0] += h[0][4].conj()
    h[4][1] += h[1][4].conj()
    h[4][2] += h[2][4].conj()
    h[4][3] += h[3][4].conj()
    h[4][4] += -M
    h[4][5] += t1*cmath.exp(1j*a*kx)
    h[5][0] += h[0][5].conj()
    h[5][1] += h[1][5].conj()
    h[5][2] += h[2][5].conj()
    h[5][3] += h[3][5].conj()
    h[5][4] += h[4][5].conj()
    h[5][5] += M
    return h

def hopping_weaken_band(band_extension,imgname):
    klist = []
    x_axis = []
    N = 20
    for i in range(int(2*N)):
        x = 4*math.pi/9/a/2/N*i
        klist.append([x,0])
        x_axis.append(x)
    for i in range(N):
        x = 4*math.pi/9/a - (math.pi/9/a/N*i)
        klist.append([x,(4*math.pi/9/a-x)*math.sqrt(3)])
        x_axis.append(4*math.pi/9/a + (4*math.pi/9/a-x)*2)
    N3 = int(math.sqrt(3)*N)
    for i in range(N3):
        x = math.pi/3/a*(1-i/N3)
        klist.append([x,x/math.sqrt(3)])
        x_axis.append(2*math.pi/3/a+ math.pi/3/a/N3*i*2/math.sqrt(3))
    if band_extension:
        for i in range(N3):
            x = math.pi/3/a*i/N3
            klist.append([x,-x/math.sqrt(3)])
            x_axis.append(2*math.pi/a/9*(3+math.sqrt(3)) + x *2/math.sqrt(3))
        for i in range(N+1):
            x = math.pi/3/a + math.pi/9/a*i/N
            klist.append([x,-(4*math.pi/9/a - x)*math.sqrt(3)])
            x_axis.append(2*math.pi/a/9*(3+2*math.sqrt(3)) + math.pi/9/a*i/N*2 )
    else:
        klist.append([0,0])
        x_axis.append(2*math.pi/9/a*(3+math.sqrt(3)))

    E = np.zeros((6,len(klist)))*0j
    for i in range(len(klist)):
        h = hopping_weaken_ft(klist[i][0],klist[i][1],alpha,M,t1,t2,phi,beta)
        #if nnn:h = fourier_2d_2(klist[i][0],klist[i][1],M,phi)
        #else: h = fourier_2d_nnnoff(klist[i][0],klist[i][1],M,phi)
        eigval,eigvec = np.linalg.eig(h)
        energy = []
        for j in range(len(eigval)):
            energy.append(eigval[j].real)
        energy.sort()
        for k in range(6):
            E[k][i] += energy[k]

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i in range(6):
        plt.plot(x_axis,E[i],color='black')
        if i==0: Emin = min(E[i])
        if i==4: Emax = max(E[i])
    plt.axvline(4*math.pi/9/a)
    plt.axvline(6*math.pi/9/a)
    if band_extension:
        plt.axvline(2*math.pi/9/a*(3+math.sqrt(3)))
        plt.axvline(2*math.pi/9/a*(3+2*math.sqrt(3)))
    ax.tick_params(bottom=False)
    if band_extension:
        ax.set_xticks([0,4*math.pi/9/a,6*math.pi/9/a,2*math.pi*(3+math.sqrt(3))/9/a,2*math.pi*(3+2*math.sqrt(3))/9/a,2*math.pi*(4+2*math.sqrt(3))/9/a])
        ax.set_xticklabels(['\N{greek capital letter gamma}', 'K', 'M' , '\N{greek capital letter gamma}','M\'','K'])
    else:
        ax.set_xticks([0,4*math.pi/9/a,6*math.pi/9/a,2*math.pi*(3+math.sqrt(3))/9/a])
        ax.set_xticklabels(['\N{greek capital letter gamma}', 'K', 'M' , '\N{greek capital letter gamma}'])
    plt.savefig(filepath + imgname)
    fig.show()
    return 0

def hopping_weaken_calc_Uprod(k1,k2,kx,ky,alpha,M,t1,t2,phi,beta):
    h1 = hopping_weaken_ft(kx,ky,alpha,M,t1,t2,phi,beta)
    h2 = hopping_weaken_ft(kx+k1[0],ky+k1[1],alpha,M,t1,t2,phi,beta)
    h3 = hopping_weaken_ft(kx+k1[0]+k2[0],ky+k1[1]+k2[1],alpha,M,t1,t2,phi,beta)
    h4 = hopping_weaken_ft(kx+k2[0],ky+k2[1],alpha,M,t1,t2,phi,beta)

    eigval1,eigvec1 = np.linalg.eig(h1)
    eigval2,eigvec2 = np.linalg.eig(h2)
    eigval3,eigvec3 = np.linalg.eig(h3)
    eigval4,eigvec4 = np.linalg.eig(h4)
    energy=[[],[],[],[]]
    eigvec = [eigvec1,eigvec2,eigvec3,eigvec4]
    for i in range(len(h1[0])):
        energy[0].append([eigval1[i].real,i])
        energy[1].append([eigval2[i].real,i])
        energy[2].append([eigval3[i].real,i])
        energy[3].append([eigval4[i].real,i])
    for i in range(4): energy[i].sort(key=lambda x:x[0])

    Uprod=[1+0j,1+0j,1+0j,1+0j,1+0j,1+0j]
    for e in range(6):#energy
        for i in range(4):#link
            U=0
            for j in range(6):#alpha
                if i<=1: U+= eigvec[i][j][energy[i][e][1]].conj() * eigvec[i+1][j][energy[i+1][e][1]]
                else: U+= eigvec[i][j][energy[i][e][1]] * eigvec[int((i+1)%4)][j][energy[int((i+1)%4)][e][1]].conj()
            if U==0:
                Uprod[e] *= 0
            else:
                if i<=1: Uprod[e] *= U/abs(U)
                else: Uprod[e] *= abs(U)/U
    return Uprod
    

def hopping_weaken_chern(N,filepath,M,phi,alpha,t1,t2,beta):
    k1 = [4*math.pi/9/a/N,0]
    k2 = [2*math.pi/9/a/N,2*math.sqrt(3)*math.pi/9/a/N]
    Kx=[]
    Ky=[]
    chern_map = [[],[],[],[],[],[]]
    chern_num = [0,0,0,0,0,0]
    for c in range(-N,N,1):
        if c < 0:
            dmin=-c-N
            dmax=N
        else:
            dmin=-N
            dmax=N-c
        for d in range(dmin,dmax,1):
            kx = 4*math.pi*d/9/a/N + 2*math.pi*c/9/a/N
            ky = 2*math.pi*c/3/math.sqrt(3)/a/N
            kx += 1/5/N
            ky += 1/5/N
            """
            if (c==0 and d==dmax-1) or (c==-1 and d==dmax-1):
                kx -= 1/2/N
            elif (c==0 and d==dmin):
                kx += 1/2/N
            elif (c==-N and d==dmin):
                kx += 1/2/N
                ky += 1/2/N
            elif (c==N-1 and d==dmin):
                kx += 1/2/N
                ky -= 1/2/N
            elif (c==-N and d==dmax-1):
                kx -= 1/2/N
                ky += 1/2/N
            elif (c==N-1 and d==dmax-1) or (c==N-1 and d==dmax-2):
                kx -= 1/2/N
                ky -= 1/2/N
            elif (c==0 or c==-1) and (d==-1):
                kx -= 1/5/N
            elif (c==0 or c==-1) and (d==0):
                kx += 1/5/N
            """
            flag = True
            while(flag):
                flag2 = True
                Uprod = hopping_weaken_calc_Uprod(k1,k2,kx,ky,alpha,M,t1,t2,phi,beta)
                for i in range(6):
                    if Uprod[i]==0:
                        kx += 1/10/N
                        flag2=False
                        break
                if flag2: flag = False
                        
            Kx.append(kx)
            Ky.append(ky)

            for e in range(6):
                chern_num[e] += cmath.log(Uprod[e]).imag/2/math.pi
                chern_map[e].append(cmath.log(Uprod[e]).imag/2/math.pi)
    
    for i in range(6):
        fig = plt.figure(figsize=(12,10))
        ax = fig.add_subplot(1,1,1)
        markermax = max(np.abs(chern_map[i]))
        mappable = ax.scatter(Kx,Ky,c=chern_map[i],cmap='bwr',s=20,alpha=1,linewidth=0.5,edgecolors='black',norm=colors.TwoSlopeNorm(0))
        fig.colorbar(mappable,ax=ax)
        fig.show()
        plt.savefig(filepath + 'berry_curv_' + str(i))

    return chern_num

if __name__ == "__main__":
    main()