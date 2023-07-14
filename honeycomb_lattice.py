import matplotlib.pyplot as plt
import math
import numpy as np


def vertex_create(a,n,m,xPBC,yPBC):
    vertexdata = []
    plaquettedata = []
    plaquettenum = 0
    N = 2*(n+1)*(m+1)
    for i in range(m+1):
        vertexmemo = 0
        for j in range(3*n+2):
            if j % 2 == 0: 
                pos = np.array([j*a/2,math.sqrt(3)*a*i])
            else : 
                pos = np.array([j*a/2,math.sqrt(3)*a*(i+0.5)])
            if j % 3 == 2: 
                neighbor = []
                dn_odd = [-1,-n,1,n-1,n,n+1]
                dn_even = [-n-1,-n,-n+1,-1,n,1]
                if (j-2)/3 % 2 == 1:
                    dn = dn_odd
                else :
                    dn = dn_even
                for k in range(6):
                    if 0 <= plaquettenum+dn[k] and plaquettenum + dn[k] < n*(m+1):
                        if (plaquettenum%n==0 and (dn[k]==-1 or dn[k]==-n-1)) or (plaquettenum%n==n-1 and (dn[k]==1 or dn[k]==-n+1 or dn[k]==n+1)):
                            continue
                        else:
                            neighbor.append(plaquettenum + dn[k])
                plaquettenum += 1
                plaquettedata.append({"pos":pos , "neighbor":neighbor})
            else:
                neighbor = []
                secondneighbor = []
                dn_0 = [-1,1,-2*n-1]
                dn_1 = [-1,1,2*n+1]
                dn_2 = [-1,1,2*n+3]
                dn_3 = [-1,1,-2*n-3]
                if vertexmemo %4 == 0: 
                    dn = dn_0
                elif vertexmemo %4 == 1:
                    dn = dn_1
                elif vertexmemo %4 == 2:
                    dn = dn_2
                else:
                    dn = dn_3
                for k in range(3):
                    if 0 <= i*2*(n+1)+vertexmemo+dn[k] and i*2*(n+1)+vertexmemo+dn[k] < 2*(n+1)*(m+1) :
                        if (vertexmemo == 0 and dn[k]==-1) or (vertexmemo == 2*n+1 and dn[k] == 1):
                            continue
                        else:
                            neighbor.append(i*2*(n+1)+vertexmemo+dn[k])
                if yPBC and i==0: 
                    if vertexmemo%4 == 0: neighbor.append(2*(n+1)*m+vertexmemo+1)
                    elif vertexmemo%4==3: neighbor.append(2*(n+1)*m+vertexmemo-1)
                if yPBC and i == m:
                    if vertexmemo%4==1: neighbor.append(vertexmemo-1)
                    elif vertexmemo%4==2: neighbor.append(vertexmemo+1)
                if xPBC and n%2==0: print("set n as odd number")
                elif xPBC and n%2==1:
                    if vertexmemo==0: neighbor.append(i*2*(n+1)+2*n+1)
                    elif vertexmemo==2*n+1: neighbor.append(i*2*(n+1))
                    
                #secondneighbor
                dn_12=[-2,2,2*n+2,-2*n-2,2*n,2*n+4]
                dn_03=[-2*n-4,-2*n,2*n+2,-2,2,-2*n-2]
                if vertexmemo%4==0 or vertexmemo%4==3:
                    dnsecond = dn_03
                else:
                    dnsecond = dn_12
                for k in range(6):
                    if 0 <= i*2*(n+1)+vertexmemo+dnsecond[k] and i*2*(n+1)+vertexmemo+dnsecond[k] < 2*(n+1)*(m+1):
                        if (vertexmemo < 2 and (dnsecond[k]==-2 or dnsecond[k]==2*n or dnsecond[k]==-2*n-4)) or (vertexmemo >=2*n and(dnsecond[k]==2 or dnsecond[k]==2*n+4 or dnsecond[k]==-2*n)):
                            continue
                        else:
                            if k < 3:
                                secondneighbor.append([i*2*(n+1)+vertexmemo+dnsecond[k],1])
                            else:
                                secondneighbor.append([i*2*(n+1)+vertexmemo+dnsecond[k],-1])
                if yPBC and i==0:
                    if vertexmemo%4==1 or vertexmemo%4==2: secondneighbor.append([2*(n+1)*m+vertexmemo,-1])
                    elif vertexmemo%4==0 or vertexmemo%4==3:
                        secondneighbor.append([2*(n+1)*m+vertexmemo,-1])
                        if j==0: secondneighbor.append([2*(n+1)*m+vertexmemo+2,1])
                        elif vertexmemo >= 2*n: secondneighbor.append([2*(n+1)*m+vertexmemo-2,1])
                        else: 
                            secondneighbor.append([2*(n+1)*m+vertexmemo+2,1])
                            secondneighbor.append([2*(n+1)*m+vertexmemo-2,1])
                if yPBC and i==m:
                    if vertexmemo%4==0 or vertexmemo%4==3: secondneighbor.append([vertexmemo,1])
                    elif vertexmemo%4== 1 or vertexmemo%4==2:
                        secondneighbor.append([vertexmemo,1])
                        if j==1: secondneighbor.append([vertexmemo+2,-1])
                        elif j>=3*n: secondneighbor.append([vertexmemo-2,-1])
                        else:
                            secondneighbor.append([vertexmemo+2,-1])
                            secondneighbor.append([vertexmemo-2,-1])
                if xPBC and n%2==1:
                    if vertexmemo == 0:
                        secondneighbor.append([i*2*(n+1)+2*n,-1])
                        if i != 0:
                            secondneighbor.append([i*2*(n+1)-2,1])
                    if vertexmemo == 1:
                        secondneighbor.append([i*2*(n+1)+2*n+1,1])
                        if i != m:
                            secondneighbor.append([(i+1)*2*(n+1)+2*n+1,-1])
                    if vertexmemo == 2*n+1:
                        secondneighbor.append([i*2*(n+1)+1,-1])
                        if i != 0:
                            secondneighbor.append([(i-1)*2*(n+1)+1,1])
                    if vertexmemo == 2*n:
                        secondneighbor.append([i*2*(n+1),1])
                        if i != m:
                            secondneighbor.append([(i+1)*2*(n+1),-1])

                vertexdata.append({"pos":pos , "neighbor":neighbor,"secondneighbor":secondneighbor})
                vertexmemo += 1
    return [vertexdata,plaquettedata]

def plot(n,m,pointdata,filepath,plaq,second):
    vertexdata = pointdata[0]
    plaquettedata = pointdata[1]
    savename = filepath + str(n) + "times" + str(m) + "honeycomb" 
    fig = plt.figure(figsize=(n,m))
    ax = fig.add_subplot(1,1,1)
    for i in range(len(vertexdata)):
        for j in range(len(vertexdata[i]["neighbor"])):
            neinum = vertexdata[i]["neighbor"][j]
            nowpos = vertexdata[i]["pos"]
            neipos = vertexdata[neinum]["pos"]
            plt.plot([nowpos[0],neipos[0]],[nowpos[1],neipos[1]],color="black")
    if plaq:
        for i in range(len(plaquettedata)):
            for j in range(len(plaquettedata[i]["neighbor"])):
                neinum = plaquettedata[i]["neighbor"][j]
                neipos = plaquettedata[neinum]["pos"]
                nowpos = plaquettedata[i]["pos"]
                plt.plot([nowpos[0],neipos[0]],[nowpos[1],neipos[1]],color="orange",linestyle="dashed",linewidth=0.5)
    if second:
        for i in range(len(vertexdata)):
            for j in range(len(vertexdata[i]["secondneighbor"])):
                nnn = vertexdata[i]["secondneighbor"][j]
                nnnpos = vertexdata[nnn[0]]["pos"]
                nowpos = vertexdata[i]["pos"]
                if nnn[1] == 1:
                    plt.plot([nowpos[0],nnnpos[0]],[nowpos[1],nnnpos[1]],color="green",linestyle="dashed",linewidth=1)
                else:
                    plt.plot([nowpos[0],nnnpos[0]],[nowpos[1],nnnpos[1]],color="green",linestyle = 'dashed',linewidth=1)
    ax.set_aspect('equal')
    fig.show()
    plt.savefig(savename)


def main():
    #lattice constant
    a = 1
    #system size 
    n = 5
    m = 4
    #plot the plaquette center and bond:True otherwise:False
    plaq = False
    second = True
    xPBC = True
    yPBC = False
    filepath = "/Users/sogenikegami/Documents/UT4S/non-crystal/honeycomb/image/"

    pointdata = vertex_create(a,n,m,xPBC,yPBC)
    print(pointdata[0][31]["secondneighbor"])
    plot(n,m,pointdata,filepath,plaq,second)

    return pointdata




if __name__ == "__main__":
    main()