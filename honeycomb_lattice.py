import matplotlib.pyplot as plt
import math
import numpy as np


def vertex_create(a,n,m):
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
                        neighbor.append(plaquettenum + dn[k])
                plaquettenum += 1
                plaquettedata.append({"pos":pos , "neighbor":neighbor})
            else:
                neighbor = []
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
                vertexdata.append({"pos":pos , "neighbor":neighbor})
                vertexmemo += 1
    return [vertexdata,plaquettedata]

def plot(n,m,pointdata,filepath):
    vertexdata = pointdata[0]
    plaquettedata = pointdata[1]
    savename = filepath + str(n) + "times" + str(m) + "honeycomb" 
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i in range(len(vertexdata)):
        for j in range(len(vertexdata[i]["neighbor"])):
            neinum = vertexdata[i]["neighbor"][j]
            nowpos = vertexdata[i]["pos"]
            neipos = vertexdata[neinum]["pos"]
            plt.plot([nowpos[0],neipos[0]],[nowpos[1],neipos[1]],color="blue")
    ax.set_aspect('equal')
    fig.show()
    plt.savefig(savename)


def main():
    #lattice constant
    a = 1
    #system size 
    n = 6
    m = 4
    filepath = "/Users/sogenikegami/Documents/UT4S/non-crystal/honeycomb/"

    pointdata = vertex_create(a,n,m)
    plot(n,m,pointdata,filepath)
    fig= plt.figure()
    ax = fig.add_subplot(1,1,1)
    for i in range(len(pointdata[0])):
        ax.scatter(pointdata[0][i]["pos"][0],pointdata[0][i]["pos"][1])
    fig.show()
    plt.savefig("/Users/sogenikegami/Documents/UT4S/non-crystal/honeycomb/test")
    print(int(13/14))



if __name__ == "__main__":
    main()