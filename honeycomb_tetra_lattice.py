import math
import numpy as np
import matplotlib.pyplot as plt
import random

def main():
    a=1
    n=12
    m=12
    xPBC=False
    yPBC=False
    mode = 0
    defectratio = 0
    secondplot=True
    nnn = True
    vertexdata = vertex_create(a,n,m,xPBC,yPBC,mode,defectratio,nnn)
    lattice_plot(vertexdata,n,m,a,xPBC,yPBC,secondplot)
    print(len(vertexdata))
    #print(vertexdata[0]["secondneighbor"][2][0])
    print(vertexdata[0])
    return 0

def vertex_create(a,n,m,xPBC,yPBC,mode,defectratio,nnn):
    vertexdata1 = position(n,m,a,xPBC,yPBC)
    if xPBC: vertexdata2 = neighbor_xPBC(vertexdata1,n,m,yPBC)
    else: vertexdata2 = neighbor_xopen(vertexdata1,n,m,yPBC)
    vertexdata3 = secondneighbor(vertexdata2,n,m,a,xPBC,yPBC)
    defectnum = defectnum_gen(vertexdata3,defectratio,xPBC,yPBC,n,m,mode)
    vertexdata4 = defect(vertexdata3,defectnum,nnn)
    #vertexdata[0] = {"pos"=[x,y],"neighbor"=[],"secondneighbor"=[],"parity"=[]}
    return vertexdata4

def position(n,m,a,xPBC,yPBC):
    N = (2*n+1) + (m-2)*(4*n+5)
    count = 0
    vertexdata=[]
    if xPBC:
        bottom=3
        bulk1=1
        bulk2=2
        top=3
    else:
        bottom=0
        bulk1=0
        bulk2=0
        if yPBC:
            top=2
        else:   #both open
            top = 0

    for i in range(2*n+1+bottom):
        if i%4==0: pos = [a/2+3*i*a/4,math.sqrt(3)/2*a]
        elif i%4==1: pos = [a+a/2+3*(i-1)*a/4,math.sqrt(3)/2*a]
        elif i%4==2: pos = [2*a+3*(i-2)*a/4,0]
        else:  pos = [3*a+3*(i-3)*a/4,0]
        vertexdata.append({"pos":pos})
        if i==0 or i==2*n+1*bottom-1 or i%4==2 or i%4==3:
            vertexdata[count]["edge"] = True
        else: vertexdata[count]["edge"] = False
        if i%2==0: vertexdata[count]["parity"] = 0
        else: vertexdata[count]["parity"] = 1
        count+=1

    for j in range(int(m/4)*3-1):
        if j%3==0:
            for i in range(2*n+3+bulk1):
                if i%4==0: pos = [a*j+3*i*a/4,(j+1)*math.sqrt(3)*a]
                elif i%4==1: pos = [a/2+a*j+3*(i-1)*a/4,(j+3/2)*math.sqrt(3)*a]
                elif i%4==2: pos = [a*3/2+a*j+3*(i-2)*a/4,(j+3/2)*math.sqrt(3)*a]
                else:  pos = [a*2+a*j+3*(i-3)*a/4,(j+1)*math.sqrt(3)*a]
                vertexdata.append({"pos":pos})
                if i==0 or i==1 or i==2*n+3*bulk1-1: vertexdata[count]["edge"] = True
                else: vertexdata[count]["edge"] = False
                if i%2==0: vertexdata[count]["parity"] = 1
                else: vertexdata[count]["parity"] = 0
                count+=1
        elif j%3==1:
            for i in range(2*n+2+bulk2):
                if i%4==0: pos = [a*3/2+(j-1)*a+3*i*a/4,math.sqrt(3)*(j+3/2)*a]
                elif i%4==1: pos = [a*2+(j-1)*a+3*(i-1)*a/4,math.sqrt(3)*(j+1)*a]
                elif i%4==2: pos = [a*3+(j-1)*a+3*(i-2)*a/4,math.sqrt(3)*(j+1)*a]
                else:  pos = [a*7/2+(j-1)*a+3*(i-3)*a/4,math.sqrt(3)*(j+3/2)*a]
                vertexdata.append({"pos":pos})
                if i==0 or i==2*n+2+bulk2-1: vertexdata[count]["edge"] = True
                else: vertexdata[count]["edge"] = False
                if i%2==0: vertexdata[count]["parity"] = 1
                else: vertexdata[count]["parity"] = 0
                count+=1
        elif j%3==2:
            for i in range(2*n+3+bulk1):
                if i%4==0: pos = [a*2+(j-2)*a+3*i*a/4,math.sqrt(3)*(j+1)*a]
                elif i%4==1: pos = [a*3+(j-2)*a+3*(i-1)*a/4,math.sqrt(3)*(j+1)*a]
                elif i%4==2: pos = [a*7/2+(j-2)*a+3*(i-2)*a/4,math.sqrt(3)*(j+3/2)*a]
                else:  pos = [a*9/2+(j-2)*a+3*(i-3)*a/4,math.sqrt(3)*(j+3/2)*a]
                vertexdata.append({"pos":pos})
                if xPBC:
                    if i==0 or i==2*n+3+bulk1-1: vertexdata[count]["edge"] = True
                    else: vertexdata[count]["edge"] = False
                else:
                    if i==0 or i==2*n+3+bulk1-1 or i==2*n+3+bulk1-2: vertexdata[count]["edge"] = True
                    else: vertexdata[count]["edge"] = False
                if i%2==0: vertexdata[count]["parity"] = 0
                else: vertexdata[count]["parity"] = 1
                count += 1
    for i in range(2*n+1+top):
        if i%4==0: pos = [a*2+(m/4*3-3)*a+3*i*a/4,math.sqrt(3)*(m/4*3)*a]
        elif i%4==1: pos = [a*3+(m/4*3-3)*a+3*(i-1)*a/4,math.sqrt(3)*(m/4*3)*a]
        elif i%4==2: pos = [a*7/2+(m/4*3-3)*a+3*(i-2)*a/4,math.sqrt(3)*(m/4*3+1/2)*a]
        else:  pos = [a*9/2+(m/4*3-3)*a+3*(i-3)*a/4,math.sqrt(3)*(m/4*3+1/2)*a]
        vertexdata.append({"pos":pos})
        if yPBC:
            if i==0 or i==2*n+top or (xPBC==False and i==2*n+top-1): vertexdata[count]["edge"] = True
            else: vertexdata[count]["edge"] = False
        else: 
            if i==0 or i==2*n+top or i%4==2 or i%4==3: vertexdata[count]["edge"] =True
            else: vertexdata[count]["edge"] = False
        if i%2==0: vertexdata[count]["parity"] = 0
        else: vertexdata[count]["parity"] = 1
        count+=1

    if yPBC:
        vertexdata.append({"pos":[a*m/4*3,(1+m/4*3)*math.sqrt(3)*a]})
        vertexdata[count]["edge"] = True
        vertexdata[count]["parity"] = 1
        count+=1
        if xPBC:
            edge=1
        else:
            vertexdata.append({"pos":[a*(m/4*3+1/2),(3/2+m/4*3)*math.sqrt(3)*a]})
            vertexdata[count]["edge"]=True
            vertexdata[count]["parity"] = 0
            count+=1
            vertexdata.append({"pos":[a*(m/4*3+3/2),(3/2+m/4*3)*math.sqrt(3)*a]})
            vertexdata[count]["edge"] = True
            vertexdata[count]["parity"] = 1
            count += 1
            edge=0
        for i in range(n+edge):
            if i%2==0:
                pos = [a*(m/4*3+2)+3*a*i/2,(1+m/4*3)*math.sqrt(3)*a]
            else:
                pos = [a*(m/4*3+3)+3*a*(i-1)/2,(1+m/4*3)*math.sqrt(3)*a]
            vertexdata.append({"pos":pos})
            if i==0: vertexdata[count]["edge"] = False
            else: vertexdata[count]["edge"] = True
            if i%2== 0: vertexdata[count]["parity"] = 0
            else: vertexdata[count]["parity"] = 1
            count+=1 
    return vertexdata



def neighbor_xopen(vertexdata,n,m,yPBC):
    count = 0
    if yPBC:
        N = int((6*n+8)*(m/4)+2*n+1+2)
        for i in range(2*n+1):
            if i==0: neighbor = [count+1,count+2*n+1,N]
            elif i== 2*n-1: neighbor = [count-1,count+1]
            elif i==2*n: neighbor = [count-1,count+2*n+1]
            else:
                if i%4==0: neighbor = [count-1,count+1,count+2*n+1]
                elif i%4==1: neighbor = [count-1,count+1,count+2*n+3]
                elif i%4==2: neighbor = [count-1,count+1,N+int(i/2)+1]
                else: neighbor = [count-1,count+1,N+int(i/2)+2]
            vertexdata[count]["neighbor"] = neighbor
            count+=1
    else:
        for i in range(2*n+1):
            if i==0: neighbor=[count+1,count+2*n+1]
            elif i==2*n: neighbor=[count-1,count+2*n+1]
            else:
                if i%4==0: neighbor = [count-1,count+1,count+2*n+1]
                elif i%4==1: neighbor = [count-1,count+1,count+2*n+3]
                else: neighbor = [count-1,count+1]
            vertexdata[count]["neighbor"] = neighbor
            count += 1

    for j in range(int(m/4)*3-1):
        if j%3==0:
            for i in range(2*n+3):
                if i==0: neighbor = [count+1,count-(2*n+1)]
                elif i==1:neighbor = [count-1,count+1]
                elif i==2*n+2: neighbor = [count-1,count+2*n+2]
                else:
                    if i%4==0: neighbor = [count-(2*n+1),count-1,count+1]
                    elif i%4==1: neighbor = [count-1,count+1,count+2*n]
                    elif i%4==2: neighbor = [count-1,count+1,count+2*n+2]
                    else: neighbor = [count-(2*n+3),count-1,count+1]
                vertexdata[count]["neighbor"] = neighbor
                count += 1
        elif j%3==1:
            for i in range(2*n+2):
                if i==0: neighbor = [count+1,count+2*n+2]
                elif i==2*n+1: neighbor = [count-(2*n+2),count-1]
                else:
                    if i%4==0: neighbor = [count-1,count+1,count+(2*n+2)]
                    elif i%4==1: neighbor = [count-(2*n+2),count-1,count+1]
                    elif i%4==2: neighbor = [count-2*n,count-1,count+1]
                    else: neighbor = [count-1,count+1,count+2*n]
                vertexdata[count]["neighbor"] = neighbor
                count+=1
        else:
            for i in range(2*n+3):
                if i==0: neighbor = [count-(2*n+2),count+1]
                elif i==2*n+2: neighbor = [count-1,count+2*n+1]
                elif i==2*n+1: neighbor = [count-1,count+1]
                else:
                    if i%4==0: neighbor = [count-(2*n+2),count-1,count+1]
                    elif i%4==1: neighbor = [count-2*n,count-1,count+1]
                    elif i%4==2: neighbor = [count-1,count+1,count+2*n+1]
                    else: neighbor = [count-1,count+1,count+2*n+3]
                vertexdata[count]["neighbor"] = neighbor
                count += 1
    if yPBC:
        for i in range(2*n+3):
            if i==0: neighbor = [count-(2*n+2),count+1]
            elif i==2*n+2: neighbor = [count-1,count+n+3]
            elif i==2*n+1: neighbor = [count-1,count+1]
            elif i==2: neighbor = [count-1,count+1,count+2*n+1]
            else:
                if i%4==0: neighbor = [count-(2*n+2),count-1,count+1]
                elif i%4==1: neighbor = [count-2*n,count-1,count+1]
                else: neighbor = [count-1,count+1,count+2*n+3-int(i/4)*2]
            vertexdata[count]["neighbor"] = neighbor
            count+=1
        for i in range(n+3):
            if i==0: neighbor = [count-(2*n+1),count+1]
            elif i==1: neighbor = [count-1,count+1]
            elif i==2: neighbor = [0,count-1,count+1]
            elif i==3: neighbor = [count-(2*n+3),count-1,count+1]
            elif i==n+2: neighbor = [count-1,2*n-2,count-(n+3)]
            else:
                if i%2==0: neighbor = [2*i-6,count-1,count-(2*n+5-int(i/2)*2)]
                else: neighbor = [count+1,2*i-7,count-(2*n+5-int(i/2)*2)]
            vertexdata[count]["neighbor"] = neighbor
            count+= 1
    else:
        for i in range(2*n+1):
            if i==0: neighbor = [count-(2*n+2),count+1]
            elif i==2*n: neighbor = [count-(2*n+2),count-1]
            else:
                if i%4==0: neighbor = [count-(2*n+2),count-1,count+1]
                elif i%4==1: neighbor = [count-2*n,count-1,count+1]
                else: neighbor = [count-1,count+1]
            vertexdata[count]["neighbor"] = neighbor
            count+= 1

    return vertexdata

def neighbor_xPBC(vertexdata,n,m,yPBC):
    count=0

    if yPBC:
        N = int((2*n+4)*(m/4*3+1))
        for i in range(2*n+4):
            if i==0: neighbor = [count+1,count+2*n+4,count+2*n+3]
            elif i==2*n+3: neighbor = [0,count-1,N+1]
            elif i==2*n+2: neighbor = [count+1,count-1,N]
            else:
                if i%4==0: neighbor = [count-1,count+1,count+2*n+4]
                elif i%4==1: neighbor = [count-1,count+1,count+2*n+6]
                elif i%4==2: neighbor = [count-1,count+1,N+int(i/4)*2+2]
                else: neighbor = [count-1,count+1,N+int(i/4)*2+3]
            vertexdata[count]["neighbor"] = neighbor
            count+=1
    else:
        for i in range(2*n+4):
            if i==0: neighbor = [count+1,count+2*n+4,count+2*n+3]
            elif i== 2*n+3: neighbor = [0,count-1]
            else:
                if i%4==0: neighbor = [count-1,count+1,count+2*n+4]
                elif i%4==1: neighbor = [count-1,count+1,count+2*n+6]
                else: neighbor = [count-1,count+1]
            vertexdata[count]["neighbor"] = neighbor
            count += 1

    for j in range(int(m/4)*3-1):
        if j%3==0 and j!=0:
            for i in range(2*(n+2)):
                if i==0 and j==0: neighbor = [0,count+1,count+2*n+3]
                elif i==0: neighbor = [count-(2*n+2),count+1,count+2*n+3]
                elif i==1: neighbor = [count-1,count+1,count+(4*n+5)]
                elif i==2*n+3 and j==0: neighbor = [2*n+1,count-1,count-(2*n+3)]
                elif i==2*n+3: neighbor = [count-(2*n+4),count-1,count-(2*n+3)]
                else:
                    if i%4==0: neighbor = [count-(2*n+2),count-1,count+1]
                    elif i%4==1: neighbor = [count+2*n+1,count-1,count+1]
                    elif i%4==2: neighbor = [count+2*n+3,count-1,count+1]
                    else: neighbor = [count-(2*n+4),count-1,count+1]
                vertexdata[count]["neighbor"] = neighbor
                count += 1
        elif j==0:
            for i in range(2*n+4):
                if i==0 and j==0: neighbor = [0,count+1,count+2*n+3]
                elif i==0: neighbor = [count-(2*n+2),count+1,count+2*n+3]
                elif i==1: neighbor = [count-1,count+1,count+(4*n+5)]
                elif i==2*n+3 and j==0: neighbor = [2*n+1,count-1,count-(2*n+3)]
                elif i==2*n+3: neighbor = [count-(2*n+4),count-1,count-(2*n+3)]
                else:
                    if i%4==0: neighbor = [count-(2*n+4),count-1,count+1]
                    elif i%4==1: neighbor = [count+2*n+1,count-1,count+1]
                    elif i%4==2: neighbor = [count+2*n+3,count-1,count+1]
                    else: neighbor = [count-(2*n+6),count-1,count+1]
                vertexdata[count]["neighbor"] = neighbor
                count += 1
        elif j%3==1:
            for i in range(2*(n+2)):
                if i==0: neighbor = [count+2*n+4,count+1,count+2*n+3]
                elif i==2*n+3: neighbor = [count-1,count-(2*n+3),count+2*n+2]
                elif i==2*n+2: neighbor = [count-1,count+1,count-(4*n+5)]
                else:
                    if i%4==0: neighbor = [count-1,count+1,count+2*n+4]
                    elif i%4==1: neighbor = [count-1,count+1,count-(2*n+3)]
                    elif i%4==2: neighbor = [count-1,count+1,count-(2*n+1)]
                    else: neighbor = [count-1,count+1,count+2*n+2]
                vertexdata[count]["neighbor"] = neighbor
                count+=1
        else:
            for i in range(2*n+4):
                if i==0: neighbor = [count-(2*n+4),count+1,count+2*n+3]
                elif i==2*n+3: neighbor = [count-1,count-(2*n+3),count+2*n+4]
                else:
                    if i%4==0: neighbor = [count-1,count+1,count-(2*n+4)]
                    elif i%4==1: neighbor = [count-1,count+1,count-(2*n+2)]
                    elif i%4==2: neighbor = [count-1,count+1,count+2*n+2]
                    else: neighbor = [count-1,count+1,count+2*n+4]
                vertexdata[count]["neighbor"] = neighbor
                count+=1

    if yPBC:
        for i in range(2*n+4):
            if i==0: neighbor = [count+1,count-(2*n+4),count+2*n+3]
            elif i==2*n+3: neighbor = [count-1,count-(2*n+3),count+n+2]
            else:
                if i%4==0: neighbor = [count-1,count+1,count-(2*n+4)]
                elif i%4==1: neighbor = [count-1,count+1,count-(2*n+2)]
                else: neighbor = [count-1,count+1,count+(2*n+2-int(i/4)*2)]
            vertexdata[count]["neighbor"] = neighbor
            count+=1 
        for i in range(n+2):
            if i==0: neighbor = [count+n+1,2*n+2,count-(2*n+2)]
            elif i==1: neighbor = [2*n+3,count+1,count-(2*n+2)]
            elif i==n+1: neighbor = [count-n-1,count-n-2,2*n-1]
            else:
                if i%2==0: neighbor = [count-1,i*2-2,count-(2*n+2-i)]
                else: neighbor = [count+1,2*i-3,count-(2*n+3-i)]
            vertexdata[count]["neighbor"] = neighbor
            count+=1
    else:
        for i in range(2*n+4):
            if i==0: neighbor = [count+1,count-(2*n+4),count+2*n+3]
            elif i==2*n+3: neighbor = [count-1,count-(2*n+3)]
            else: 
                if i%4==0: neighbor = [count-1,count+1,count-(2*n+4)]
                elif i%4==1: neighbor = [count-1,count+1,count-(2*n+2)]
                else: neighbor = [count-1,count+1]
            vertexdata[count]["neighbor"] = neighbor
            count+=1

    return vertexdata

def secondneighbor(vertexdata,n,m,a,xPBC,yPBC):
    Lplus1 = [0,math.sqrt(3)*a]
    Lplus2 = [3/2*a,-math.sqrt(3)*a/2]
    Lplus3 = [-3/2*a,-math.sqrt(3)*a/2]
    Lminus1 = [0,-math.sqrt(3)*a]
    Lminus2 = [3/2*a,math.sqrt(3)*a/2]
    Lminus3 = [-3/2*a,math.sqrt(3)*a/2]
    Lplus = [Lplus1,Lplus2,Lplus3]
    Lminus = [Lminus1,Lminus2,Lminus3]
    Rplus = [Lminus1,Lminus2,Lminus3]
    Rminus = [Lplus1,Lplus2,Lplus3]
    if xPBC: shiftx = 3*a*(n/2+1)
    else: shiftx = 0
    if yPBC: shifty = [3/2*a*(m/2+1),3/2*math.sqrt(3)*a*(m/2+1)]
    else: shifty=[0,0]

    for i in range(len(vertexdata)):
        if vertexdata[i]["parity"] ==0 :
            plus = Lplus
            minus = Lminus
        else: 
            plus = Rplus
            minus = Rminus
        nowpos = vertexdata[i]["pos"]
        vertexdata[i]["secondneighbor"] = []
        for j in range(len(vertexdata[i]["neighbor"])):
            neinum = vertexdata[i]["neighbor"][j]
            for k in range(len(vertexdata[neinum]["neighbor"])):
                nnnnum = vertexdata[neinum]["neighbor"][k]
                nnnpos0 = vertexdata[nnnnum]["pos"][0]
                nnnpos1 = vertexdata[nnnnum]["pos"][1]
                if nowpos[1]-nnnpos1 > 2*a:
                    nnnpos0 += shifty[0]
                    nnnpos1 += shifty[1]
                elif nnnpos1-nowpos[1] > 2*a:
                    nnnpos0 -= shifty[0]
                    nnnpos1 -= shifty[1]
            
                if nowpos[0]-nnnpos0 > 3/2*a*(m/2+1):  
                    nnnpos0+=shiftx
                elif nnnpos0-nowpos[0] > 3/2*a*(m/2+1): 
                    nnnpos0-=shiftx
                for p in range(3):
                    if abs((nnnpos0-nowpos[0]-plus[p][0])**2 + (nnnpos1-nowpos[1]-plus[p][1])**2) < 0.0001:
                        vertexdata[i]["secondneighbor"].append([nnnnum,1])
                        break
                for m in range(3):
                    if abs((nnnpos0-nowpos[0]-minus[m][0])**2 + (nnnpos1-nowpos[1]-minus[m][1])**2) < 0.0001:
                        vertexdata[i]["secondneighbor"].append([nnnnum,-1])
                        break
    return vertexdata

def lattice_plot(vertexdata,n,m,a,xPBC,yPBC,secondplot):
    if xPBC: shiftx = 3*a*(n/2+1)
    else: shiftx = 0
    if yPBC: shifty = [3/2*a*(m/2+1),3/2*math.sqrt(3)*a*(m/2+1)]
    else: shifty=[0,0]
    candidate = candidate_gen(xPBC,yPBC,n,m)
    fig = plt.figure(figsize=(n,m))
    ax = fig.add_subplot(1,1,1)
    filepath = '/Users/sogenikegami/Documents/UT4S/non-crystal/honeycomb/storage/'
    savename = str(n) + '_' + str(m) + 'lattice'
    for i in range(len(vertexdata)):
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
            else: 
                if (i in candidate) or (neinum in candidate):
                    plt.plot([nowpos[0],neipos0],[nowpos[1],neipos1],color="red")
                else:
                    plt.plot([nowpos[0],neipos0],[nowpos[1],neipos1],color="black")
    ax.set_aspect('equal')
    for i in range(len(vertexdata)):
        nowpos = vertexdata[i]["pos"]
        if vertexdata[i]["parity"] == 0:
            plt.scatter(nowpos[0],nowpos[1],c='red')
        else:
            plt.scatter(nowpos[0],nowpos[1],c='blue')

    if secondplot:
        for i in range(len(vertexdata)):
            nowpos = vertexdata[i]["pos"]
            for j in range(len(vertexdata[i]["secondneighbor"])):
                nnn = vertexdata[i]["secondneighbor"][j]
                nnnpos = vertexdata[nnn[0]]["pos"]
                if (i in candidate) or (nnn[0] in candidate):
                    plt.plot([nowpos[0],nnnpos[0]],[nowpos[1],nnnpos[1]],color="green")
                else:
                    plt.plot([nowpos[0],nnnpos[0]],[nowpos[1],nnnpos[1]],color="black")
    fig.show()
    plt.savefig(filepath + savename)
    return 0

def candidate_gen(xPBC,yPBC,n,m):
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

def defectnum_gen(vertexdata,defectratio,xPBC,yPBC,n,m,mode):
    candidate = []
    defectnum=[]
    if mode==0:
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
        n = int(len(candidate)*defectratio)
        randlist=[]
        while len(randlist)<n:
            temp = random.randint(0,len(candidate)-1)
            if temp in randlist:continue
            else: randlist.append(temp)

        for i in range(len(randlist)):
            defectnum.append(candidate[randlist[i]])
    
    elif mode==1:
        if xPBC:
            if yPBC: L = int((n/2+1)*(m/2+1))
            else: L=int((n/2+1)*m/2)
        else:
            if yPBC: L=int(n/2*(m/2+1))
            else:   L=int(n/2*m/2)
        n = int(L*defectratio)
        while len(defectnum)<n:
            temp = random.randint(0,len(vertexdata)-1)
            if temp in defectnum:continue
            else: defectnum.append(temp)

    return defectnum

def defect(vertexdata,defectnum,nnn):
    subtract_list=[]
    if 0 in defectnum: subtract_list.append(1)
    else: subtract_list.append(0)
    for i in range(len(vertexdata)-1):
        if i+1 in defectnum: subtract_list.append(subtract_list[i]+1)
        else: subtract_list.append(subtract_list[i])
    
    if nnn: 
        for i in range(len(vertexdata)):
            count1=0
            count2=0
            for j in range(len(vertexdata[i]["neighbor"])):
                neinum = vertexdata[i]["neighbor"][j-count1]
                if neinum in defectnum: 
                    del vertexdata[i]["neighbor"][j-count1]
                    count1+=1
                else: vertexdata[i]["neighbor"][j-count1] -= subtract_list[neinum]
            for j in range(len(vertexdata[i]["secondneighbor"])):
                nnnnum = vertexdata[i]["secondneighbor"][j-count2][0]
                if nnnnum in defectnum: 
                    del vertexdata[i]["secondneighbor"][j-count2]
                    count2+=1
                else: vertexdata[i]["secondneighbor"][j-count2][0] -= subtract_list[nnnnum]
    else:
        for i in range(len(vertexdata)):
            count1=0
            count2=0
            for j in range(len(vertexdata[i]["neighbor"])):
                neinum = vertexdata[i]["neighbor"][j]
                for k in range(len(vertexdata[neinum]["neighbor"])):
                    nnnnum = vertexdata[neinum]["neighbor"][k]
                    if nnnnum!=i and (neinum in defectnum):
                        for l in range(len(vertexdata[i]["secondneighbor"])):
                            if vertexdata[i]["secondneighbor"][l][0] == nnnnum:
                               del vertexdata[i]["secondneighbor"][l]
                               break
        for i in range(len(vertexdata)):
            count1=0
            count2=0
            for j in range(len(vertexdata[i]["neighbor"])):
                neinum = vertexdata[i]["neighbor"][j-count1]
                if neinum in defectnum: 
                    del vertexdata[i]["neighbor"][j-count1]
                    count1+=1
                else: vertexdata[i]["neighbor"][j-count1] -= subtract_list[neinum]
            for j in range(len(vertexdata[i]["secondneighbor"])):
                nnnnum = vertexdata[i]["secondneighbor"][j-count2][0]
                if nnnnum in defectnum: 
                    del vertexdata[i]["secondneighbor"][j-count2]
                    count2+=1
                else: vertexdata[i]["secondneighbor"][j-count2][0] -= subtract_list[nnnnum]                   

    
    defectset=[]
    for i in range(len(defectnum)):
        if defectnum[i] in defectset: continue
        else: defectset.append(defectnum[i])
    defectset.sort(reverse=True)
    for i in range(len(defectset)):
        del vertexdata[defectset[i]]

    return vertexdata


if __name__ == "__main__":
    main()