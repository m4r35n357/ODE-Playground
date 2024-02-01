import math
import numpy as np
import time
from matplotlib import pyplot as plt

def f(x):
    # problem 1 : Find the optimum value of f (minima function)
    #g = (x[0]**4-16*(x[0]**2)+5*x[0])/2+(x[1]**4-16*(x[1]**2)+5*x[1])/2

    #problem 2 :
    #d = 6
    #exp=math.exp(1)
    #p = np.zeros(d)
    #p[0] = x[0] + ((x[1]**2)*x[3]*x[5])*0.25 + 0.75
    #p[1] = x[1] + 0.405*exp**(1 + x[0]*x[1]) - 1.405
    #p[2] = x[2] - (x[3]*x[5])*0.5 + 1.5
    #p[3] = x[3] - 0.605*exp**(1 - x[2]**2) - 0.395
    #p[4] = x[4] - (x[1]*x[5])*0.5 + 1.5
    #p[5] = x[5] - x[0]*x[4]
    #pp=0
    #for i in range(d):
    #    pp=pp+p[i]**2
    ##g=-1/(1+pp)
    #problem 3 : sprinter finding problem
    #function to be minimized
    tau=[[12.27,11.57,11.54,12.07],[11.34,11.45,12.45,12.34],
         [11.29,11.50,11.45,11.52],[12.54,12.34,12.32,11.57],
         [12.20,11.22,12.07,12.03],[11.54,11.48,11.56,12.30]]
    F=0
    i=0
    for t in tau:
        for tt in t:
            F=F+tt*x[i]
            i=i+1

    #equality constraints
    e=np.zeros(4)
    e[0] = x[0]+x[4]+x[8]+x[12]+x[16]+x[20]-1
    e[1] = x[1]+x[5]+x[9]+x[13]+x[17]+x[21]-1
    e[2] = x[2]+x[6]+x[10]+x[14]+x[18]+x[22]-1
    e[3] = x[3]+x[7]+x[11]+x[15]+x[19]+x[23]-1
    ee=0
    for i in e:
        ee=ee+i**2

    #inequality constraints
    i=np.zeros(6)
    i[0] = x[0]+x[1]+x[2]+x[3]-1
    i[1] = x[4]+x[5]+x[6]+x[7]-1
    i[2] = x[8]+x[9]+x[10]+x[11]-1
    i[3] = x[12]+x[13]+x[14]+x[15]-1
    i[4] = x[16]+x[17]+x[18]+x[19]-1
    i[5] = x[20]+x[21]+x[22]+x[23]-1
    ii=0
    for j in i:
        ii=ii+np.max([j,0])**2

    #penalty function
    alpha = 10**15
    beta = 10**15
    g = F + alpha*ee + beta*ii
    return g

def val(x):
    d = 6
    exp=math.exp(1)
    p = np.zeros(d)
    p[0] = x[0] + ((x[1]**2)*x[3]*x[5])*0.25 + 0.75
    p[1] = x[1] + 0.405*exp**(1 + x[0]*x[1]) - 1.405
    p[2] = x[2] - (x[3]*x[5])*0.5 + 1.5
    p[3] = x[3] - 0.605*exp**(1 - x[2]**2) - 0.395
    p[4] = x[4] - (x[1]*x[5])*0.5 + 1.5
    p[5] = x[5] - x[0]*x[4]
    return p

def SOA_nd(xlow,xup,m,kmax,theta,r,d):
    ## Rotation matrix for n dimensional spiral
    R=np.identity(d)
    for i in range(d-1):
        for j in range(i):
            Rr=np.identity(d)
            Rr[d-2-i][d-2-i]=math.cos(theta)
            Rr[d-2-i][d-2+1-j]=-math.sin(theta)
            Rr[d-2+1-j][d-2-i]=math.sin(theta)
            Rr[d-2+1-j][d-2+1-j]=math.cos(theta)
            R=np.matmul(R,Rr)
    I=np.identity(d)
    D=r*I
    S=np.matmul(D,R)
    tm = time.time()
    x = np.random.uniform(xlow,xup,[m,d])                    #Generate random numbers
    #x = xlow+(xup-xlow)*np.random.rand(m,d)

    #spiral algorithm
    xstar=x[0]
    for j in range(kmax):                                   #finding the next center of spiral
        for i in range (m):
            if f(x[i])<f(xstar):
                found=True
                for k in x[i]:
                    if xlow < k < xup:
                        found=True
                    else:
                        found=False
                        break
                if found:
                    xstar=x[i]

    #    spiral process
        v=np.zeros(shape=[m,d])
        for i in range(m):
            v[i]=np.matmul(S,x[i])-np.matmul((S-I),xstar)
        x=v
        plt.plot([X[2] for X in x],[X[5] for X in x],'ro')
        plt.axis([xlow,xup,xlow,xup])
        plt.show()

    elapsed = time.time() - tm
    print('execution time:',elapsed,'seconds')
    print('xstar:', xstar)
    print('feval:',f(xstar))
    print('value:', val(xstar))

#SOA_nd(0,1,100,100,math.pi/4,0.95,24)
SOA_nd(0,1,100,1,math.pi/4,0.95,24)
