from math import *
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as si


#---------------------------------------------------------------
#          definition de fonction
#---------------------------------------------------------------



def teta(t):
    return (pi/2)*cos(3.132*t)

def Pendule(Y):
    g = 9.81
    L = 1
    w = sqrt(g/L)

    return np.array([Y[1],-(w**2)*sin(Y[0])])

def Euler(f,y0,h):
    N =100
    EulerE = np.zeros((N,2))
    EulerE[0,:] = y0
    
    for k in range(N-1):
        EulerE[k+1,:] = EulerE[k,:]+ h*f(EulerE[k,:])

    return EulerE


def RK(f,y0,h):    
    N =100
    Ykr = np.zeros((N,2))
    Ykr[0,:] = y0
    
    for k in range(N-1):
        k1 = f(Ykr[k,:])
        k2 = f(Ykr[k,:] + h/2*k1)
        k3 = f(Ykr[k,:] + h/2*k2)
        k4 = f(Ykr[k,:] + h*k3)
        
        Ykr[k+1,:] = Ykr[k,:] + h/6*(k1 + 2*k2 + 2*k3 + k4)
        
    return Ykr

def RK2(f,y0,h):    
    N =100
    Ykr2 = np.zeros((N,2))
    Ykr2[0,:] = y0
    
    for k in range(N-1):
        ynn  = Ykr2[k,:] + (h/2)*f(Ykr2[k,:])
        Ykr2[k+1,:] = Ykr2[k,:] + h*f(ynn)       
    return Ykr2

   
#---------------------------------------------------------------
#              programme principal
#---------------------------------------------------------------

#4.1 q4
t= []
Lteta = []
t = np.linspace(0,4,100)
for i in t: 
    Lteta.append(teta(i))

plt.figure(1)
plt.title('theta(t) = Acos(wt+phi)')
plt.plot(t,Lteta,label="theta(t)")

plt.xlabel("temps (s)")
plt.ylabel("theta")
plt.grid()
plt.legend()
plt.show()




#4.2.2

y0 = np.array([pi/2 , 0])
h = 0.04
E = Euler(Pendule,y0,h)


plt.figure(2)
plt.title('theta(t) = Acos(wt+phi)')
plt.plot(t,Lteta,label="theta(t)")
plt.plot(t,E[:,0],label = 'Euler Explicite')
plt.xlabel("temps (s)")
plt.ylabel("theta")
plt.grid()
plt.legend()
plt.show()

#4.2.3


Yrk=RK(Pendule,y0,h)

plt.figure(3)
plt.title('theta(t) = Acos(wt+phi)')
plt.plot(t,Lteta,label="theta(t)")
plt.plot(t,E[:,0],label = 'Euler Explicite')
plt.plot(t,Yrk[:,0], label="Runge Kutta")
plt.xlabel("temps (s)")
plt.ylabel("theta")
plt.grid()
plt.legend()
plt.show()

#4.2.4

def pendule2(Y,t):
    g=9.81
    L=1
    yprime=np.zeros(2)
    yprime[0]=Y[1]
    yprime[1]=-(g/L)*sin(Y[0])
    return yprime

Yode=si.odeint(pendule2,y0,t)

plt.figure(4)
plt.title('theta(t) = Acos(wt+phi)')
plt.plot(t,Lteta,label="theta(t)")
plt.plot(t,E[:,0],label = 'Euler Explicite')
plt.plot(t,Yrk[:,0], label="Runge Kutta")
plt.plot(t,Yode[:,0],label="odeint",color="red")
plt.xlabel("temps (s)")
plt.ylabel("theta")
plt.grid()
plt.legend(loc = 3)
plt.show()


#---------------------------------------------------------------
#                       Partie 5
#---------------------------------------------------------------

#1
theta1, theta_p = si.odeint(pendule2,y0,t).T
theta2, theta_pp = Euler(Pendule,y0,h).T
theta3, theta_ppp = RK(Pendule,y0,h).T

plt.figure(5)
plt.plot(theta1, theta_p,label="odeint",color='red')
plt.plot(theta3, theta_ppp,label="Runge-Kutta")
plt.plot(theta2, theta_pp,label="Euler Explicite")
plt.xlabel('Theta')
plt.ylabel('Theta_point')
plt.legend(loc=3)
plt.grid()
plt.title('portrait de phase')
plt.show()

#2
Yrk2=RK2(Pendule,y0,h)

plt.figure(6)
plt.title('theta(t) = Acos(wt+phi)')
plt.plot(t,Lteta,label="theta(t)")
plt.plot(t,Yrk[:,0], label="Runge Kutta4")
plt.plot(t,Yrk2[:,0], label="Runge Kutta2")
plt.xlabel("temps (s)")
plt.ylabel("theta")
plt.grid()
plt.legend()
plt.show()

#---------------------------------------------------------------
#                       Partie 6
#---------------------------------------------------------------


def suspension(Y,t):
    yp = np.zeros(4)
    M1 = 15
    M2 = 200
    C2 = 1200
    K1 = 50000
    K2 = 5000

    yp[0] = Y[1]
    yp[1] = (1/M1)*(C2*(Y[3]-Y[1]) + (-K1 + K2)*Y[0] + K2*Y[2])
    yp[2] = Y[3]
    yp[3] = (1/M2)*(C2*(Y[1]-Y[3]) + K2*(Y[0]-Y[2]) -1000)
    return yp

y00 = np.array([0 , 0 , 0 , 0])
tfin = 3

t2 = np.linspace(0,3,100)
Yode2 = si.odeint(suspension,y00,t2)

plt.figure(7)
plt.title("x1(t) et x2(t)")
plt.plot(t2,Yode2[:,0],label="x1(t)")
plt.plot(t2,Yode2[:,2],label="x2(t)")
plt.xlabel("Temps (s)")
plt.ylabel("x(t)")
plt.legend()
plt.show()  




