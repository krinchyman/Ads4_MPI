#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import numpy as np

from matplotlib import pyplot as plt

import logging
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)

def derivada1a(dx, vector):

    # Simétrica en 5 puntos:
    der = np.roll(vector,2) - 8*np.roll(vector,1) + 8*np.roll(vector,-1) - np.roll(vector,-2)
    
    # La formulación anterior no es válida en i=[0,1,Nx-1,Nx]:
    der[0] = -25*vector[0] + 48*vector[1] - 36*vector[2] + 16*vector[3] -3*vector[4]
    der[1] = -10*vector[1] - 3*vector[0] + 18*vector[2] - 6*vector[3] + vector[4]
    
    der[-1] = 25*vector[-1] - 48*vector[-2] + 36*vector[-3] - 16*vector[-4] + 3*vector[-5]
    der[-2] = 10*vector[-2] + 3*vector[-1] - 18*vector[-3] + 6*vector[-4] - vector[-5]

    der /= 12*dx
    
    return der
    
def derivada2a(dx, vector):
    
    ''' Función para calcular la derivada segunda. Nótese que la prescripción en [0,1] es la misma que en [Nx-1,Nx],
        ya que cambiar h por -h no tiene efecto como en la derivada primera. '''

    # Simétrica en 5 puntos:
    der = -np.roll(vector,2) + 16*np.roll(vector,1) - 30*vector + 16*np.roll(vector,-1) - np.roll(vector,-2)

    # La formulación anterior no es válida en i=[0,1,Nx-1,Nx]:
    der[0] =  35*vector[0] - 104*vector[1] + 114*vector[2] - 56*vector[3] + 11*vector[4]
    der[1] = -20*vector[1] +  11*vector[0] +   6*vector[2] +  4*vector[3] -    vector[4]
    
    #der[-1] =  35*vector[-1] - 104*vector[-2] + 114*vector[-3] -  56*vector[-4] + 11*vector[-5]
    der[-1] =  45*vector[-1] - 154*vector[-2] + 214*vector[-3] - 156*vector[-4] + 61*vector[-5] - 10*vector[-6]
    der[-2] = -20*vector[-2] +  11*vector[-1] +   6*vector[-3] +  4*vector[-4] -    vector[-5]

    der /= 12*dx**2

    return der

def derivada(dx, vector):

    # Simétrica en 5 puntos:
    der = np.roll(vector,2) - 8*np.roll(vector,1) + 8*np.roll(vector,-1) - np.roll(vector,-2)
    
    # La formulación anterior no es válida en i=[0,1,Nx-1,Nx].
    # Ojo aquí!!! Lo que sigue solo es válido para el caso que nos ocupa:

    der[0] = 0 # Este término es idénticamente cero por ser función par.
    der[1] = vector[1] - 8*vector[0] + 8*vector[2] -vector[3] # Nóta: El primer dato es en relidad vector[-1]=vector[1] si f es par en el origen

    der /= 12*dx
    
    return der

def euler(dx, vector, y0):

    integral = dx*vector
    integral[0] = y0 # Esta es en realidad la constante de integración

    return integral.cumsum()

def trapecio(dx, vector, y0):

    integral = np.roll(vector,1) + vector

    integral *= dx/2

    integral[0] = y0

    return integral.cumsum()

def rk4_novec(dx, vector, y0):

    integral    = np.zeros(vector.shape)
    integral[0] = y0

    for i in range(1,len(vector)):

        # RK4 es lo mismo que Simpson si no hay dependencia en y:
        integral[i] = integral[i-1] + dx*(vector[i-1] + 4*(vector[i-1]+vector[i])/2 + vector[i])/6

    return integral

def rk4(dx, vector, y0):

    # No existe dependencia explícita en y, por lo tanto RK4 es la regla de Simpson:
    integral  = np.roll(vector,1) + 4*(np.roll(vector,1)+vector)/2 + vector
    integral *= dx/6

    integral[0] = y0 # Esta es en realidad la constante de integración
   
    return integral.cumsum()

def alex(dx, vector, y0):

    Nx, = vector.shape
    Nx -= 1

    integral     = -np.roll(vector,2) + 13*np.roll(vector,1) + 13*vector - np.roll(vector,-1)
    integral[1]  = 9*vector[0] + 19*vector[1] - 5*vector[2] + vector[3]
    integral[Nx] = vector[Nx-3] - 5*vector[Nx-2] + 19*vector[Nx-1] + 9*vector[Nx]

    integral   *= dx/24

    integral[0] = y0 # Esta es en realidad la constante de integración

    return integral.cumsum()


def sim5pto(dx, vector, y0):

    Nx, = vector.shape
    Nx -= 1

    integral       = -19*np.roll(vector,2) + 346*np.roll(vector,1) + 456*vector - 74*np.roll(vector,-1) + 11*np.roll(vector,-2)
    integral[Nx-1] =  11*vector[Nx-4] -  74*vector[Nx-3] + 456*vector[Nx-2] + 346*vector[Nx-1] -  19*vector[Nx]
    integral[Nx]   = -19*vector[Nx-4] + 106*vector[Nx-3] - 264*vector[Nx-2] + 646*vector[Nx-1] + 251*vector[Nx]

    integral[1]    = -19*vector[4] + 106*vector[3] - 264*vector[2] + 646*vector[1] + 251*vector[0]

    integral   *= dx/720

    integral[0] = y0 # Esta es en realidad la constante de integración

    return integral.cumsum()


if __name__=='__main__':


    #####################################
    # Sección de test de las subrutinas #
    #####################################

    x = np.linspace(0, np.pi/2, 1000, dtype=np.longdouble)
        
    dx = x[1] - x[0]

    integrando  = np.sin(x)
    analitica  = 1 - np.cos(x)

    plt.plot(x,analitica,'b-')

    solucion = euler(dx, integrando, 0)
    plt.plot(x,solucion,'g-')

    solucion = trapecio(dx, integrando, 0)
    plt.plot(x,solucion,'c-')

    solucion = rk4(dx, integrando, 0)
    plt.plot(x,solucion,'r-')

    solucion = alex(dx, integrando, 0)
    plt.plot(x,solucion,'k-')

    plt.grid(True)

    plt.show()

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # OJO!!!! IMPORTANTE!!!! No! Mejor: IMPORTANTISISIMO !!!!
    # Con el algoritmo de resolción de Alex, los reales  !!!!
    # se desbordan! Comparar las dos versiones!!         !!!!
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    intervalo = np.pi/2

    N1 = 1000;
    N2 = 2000;
    N4 = 4000;

    grid1 = np.longdouble(intervalo/N1)*np.arange(N1+1)
    grid2 = np.longdouble(intervalo/N2)*np.arange(N2+1)
    grid4 = np.longdouble(intervalo/N4)*np.arange(N4+1)

    grid1 = (intervalo/N1)*np.arange(N1+1)
    grid2 = (intervalo/N2)*np.arange(N2+1)
    grid4 = (intervalo/N4)*np.arange(N4+1)

    integrando1 = np.sin(grid1)
    integrando2 = np.sin(grid2)
    integrando4 = np.sin(grid4)

    analitica1  = 1 - np.cos(grid1)
    analitica2  = 1 - np.cos(grid2)
    analitica4  = 1 - np.cos(grid4)

    res1 = alex(intervalo/N1,integrando1,0)
    res2 = alex(intervalo/N2,integrando2,0)
    res4 = alex(intervalo/N4,integrando4,0)

    plt.plot(grid1,analitica1-res1)
    plt.plot(grid2,analitica2-res2)
    plt.plot(grid4,analitica4-res4)
    plt.grid(True)
    plt.show()

    # A continuación hacemos lo mismo quitando el valor en el origen:
    analitica  = 1 - np.cos(grid1)
    analitica  = analitica[1:]

    res1 = alex(intervalo/N1,integrando1,0)[::1][1:]
    res2 = alex(intervalo/N2,integrando2,0)[::2][1:]
    res4 = alex(intervalo/N4,integrando4,0)[::4][1:]

    plt.plot(analitica-res1,'b')
    plt.plot(analitica-res2,'g')
    plt.plot(analitica-res4,'r')
    plt.grid(True)
    plt.show()

    plt.plot((res1-res2)/(res2-res4),'b.')
#    plt.axis([0,1000,14.5,17.5])
    plt.grid(True)
    plt.show()

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # OJO!!!! IMPORTANTE!!!! No! Mejor: IMPORTANTISISIMO !!!!
    # Con el algoritmo de resolción de Alex, los reales  !!!!
    # se desbordan! Comparar las dos versiones!!         !!!!
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    resoluciones = {1000: 4, 2000 : 2, 4000 : 1}
    x = np.linspace(0, np.pi/2, 4000, dtype=np.float64)
    x = np.linspace(0, np.pi/2, 4000, dtype=np.longdouble) # np.logdouble==np.float96 e mi máquina que no es 64bits. Muy costosa la aritmética!!!

    soluciones = {} # Guardamos las soluciones para las distintas resoluciones
    X          = {} # Guardamos las x solo para comprobar que lo estamos haciendo bien

    for resolucion,idx in resoluciones.iteritems():

        tmp      = x[::idx]
        dx       = tmp[1] - tmp[0]

        integrando  = np.sin(tmp)

        soluciones[resolucion] = sim5pto(dx, integrando, 0)
        X[resolucion]          = tmp

    # Quitamos el primer punto para evitar la divisón por 0:
    s1 = soluciones[1000][::1][1::]
    s2 = soluciones[2000][::2][1::]
    s3 = soluciones[4000][::4][1::]

    x1 = X[1000][::1][1::]
    x2 = X[2000][::2][1::]
    x3 = X[4000][::4][1::]

    analitica  = 1 - np.cos(x1)

    # Ploteado:
    estilos = {1000: 'b-', 2000 : 'g-', 4000 : 'r-'}
    plt.figure()

    plt.plot(s1-analitica,estilos[1000], label='1000')
    plt.plot(s2-analitica,estilos[2000], label='2000')
    plt.plot(s3-analitica,estilos[4000], label='4000')

    plt.grid(True)

    plt.figure()
    plt.plot((s1-s2)/(s2-s3))

    plt.axis([0,1000,-1,20])

    plt.grid(True)

    plt.show()


    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # OJO!!!! IMPORTANTE!!!! No! Mejor: IMPORTANTISISIMO !!!!
    # Con el algoritmo de resolción de Alex, los reales  !!!!
    # se desbordan! Comparar las dos versiones!!         !!!!
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    resoluciones = {100: 4, 200 : 2, 400 : 1}
    x = np.linspace(0, np.pi/2, 4000, dtype=np.float64)
    x = np.linspace(0, np.pi/2, 4000, dtype=np.longdouble) # np.logdouble==np.float96 e mi máquina que no es 64bits. Muy costosa la aritmética!!!

    soluciones = {} # Guardamos las soluciones para las distintas resoluciones
    X          = {} # Guardamos las x solo para comprobar que lo estamos haciendo bien

    for resolucion,idx in resoluciones.iteritems():

        tmp      = x[::idx]
        dx       = tmp[1] - tmp[0]

        derivando  = np.sin(tmp)

        soluciones[resolucion] = derivada1a(dx, derivando)
        X[resolucion]          = tmp

    # Quitamos el primer punto para evitar la divisón por 0:
    s1 = soluciones[100][::1][1::]
    s2 = soluciones[200][::2][1::]
    s3 = soluciones[400][::4][1::]

    x1 = X[100][::1][1::]
    x2 = X[200][::2][1::]
    x3 = X[400][::4][1::]

    analitica  = np.cos(x1)

    # Ploteado:
    estilos = {100: 'b-', 200 : 'g-', 400 : 'r-'}
    plt.figure()

    plt.plot(s1-analitica,estilos[100], label='100')
    plt.plot(s2-analitica,estilos[200], label='200')
    plt.plot(s3-analitica,estilos[400], label='400')

    plt.legend()

    plt.grid(True)

    plt.figure()
    plt.plot((s1-s2)/(s2-s3))

    plt.grid(True)

    plt.show()


#   New:

    resoluciones = 2**np.arange(3,14)

    integral = []
    h        = []

    for resolucion in resoluciones:

        x  = np.linspace(0, np.pi/2, resolucion, dtype=np.float96)
        dx = x[1] - x[0]

        h.append( dx )

        integrando = np.cos(x)
        val1 = sim5pto(dx, integrando, 0)[-1]
        integral.append( val1 )
    
    h        = np.array(h)
    integral = np.array(integral)

    plt.loglog(h,h,'k--')
    plt.loglog(h,h**2,'b--')
    plt.loglog(h,h**3,'g--')
    plt.loglog(h,h**4,'r--')
    plt.loglog(h,h**5,'k--')
    plt.loglog(h,np.abs(integral-1))

    plt.grid()

    plt.show()

    nd = len(integral)

    # Conociendo la solución analítica:
    for i in range(nd-1):
    
        print np.log2((integral[i]-1)/(integral[i+1]-1))


    resoluciones = 2**np.arange(3,14)

    derivada = []
    h        = []

    for resolucion in resoluciones:

        x  = np.linspace(0, np.pi/2, resolucion+1, dtype=np.float64)
        dx = x[1] - x[0]

        h.append( dx )

        derivando = np.sin(x)
        derivada.append( derivada1a(dx, derivando)[-1] )
    
    h        = np.array(h)
    derivada = np.array(derivada)

    plt.loglog(h,h,'k--')
    plt.loglog(h,h**2,'b--')
    plt.loglog(h,h**3,'g--')
    plt.loglog(h,h**4,'r--')
    plt.loglog(h,h**5,'k--')
    plt.loglog(h,np.abs(derivada))

    plt.grid()

    plt.show()

    nd = len(derivada)

    # Si no conociesemos la solución analítica:
    for i in range(nd-2):
    
        print np.log2(np.abs((derivada[i]-derivada[i+1])/(derivada[i+1]-derivada[i+2])))

    # Conociendo la solución analítica:
    for i in range(nd-1):
    
        print np.log2(np.abs((derivada[i])/(derivada[i+1])))
