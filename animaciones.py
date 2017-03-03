#! /usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import platform

import argparse

import numpy as np

import matplotlib

if platform.system() == 'Darwin':
    matplotlib.use('TkAgg')

from matplotlib import pyplot as plt
import matplotlib.animation as animation

from derivacion_integracion import *

from netCDF4 import Dataset

import scipy.optimize as sc

import logging
logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)

epsilon = 0.5 
tau = 5.

def Pi0(t):
    
    if t==0:
        
        return 0
        
    else:    
  
#       No analítico:
        return epsilon*np.exp(-(tau/t)**2)

def dPi0(t):
    
    if t==0:
        
        return 0
        
    else:    
  
#       No analítico:
        return 2*epsilon*np.exp(-(tau/t)**2)*tau**2/t**3

def ddPi0(t):
    
    if t==0:
        
        return 0
        
    else:    
  
#       No analítico:
        return 2*epsilon*tau**2*(2*tau**2-3*t**2)*np.exp(-(tau/t)**2)/t**6

if __name__=='__main__':

    # Parseado de la linea de comandos:
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--fichero', help='Fichero con los datos a procesar', type=str)
    parser.add_argument('-i','--iter', help='Escoger instantes de tiempo', type=int)
    parser.add_argument('-v','--video', help='Extraer animaciones', action='store_true')
    parser.add_argument('-t','--tiempo', help='Tiempo final', type=float)
    args = parser.parse_args()

    fichero = args.fichero
    Video   = args.video
    n       = args.iter
    fin     = args.tiempo

    a = [0,1.7,-1.5,1.5]

    datos = Dataset(fichero)

    x     = datos.variables['x'][:]
    dx    = x[1] - x[0]
    Nx,   = x.shape
    Nx   -= 1

    t       = datos.variables['time'][:]
    dt      = t[11] - t[10]

    if fin>t[-1]:
        fin = t[-1]


    Frames  = int(fin/dt)

    t     = datos.variables['time'][0:Frames:n]
    Pi    = datos.variables['Pi'][0:Frames:n,:]
    Phi   = datos.variables['Phi'][0:Frames:n,:]
    A     = datos.variables['A'][0:Frames:n,:]
    delta = datos.variables['delta'][0:Frames:n,:]

    Frames, = t.shape

    intervalo =50 

    # Comprobamos si hay video encoder disponible:
    if Video == True: 

        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=intervalo, metadata=dict(artist='Pablo E. Carracedo'), bitrate=1800)

    # Phi:
    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(x.min(), x.max()), ylim=(Phi.min(), Phi.max()))
    ax.set_xlabel('x')
    ax.set_ylabel('$\Phi$')
    ax.grid()

    ims = []

    for Frame in range(Frames):

        line, = ax.plot(x,Phi[Frame], 'b')
        texto = ax.text(0.75, 0.9, 'Time: %5.3f' % t[Frame], transform=ax.transAxes)

        ims.append([line, texto])

    if Video==True:

        ani = animation.ArtistAnimation(fig, ims, interval=intervalo, blit=True)
        ani.save('Phi.mp4', writer=writer)

    else:

        ani = animation.ArtistAnimation(fig, ims, interval=intervalo, blit=True, repeat_delay=1000)
        plt.show()

    # VeV:
    ddPhi = np.empty_like(Phi)

    for Frame in range(Frames):

        ddPhi[Frame,:] = derivada2a(dx,Phi[Frame,:])

    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(x.min()-0.1, x.max()+0.1), ylim=(ddPhi.min(), ddPhi.max()+0.005))
    #ax = fig.add_subplot(111, autoscale_on=False, xlim=(x.min(), x.max()), ylim=(-10,10))
    #ax = fig.add_subplot(111, autoscale_on=False, xlim=(1.50, 1.58), ylim=(-0.005,0.005))
    ax.set_xlabel('x')
    ax.set_ylabel('$d^2\Phi/dx^2$')
    ax.grid()

    ims = []

    for Frame in range(Frames):

        line, = ax.plot(x,ddPhi[Frame], 'b')
        texto = ax.text(0.75, 0.9, 'Time: %5.3f' % t[Frame], transform=ax.transAxes)

        ims.append([line, texto])


    if Video==True:

        ani = animation.ArtistAnimation(fig, ims, interval=intervalo, blit=True)
        ani.save('dPhi.mp4', writer=writer)

    else:

        ani = animation.ArtistAnimation(fig, ims, interval=intervalo, blit=True, repeat_delay=1000)
        plt.show()

#   Pi:
    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(x.min(), x.max()), ylim=(Pi.min(), Pi.max()))
    ax.set_xlabel('x')
    ax.set_ylabel('$\Pi$')
    ax.grid()

    ims = []

    for Frame in range(Frames):

        line, = ax.plot(x,Pi[Frame], 'b')
        texto = ax.text(0.75, 0.9, 'Time: %5.3f' % t[Frame], transform=ax.transAxes)

        ims.append([line, texto])


    if Video==True:

        ani = animation.ArtistAnimation(fig, ims, interval=intervalo, blit=True)
        ani.save('Pi.mp4', writer=writer)

    else:

        ani = animation.ArtistAnimation(fig, ims, interval=intervalo, blit=True, repeat_delay=1000)
        plt.show()

#   Integrando:
    integrando = np.empty_like(Phi)
    integral   = np.empty_like(Phi)


    for Frame in range(Frames):

        rho                  = Phi[Frame,:]**2 + Pi[Frame,:]**2
        integrando[Frame,:]  = rho*np.exp(-delta[Frame,:])  - np.exp(-delta[Frame,Nx])*Pi[Frame,Nx]**2
        integrando[Frame,:] *= np.tan(x)**2
        integrando[Frame,-1] = 5*Pi0(t[Frame])**4/2 + dPi0(t[Frame])**2 - Pi0(t[Frame])*ddPi0(t[Frame])
        integrando[Frame,-1] = 5*Pi[Frame,Nx]**4/2 # Para las simulaciones estáticas
        integral[Frame,:]    = alex(dx, integrando[Frame,:], 0)

    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(x.min(), x.max()), ylim=(integrando.min(), integrando.max()))
    ax.set_xlabel('x')
    ax.set_ylabel('Integrando')
    ax.grid()

    ims = []

    for Frame in range(Frames):

        line,  = ax.plot(x,integrando[Frame], 'b')
        estrella = 5*Pi0(t[Frame])**4/2 + dPi0(t[Frame])**2 - Pi0(t[Frame])*ddPi0(t[Frame])
        line2, = ax.plot(x[-1], estrella, 'r*')
        line3,  = ax.plot(x,integral[Frame], 'y')
        texto = ax.text(0.75, 0.9, 'Time: %5.3f' % t[Frame], transform=ax.transAxes)

        ims.append([line, line2, line3, texto])


    if Video==True:

        ani = animation.ArtistAnimation(fig, ims, interval=intervalo, blit=True)
        ani.save('Integral.mp4', writer=writer)

    else:

        ani = animation.ArtistAnimation(fig, ims, interval=intervalo, blit=True, repeat_delay=1000)
        plt.show()


#   delta:
    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(x.min(), x.max()), ylim=(delta.min(), delta.max()))
    ax.set_xlabel('x')
    ax.set_ylabel('$\delta$')
    ax.grid()

    ims = []

    for Frame in range(Frames):

        line, = ax.plot(x,delta[Frame], 'b')
        texto = ax.text(0.75, 0.9, 'Time: %5.3f' % t[Frame], transform=ax.transAxes)

        ims.append([line, texto])


    if Video==True:

        ani = animation.ArtistAnimation(fig, ims, interval=intervalo, blit=True)
        ani.save('delta.mp4', writer=writer)

    else:

        ani = animation.ArtistAnimation(fig, ims, interval=intervalo, blit=True, repeat_delay=1000)
        plt.show()


#   alpha:
    alpha = np.empty_like(Phi)

    for Frame in range(Frames):

        alpha[Frame,:] = Pi[Frame]/(np.exp(delta[Frame])/A[Frame])

    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(x.min(), x.max()), ylim=(alpha.min(), alpha.max()))
    ax.set_xlabel('x')
    ax.set_ylabel(r'$\alpha=\Pi A /e^{\delta}$')
    ax.grid()

    ims = []

    for Frame in range(Frames):

        line, = ax.plot(x,alpha[Frame], 'b')
        texto = ax.text(0.75, 0.9, 'Time: %5.3f' % t[Frame], transform=ax.transAxes)

        ims.append([line, texto])


    if Video==True:

        ani = animation.ArtistAnimation(fig, ims, interval=intervalo, blit=True)
        ani.save('aPi0.mp4', writer=writer)

    else:

        ani = animation.ArtistAnimation(fig, ims, interval=intervalo, blit=True, repeat_delay=1000)
        plt.show()


#   A:
    fig = plt.figure()
    ax = fig.add_subplot(111, autoscale_on=False, xlim=(x.min(), x.max()), ylim=(A.min(), A.max()))
    ax.set_xlabel('x')
    ax.set_ylabel('$A$')
    ax.grid()

    ims = []

    Na = len(datos.ncattrs())

    #if Na!=0:
    if False:

        M0 = datos.Masa
        R  = datos.Radio

        def f(x):
            return 1 - M0*np.cos(x)**3/np.sin(x)

        A_Sch = f(x)
        R_Sch = sc.newton(f,0.01)

    for Frame in range(Frames):

        #if Na!=0:
        if False:
            line, = ax.plot(x, A[Frame], 'b', x, A_Sch, 'k--')
        else:
            line, = ax.plot(x, A[Frame], 'b')

        texto = ax.text(0.75, 0.9, 'Time: %5.3f' % t[Frame], transform=ax.transAxes)

        ims.append([line, texto])

        a = plt.axis()

        #if Na!=0:
        if False:
            plt.axis([a[0],0.2,0,1.2])
            plt.xticks([R, R_Sch, 0.1, 0.2], ['$R$', '$R_{Sch}$', 0.1, 0.2])

    for i in range(2*intervalo):
        ims.append(ims[-1])


    if Video==True:

        ani = animation.ArtistAnimation(fig, ims, interval=intervalo, blit=True)
        ani.save('A.mp4', writer=writer)

    else:

        ani = animation.ArtistAnimation(fig, ims, interval=intervalo, blit=True, repeat_delay=1000)
        plt.show()
