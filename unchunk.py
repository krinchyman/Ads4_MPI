#!/usr/bin/env python
# -*- coding: utf-8 -*-

from netCDF4 import Dataset
from matplotlib import pyplot as plt

import numpy as np

import argparse

def crea_evolucion(fout, x):

    destino = Dataset(filename=fout, mode='w', clobber=True, format='NETCDF3_CLASSIC')

    destino.createDimension('time', None)
    destino.createDimension('x', len(x))

    variable    = destino.createVariable('x','f4',('x',))
    variable[:] = x[:]
    variable.standard_name = "x"
    variable.long_name = "Coordenada holografica"

    variable         = destino.createVariable('time','f4',('time',))
    variable.standard_name = "time"
    variable.long_name = "time"

    variables = ['A', 'delta', 'Pi', 'Phi', 'MomentumConstraint']
    #variables = ['A', 'delta', 'Pi', 'Phi']

    for nombre in variables:

        variable    = destino.createVariable(nombre,'f4',('time','x'))

    destino.close()

def buffer_out(fout, t, Phi, Pi, A, delta, constraint, **kwargs):
#def buffer_out(fout, t, Phi, Pi, A, delta, **kwargs):

#    variables = {'A' : A, 'delta' : delta, 'Pi' : Pi, 'Phi' : Phi}
    variables = {'A' : A, 'delta' : delta, 'Pi' : Pi, 'Phi' : Phi, 'MomentumConstraint' : constraint}

    destino = Dataset(filename=fout, mode='a')

    Nt = len(destino.dimensions['time'])

    destino.variables['time'][Nt] = t

    for n, v in variables.iteritems():

        variable       = destino.variables[n]
        variable[Nt,:] = v

    # Atribulos del run:
    if kwargs:

        destino.setncatts(kwargs)

    destino.close()

    return

if __name__=='__main__':

    # Procesado de la linea de comandos:
    parser = argparse.ArgumentParser()
    parser.add_argument('-n','--procs',help='Número de procesadores',type=int)
    parser.add_argument('-a','--alpha',help='Parámetro alpha',type=float)
    parser.add_argument('-t','--tau',help='Parámetro tau',type=float)
    parser.add_argument('-p','--puntos',help='Número de puntos',type=int)
    parser.add_argument('-d','--derivadas',help='Estilo derivadas',type=str)

    args   = parser.parse_args()

    # Configuración:
    nChunks   = args.procs
    alpha     = args.alpha
    tau       = args.tau
    points    = args.puntos
    derivadas = args.derivadas

    # Ficheros de entrada:
    ficheros = ['%s_evolucion_Source_a%4.2f_t%4.2fx%i.%i.nc' % (derivadas, alpha, tau, points, i) for i in range(nChunks)]

    X = []
    for fichero in ficheros:
        datos = Dataset(fichero)
        X.append( datos.variables['x'][:] )
        Nt = len(datos.dimensions['time'])
        datos.close()

    x = np.concatenate(X)

    # Fichero de salida:
    fout = '%s_evolucion_Source_a%4.2f_t%4.2fx%i.nc' % (derivadas, alpha, tau, points)
    crea_evolucion(fout, x)

    for t in range(Nt):

        variables = {'A' : [], 'delta' : [], 'Pi' : [], 'Phi' : [], 'MomentumConstraint' : []}
    
        for nombre, lista in variables.iteritems():

            for fichero in ficheros:

                datos = Dataset(fichero)
                time  = datos.variables['time'][t]
                lista.append( datos.variables[nombre][t,:] )
                datos.close()                


        buffer_out(fout, time, np.concatenate(variables['Phi']),
                            np.concatenate(variables['Pi']),
                            np.concatenate(variables['A']),
                            np.concatenate(variables['delta']),
                            np.concatenate(variables['MomentumConstraint']))

