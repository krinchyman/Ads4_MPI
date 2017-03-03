# Ads4_MPI

Este proyecto tiene los códigos MPI para la resolución de AdS4:

- RK4_Source_estatico_MPI_DC/DD.py: Son las versiones viejas para el caso estacionario.
- RK4_Source_MPI_DC/DD.py: Son las nuevas versiones para source dependiente del tiempo.
- derivacion_integracion.py: Es una libreria con funciones de derivación/integración utilizadad por el programa animaciones.py
- animaciones.py: Programa para generar animaciones a partir de los ficheros de evolución en formato NetCDF
- unchunk.py: Programa para unir varios chunks NetCDF de la ejecución en paralelo en un solo fichero.
