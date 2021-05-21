# -*- coding: utf-8 -*-
import numpy as np
from scipy.fftpack import hilbert
from obspy.core import read    
import matplotlib.pyplot as plt
from pylab import *


st = read('C:/Users/Propietario/Desktop/tesis_de_grado/ATENUACION/datos_prueba/Trazas_Fuego/GI.FG16.00.BHZ.D.2020.218')

tr = st[0]
data = tr.data

def envelope(data):
    """
    Envelope of a function.
    Computes the envelope of the given function. The envelope is determined by
    adding the squared amplitudes of the function and it's Hilbert-Transform
    and then taking the square-root. (See [Kanasewich1981]_)
    The envelope at the start/end should not be taken too seriously.
    :type data: numpy.ndarray
    :param data: Data to make envelope of.
    :return: Envelope of input data.
    """
    hilb = hilbert(data)
    data = (data ** 2 + hilb ** 2) ** 0.5
    return data

# Calculo envolvente a la variable data que es del tipo numpy.ndarray
a = envelope(data)

# Ploteo a
plot(a)
