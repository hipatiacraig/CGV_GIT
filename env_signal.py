# -*- coding: utf-8 -*-

import numpy as np
from scipy import signal
from scipy.integrate import cumtrapz
from obspy.core import read                   # Importo para leer

#from . import util
from obspy.signal import util

from pylab import *                           # Para graficar


st = read('C:/Users/Propietario/Desktop/tesis_de_grado/ATENUACION/datos_prueba/Trazas_Fuego/GI.FG16.00.BHZ.D.2020.218')

# Lo paso a tipo trace
tr = st[0]

# Datos a los que le voy a aplicar la función envolvente
# (son tipo numpy.ndarray)????
trdat = tr.data

def envelope(data):
    """
    Envelope of a signal.
    Computes the envelope of the given data which can be windowed or
    not. The envelope is determined by the absolute value of the analytic
    signal of the given data.
    
    If data are windowed the analytic signal and the envelope of each
    window is returned.
    
    :type data: :class:`~numpy.ndarray`
    :param data: Data to make envelope of.
    :return: **A_cpx, A_abs** - Analytic signal of input data, Envelope of
        input data.
    """
    nfft = util.next_pow_2(data.shape[-1]) # valor de pot de 2 más cercana a
                                           # data.shape[-1]
    a_cpx = np.zeros((data.shape), dtype=np.complex64)
    a_abs = np.zeros((data.shape), dtype=np.float64)
    if len(data.shape) > 1:
        i = 0
        for row in data:
            a_cpx[i, :] = signal.hilbert(row, nfft)
            a_abs[i, :] = abs(signal.hilbert(row, nfft))
            i = i + 1
    else:
        a_cpx = signal.hilbert(data, nfft)
        a_abs = abs(signal.hilbert(data, nfft))
    return a_cpx, a_abs

# Calculo envolvente para los datos que quiero --> trdat que ya estarían en la
# forma numpy.ndarray
b = envelope(trdat)

# Grafico lo que me interesa --> la parte del valor absoluto
plot(b[1])
