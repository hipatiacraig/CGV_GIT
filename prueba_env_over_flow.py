# -*- coding: utf-8 -*-
'''
Código adapado al que está en https://stackoverflow.com/questions/37968221/similar-function-envelope-matlab-in-python
Las líneas comentadas con # están así ya que en nuestro caso cargamos el dato 
de interés y no es que lo "creamos" para hacer la prueba.
La parte de la frecuencia y fase también las comenté. --> VER SI ES NECESARIO
O NO CALCULARLO
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert, chirp
from obspy.core import read 

#duration = 1.0
#fs = 400.0
#samples = int(fs*duration)
#t = np.arange(samples) / fs

# Dato de interés
st = read('d:/Victoria/Desktop/PHD/FUEGO_INSIVUMEH/DATOS_05052021/FG8/GI.FG8.00.BHE.D.2020.061')

# Convierto steam a trace
tr = st[0]
# Dato de tipo numpy.ndarray
data = tr.data

analytic_data = hilbert(data)
amplitude_envelope = np.abs(analytic_data)
#instantaneous_phase = np.unwrap(np.angle(analytic_data))
#instantaneous_frequency = np.diff(instantaneous_phase) / (2.0*np.pi) * fs

# Grafico
fig = plt.figure()
ax0 = fig.add_subplot(211)
ax0.plot(data, label='signal')
ax0.plot(amplitude_envelope, label='envelope')
ax0.set_xlim(6000,10000)
ax0.set_ylim(-10000,+10000)
ax0.legend()


#ax0.set_xlabel("time in seconds")
#ax0.legend()
#ax1 = fig.add_subplot(212)
#ax1.plot(t[1:], instantaneous_frequency)
#ax1.set_xlabel("time in seconds")
#ax1.set_ylim(0.0, 120.0)
