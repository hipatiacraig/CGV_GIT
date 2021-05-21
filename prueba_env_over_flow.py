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

st = read('C:/Users/Propietario/Desktop/tesis_de_grado/ATENUACION/datos_prueba/Trazas_Fuego/GI.FG16.00.BHZ.D.2020.218')

tr = st[0]
data = tr.data

analytic_data = hilbert(data)
amplitude_envelope = np.abs(analytic_data)
#instantaneous_phase = np.unwrap(np.angle(analytic_data))
#instantaneous_frequency = np.diff(instantaneous_phase) / (2.0*np.pi) * fs

fig = plt.figure()
ax0 = fig.add_subplot(211)
ax0.plot(data, label='signal')
ax0.plot(amplitude_envelope, label='envelope')
ax0.legend()


#ax0.set_xlabel("time in seconds")
#ax0.legend()
#ax1 = fig.add_subplot(212)
#ax1.plot(t[1:], instantaneous_frequency)
#ax1.set_xlabel("time in seconds")
#ax1.set_ylim(0.0, 120.0)
