# -*- coding: utf-8 -*-
"""
Autora: María Celeste Novak Merquel

Descripción:
    1. abrir la señal 
    2. cortar la parte que quiero analizar --> evento
    3. quitar tendencia al evento
    4. quitar valor medio al punto 3
    5. aplicar al punto 4 taper hamming de 5% y 10% (graficar ambos)
    6. filtrar con pasabanda
    7. calcular fft (graficar)
"""

from obspy.core import read
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import rfft, rfftfreq
from obspy.core import Trace, Stream


'''
from obspy.core import Trace, Stream
import pandas as pd
'''

#-----------------------------------------------------------------------------
# Leo la señal que quiero analizar
#-----------------------------------------------------------------------------
st = read('C:/Users/Propietario/Desktop/tesis_de_grado/ATENUACION/DATOS/FG16/GI.FG16.00.BHZ.D.2020.064')

tr = st[0]
tr_copy = tr


#-----------------------------------------------------------------------------
# Corto la parte de la señal que me interesa analizar
#-----------------------------------------------------------------------------
inicio = UTCDateTime('2020-03-04T17:20:10')
fin = UTCDateTime('2020-03-04T17:21:25')
evento = tr_copy.slice(starttime = inicio, endtime = fin)
'''
plt.figure(1)
plt.plot(evento)
'''


#-----------------------------------------------------------------------------
# Se quita tendencia al evento
#-----------------------------------------------------------------------------
evento_sin_tendencia = evento.detrend()


#-----------------------------------------------------------------------------
# Calculo el valor medio y se lo quito al evento
#-----------------------------------------------------------------------------
val_medio_sin_tend = np.mean(evento_sin_tendencia.data)
evento_sin_medio = evento_sin_tendencia - val_medio_sin_tend


#-----------------------------------------------------------------------------
# Agrego los valores a una variable tipo trace
#-----------------------------------------------------------------------------
st_edit = Stream(Trace())
st_edit.append(Trace(data=evento_sin_medio))
tr_edit = st_edit[1]


#-----------------------------------------------------------------------------
# Aplico taper
#-----------------------------------------------------------------------------
# Taper del 5%
st_edit_taper5 = st_edit.taper(max_percentage=0.05, type='hamming', max_length=None, side='both')
tr_edit_taper5 = st_edit_taper5[1]

plt.figure(1)
plt.plot(evento, color= "black", label="evento sin taper")
plt.plot(tr_edit_taper5, color="violet", label="evento con taper del 5%")
plt.legend()


# Taper del 10%
st_edit_taper10 = st_edit.taper(max_percentage=0.1, type='hamming', max_length=None, side='both')
tr_edit_taper10 = st_edit_taper10[1]

plt.figure(2)
plt.plot(evento, color= "black", label="evento sin taper")
plt.plot(tr_edit_taper10, color="violet", label="evento con taper del 10%")
plt.legend()


#--------------
# Me quedo con taper de 10% porque no afecta a la señal
#--------------


#-----------------------------------------------------------------------------
# Aplico filtro pasabanda
#-----------------------------------------------------------------------------
# Defino frec angular digital para el filtro (frec ang dig = frec[Hz]/fm)
dt = tr.meta.delta
fm = 1/dt                   # frec de muestreo
w1 = 0.5/fm
w2 = 5/fm

evento_filtrado = tr_edit_taper10.filter("bandpass", freqmin=w1, freqmax=w2, corners=2, zerophase=False)


#-----------------------------------------------------------------------------
# Calculo FFT del evento previamente trabajado
#-----------------------------------------------------------------------------
ns = len(tr_edit_taper10)
dt = tr.meta.delta
evento_edit_fft = rfft(tr_edit_taper10)   # fft del evento con taper
f = rfftfreq(ns, dt)
evento_fft = rfft(evento_sin_medio)       # fft del evento sin taper

plt.figure(3)
plt.plot(f,np.abs(evento_fft), color="black", label="fft del evento original")
plt.plot(f,np.abs(evento_edit_fft), color="violet", label="fft del evento con taper y pasabanda")
plt.title("Espectro",size=20)
plt.xlabel("f [Hz]",size=15)
plt.ylabel("Amplitud",size=15)
plt.legend()

