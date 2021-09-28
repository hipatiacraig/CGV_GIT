# -*- coding: utf-8 -*-
"""
Autora: María Celeste Novak Merquel

Descripción:
    1. abrir la señal 
    2. cortar la parte que quiero analizar --> evento
    3. quitar tendencia al evento
    4. quitar valor medio al punto 3
    5. aplicar al punto 4 taper hamming de 10% y 15% (graficar ambos)
    6. filtrar con pasabanda de 1Hz. Probar con distintos rangos
    7. calcular fft (graficar) y comparar con fft de la señal sin filtrar
"""

from obspy.core import read
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import rfft, rfftfreq
from obspy.core import Trace, Stream
import copy

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
fin = UTCDateTime('2020-03-04T17:21:10')
evento = tr_copy.slice(starttime = inicio, endtime = fin)
evento.plot()


#-----------------------------------------------------------------------------
# Se quita tendencia al evento
#-----------------------------------------------------------------------------
evento_sin_tendencia = evento.detrend()


#-----------------------------------------------------------------------------
# Calculo el valor medio y se lo quito al evento
#-----------------------------------------------------------------------------
val_medio_sin_tend = np.mean(evento_sin_tendencia.data)
evento_sin_medio = evento_sin_tendencia - val_medio_sin_tend
evento_sin_medio_copy = copy.copy(evento_sin_medio)



#-----------------------------------------------------------------------------
# Agrego los valores a una variable tipo trace
#-----------------------------------------------------------------------------
st_edit = Stream(Trace())
st_edit.append(Trace(data=evento_sin_medio_copy))
tr_edit = st_edit[1]


#-----------------------------------------------------------------------------
# Aplico taper
#-----------------------------------------------------------------------------
# El pre y post evento equivalen a un 18,18% del total de la señal a analizar.
# Es por eso que si aplico un taper del 10%, éste no me va a estar afectando 
# al evento en sí, sino que va a estar afectando solo al pre y post evento 
# y todavía queda parte sin ser modificado.
# Taper del 10%
tr_edit_copy10 = tr_edit.copy()   # hago una copia para no sobreescribir
tr_taper10 = tr_edit_copy10.taper(max_percentage=0.1, type='hamming', max_length=None, side='both')

# Taper del 15%
tr_edit_copy15 = tr_edit.copy()   # hago una copia para no sobreescribir
tr_taper15 = tr_edit_copy15.taper(max_percentage=0.15, type='hamming', max_length=None, side='both')

# Grafico para comparar
fig, (ax1, ax2) = plt.subplots(2)
fig.suptitle('Comparaciones')
ax1.plot(evento_sin_medio, color="black", label="señal sin tendencia ni valor medio")
ax1.plot(tr_taper10, color="violet", label="señal con taper del 10%")
ax1.legend()
ax2.plot(evento_sin_medio, color="black", label="señal sin tendencia ni valor medio")
ax2.plot(tr_taper15, color="violet", label="señal con taper del 15%")
ax2.legend()


#-----------------------------------------------------------------------------
# Calculo el espectro de la señal con taper y sin taper y comparo
#-----------------------------------------------------------------------------
# Taper del 10%
ns = len(tr_taper10)
dt = tr.meta.delta
f = rfftfreq(ns, dt)
fft_tr_taper10 = rfft(tr_taper10)   # fft de la señal con taper
fft_evento_sin_medio = rfft(evento_sin_medio_copy)   # fft de la señal evento sin taper

# Taper del 15%
fft_tr_taper15 = rfft(tr_taper15)


fig2, (ax1, ax2) = plt.subplots(2)
fig2.suptitle('Comparaciones de espectros')
ax1.plot(f,np.abs(fft_evento_sin_medio), color="black", label="FFT señal sin taper")
ax1.plot(f,np.abs(fft_tr_taper10), color="violet", label="FFT señal con taper del 10%")
#ax1.set_xlim(0,1)
ax1.legend()
ax2.plot(f,np.abs(fft_evento_sin_medio), color="black", label="FFT señal sin taper")
ax2.plot(f,np.abs(fft_tr_taper15), color="violet", label="FFT señal con taper del 15%")
#ax2.set_xlim(0,1)
ax2.legend()


#-----------------------------------------------------------------------------
# Aplico filtros pasabanda de 1 Hz
#-----------------------------------------------------------------------------
# Defino frec angular digital para el filtro (frec ang dig = frec[Hz]/fm)

# 01
dt = tr.meta.delta
fm = 1/dt                   # frec de muestreo
w1 = 0.5/fm
w2 = 1.5/fm

tr_taper10_copy01 = tr_taper10.copy()   # hago una copia para no sobreescribir
evento_filtrado_01 = tr_taper10_copy01.filter("bandpass", freqmin=w1, freqmax=w2, corners=2, zerophase=False)


# ----------------------------
# 02
w3 = 1.0/fm
w4 = 2.0/fm
tr_taper10_copy02 = tr_taper10.copy()   # hago una copia para no sobreescribir
evento_filtrado_02 = tr_taper10_copy02.filter("bandpass", freqmin=w3, freqmax=w4, corners=2, zerophase=False)

# ----------------------------
# 03
w5 = 1.5/fm
w6 = 2.5/fm
tr_taper10_copy03 = tr_taper10.copy()   # hago una copia para no sobreescribir
evento_filtrado_03 = tr_taper10_copy03.filter("bandpass", freqmin=w5, freqmax=w6, corners=2, zerophase=False)

# ----------------------------
# 04
w7 = 2.0/fm
w8 = 3.0/fm
tr_taper10_copy04 = tr_taper10.copy()   # hago una copia para no sobreescribir
evento_filtrado_04 = tr_taper10_copy04.filter("bandpass", freqmin=w7, freqmax=w8, corners=2, zerophase=False)

# ----------------------------
# 05
w9 = 2.5/fm
w10 = 3.5/fm
tr_taper10_copy05 = tr_taper10.copy()   # hago una copia para no sobreescribir
evento_filtrado_05 = tr_taper10_copy05.filter("bandpass", freqmin=w9, freqmax=w10, corners=2, zerophase=False)

# ----------------------------
# 06
w11 = 3.0/fm
w12 = 4.0/fm
tr_taper10_copy06 = tr_taper10.copy()   # hago una copia para no sobreescribir
evento_filtrado_06 = tr_taper10_copy06.filter("bandpass", freqmin=w11, freqmax=w12, corners=2, zerophase=False)

# ----------------------------
# 07
w13 = 3.5/fm
w14 = 4.5/fm
tr_taper10_copy07 = tr_taper10.copy()   # hago una copia para no sobreescribir
evento_filtrado_07 = tr_taper10_copy07.filter("bandpass", freqmin=w13, freqmax=w14, corners=2, zerophase=False)

# ----------------------------
# 08
w15 = 4.0/fm
w16 = 5.0/fm
tr_taper10_copy08 = tr_taper10.copy()   # hago una copia para no sobreescribir
evento_filtrado_08 = tr_taper10_copy08.filter("bandpass", freqmin=w15, freqmax=w16, corners=2, zerophase=False)


#-----------------------------------------------------------------------------
# Calculo FFT del evento previamente trabajado
#-----------------------------------------------------------------------------
fft_01 = rfft(evento_filtrado_01)   # fft de BP e/ 0.5 - 1.5 Hz
f = rfftfreq(ns, dt)

fft_02 = rfft(evento_filtrado_02)   # fft de BP e/ 1.0 - 2.0 Hz

fft_03 = rfft(evento_filtrado_03)   # fft de BP e/ 1.5 - 2.5 Hz

fft_04 = rfft(evento_filtrado_04)   # fft de BP e/ 2.0 - 3.0 Hz

fft_05 = rfft(evento_filtrado_05)   # fft de BP e/ 2.5 - 3.5 Hz

fft_06 = rfft(evento_filtrado_06)   # fft de BP e/ 3.0 - 4.0 Hz

fft_07 = rfft(evento_filtrado_07)   # fft de BP e/ 3.5 - 4.5 Hz

fft_08 = rfft(evento_filtrado_08)   # fft de BP e/ 4.0 - 5.0 Hz

plt.figure(3)
plt.plot(f,np.abs(fft_evento_sin_medio), color="black", label="FFT señal sin taper")
plt.plot(f,np.abs(fft_01), color="violet", label="FFT señal con BP 0.5 - 1.5 Hz")
plt.legend()


plt.figure(4)
plt.plot(f,np.abs(fft_evento_sin_medio), color="black", label="FFT señal sin taper")
plt.plot(f,np.abs(fft_02), color="violet", label="FFT señal con BP 1.0 - 2.0 Hz")
plt.legend()

plt.figure(5)
plt.plot(f,np.abs(fft_evento_sin_medio), color="black", label="FFT señal sin taper")
plt.plot(f,np.abs(fft_03), color="violet", label="FFT señal con BP 1.5 - 2.5 Hz")
plt.legend()

plt.figure(6)
plt.plot(f,np.abs(fft_evento_sin_medio), color="black", label="FFT señal sin taper")
plt.plot(f,np.abs(fft_04), color="violet", label="FFT señal con BP 2.0 - 3.0 Hz")
plt.legend()

plt.figure(7)
plt.plot(f,np.abs(fft_evento_sin_medio), color="black", label="FFT señal sin taper")
plt.plot(f,np.abs(fft_05), color="violet", label="FFT señal con BP 2.5 - 3.5 Hz")
plt.legend()

plt.figure(8)
plt.plot(f,np.abs(fft_evento_sin_medio), color="black", label="FFT señal sin taper")
plt.plot(f,np.abs(fft_06), color="violet", label="FFT señal con BP 3.0 - 4.0 Hz")
plt.legend()

plt.figure(9)
plt.plot(f,np.abs(fft_evento_sin_medio), color="black", label="FFT señal sin taper")
plt.plot(f,np.abs(fft_07), color="violet", label="FFT señal con BP 4.0 - 5.0 Hz")
plt.legend()

plt.figure(10)
plt.plot(f,np.abs(fft_evento_sin_medio), color="black", label="FFT señal sin taper")
plt.plot(f,np.abs(fft_08), color="violet", label="FFT señal con BP 3.5 - 4.5 Hz")
plt.legend()
