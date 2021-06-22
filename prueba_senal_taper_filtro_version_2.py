# -*- coding: utf-8 -*-
"""
Título: Prueba de taper y filtro en una señal creada.

Autora: María Celeste Novak Merquel


Resumen del código:
    
1er bloque:
    - Se define una señal a partir de la suma de 2 senos y 1 coseno, cada uno con
frecuencia y amplitud diferente.
    - Se le agrega un header a la señal.
    - Se plotean cada señal por separado y también la suma.

2do bloque:
    - Se calcula la FFT se la señal.
    - Se grafica el espectro 

3er bloque:
    - Se aplica una ventana de Hamming en frecuencia.
    - Se grafica.
    
4to bloque:
    - Se aplica un filtro pasabanda.
    - Se grafica el resultado.
    --> Notar que en esta sección debo modificar manualmente las frec mín y máx
    para recuperar la señal deseada.

ACLARACIÓN: para que funcione el método .plot() de obspy es necesario que la 
salida de los gráficos sea en la misma interfaz que Spyder y no en ventanas 
aparte. (Cambiar la preferencia a EN LÍNEA)
"""

import numpy as np
import matplotlib.pyplot as plt
import obspy
from obspy.core import Trace, Stream

#------------------------------------------------------------------------------
# Defino señal en tiempo
#------------------------------------------------------------------------------

# defino parámetros
fm = 100                    # num de muestras por seg (frec de muestreo)
dt = 1/fm             
t = np.arange(0,10,dt)                                                   
a1 = 3                      # amplitud [m]
a2 = 1.5                    # amplitud [m]
a3 = 2                      # amplitud [m]
f1 = 5                      # frec [Hz]
f2 = 10                     # frec [Hz]
f3 = 2                      # frec [Hz]

# señal:
s1 = a1*np.sin(2*np.pi*f1*t)
s2 = a2*np.sin(2*np.pi*f2*t)
s3 = a3*np.cos(2*np.pi*f3*t)

s = s1 + s2 + s3


# convierto la señal en una variable tipo stream
st = Stream(Trace())  


# defino header de st
start_time = t[0] 
header = {'network': 'CE', 'station': 'VIC', 'location': 'GA',
                      'channel': 'HHZ', 'starttime': start_time}
st.append(Trace(data=s, header=header))


# gráficos en tiempo
figure, axis = plt.subplots(3, 1)
plt.subplots_adjust(hspace=1)
# s1:
axis[0].set_title('Seno de 5 Hz')
axis[0].plot(t, s1)
axis[0].set_xlabel('t [seg]')
axis[0].set_ylabel('Amplitud [m]')
# s2:
axis[1].set_title('Sino de 10 Hz')
axis[1].plot(t, s2)
axis[1].set_xlabel('t [seg]')
axis[1].set_ylabel('Amplitud [m]')
# s3:
axis[2].set_title('Coseno de 2 Hz')
axis[2].plot(t, s3)
axis[2].set_xlabel('t [seg]')
axis[2].set_ylabel('Amplitud [m]')


fig, ax = plt.subplots()
# s:
ax.plot(t,s)
ax.set_title("Señal",size=20)
ax.set_xlabel('t [seg]',size=15)
ax.set_ylabel('Amplitud [m]',size=15)



#------------------------------------------------------------------------------
# Defino señal en frecuencia
#------------------------------------------------------------------------------

# Hago zero padding
#s0 = np.pad(s, (1024,1024), 'constant') # np.pad(señal, (0 al ppio,0 al final), 'constant')


# hago fft de la señal con zero padding
#S = np.fft.fft(s)


# hago FFT
ns = len(s)
S = np.fft.fft(s)/ns          # normalizo, por eso divido. ESTPA BIEN?


f = np.linspace(0,100,ns)          # vector frec [Hz] np.linspace(start,stop,num samples)

# gráficos
plt.figure(3)
plt.plot(f,np.abs(S))
plt.title("Espectro",size=20)
plt.xlabel("f [Hz]",size=15)
plt.ylabel("Amplitud [m]",size=15)



#------------------------------------------------------------------------------
# Aplico taper
#------------------------------------------------------------------------------

# diferentes tipos de ventanas:
# https://docs.obspy.org/master/packages/autogen/obspy.core.trace.Trace.taper.html

# VER SI ES LO MISMO APLICARLO AL ST O AL TR 
# --> si lo aplico a TR me queda una variable tipo NoneType y por ende no puedo
# aplicar .plot()

st_taper = st.taper(max_percentage=0.1, type='hamming', max_length=None, side='both')

# grafico cómo queda la señal con el taper
st_taper.plot()



#------------------------------------------------------------------------------
# Aplico filtro
#------------------------------------------------------------------------------
tr = st[1]              # [1] porque por alguna razón en [0] tengo vacío
# genero una copia de t
tr_filt = tr.copy()


# RECUPERO SEÑAL DE 10 HZ
# defino frec angular digital
w1 = 9/fm              # frec ang dig = frec[Hz]/fm
w2 = 11/fm

# aplico filtro pasabanda a la copia que hice anteriormente
tr_filt.filter("bandpass", freqmin=w1, freqmax=w2, corners=2, zerophase=False)

tr_filt.plot()