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
    - Se grafica el resultado para cada frecuencia que se quiere recuperar.
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
axis[1].set_title('Seno de 10 Hz')
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
S = np.fft.fft(s)/ns               # normalizo, por eso divido.


f = np.linspace(0,100,ns)          # vector frec [Hz] np.linspace(start,stop,num samples)


# gráficos
plt.figure(3)
plt.plot(f,2*np.abs(S))
plt.title("Espectro",size=20)
plt.xlabel("f [Hz]",size=15)
plt.ylabel("Amplitud [m]",size=15)


'''
############### VER GABY ###############
# Espectro de amplitud de la señal suma de senos y coseno

f = np.arange(0,ns/2+1,1)
fg = fm/f

plt.figure(3)
plt.plot(fg,2*np.abs(S[:int(ns/2+1)]))
plt.title("Espectro",size=20)
plt.xlabel("f [Hz]",size=15)
plt.ylabel("Amplitud [m]",size=15)
'''



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
plt.figure(4)
plt.plot(t,st_taper[1])
plt.title("Señal con taper Hamming",size=20)
plt.xlabel("t [seg]",size=15)
plt.ylabel("Amplitud [m]",size=15)



#------------------------------------------------------------------------------
# Aplico filtro
#------------------------------------------------------------------------------
tr = st[1]              # [1] porque por alguna razón en [0] tengo vacío
# genero una copia de t
tr_filt_10 = tr.copy()


# RECUPERO SEÑAL DE 10 HZ
# defino frec angular digital
w1_10 = 9/fm              # frec ang dig = frec[Hz]/fm
w2_10 = 11/fm

# aplico filtro pasabanda a la copia que hice anteriormente
tr_filt_10 = tr_filt_10.filter("bandpass", freqmin=w1_10, freqmax=w2_10, corners=2, zerophase=False)

# grafico la señal filtrada en tiempo
plt.figure(5)
plt.plot(t,tr_filt_10)
plt.title("Recupero 10 Hz",size=20)
plt.xlabel("t [seg]",size=15)
plt.ylabel("Amplitud [m]",size=15)


# hago FFT
ntr_filt = len(tr_filt_10)
TR_FILT_10 = np.fft.fft(tr_filt_10)/ntr_filt

# grafico el espectro de la señal filtrada y comparo con original de 10 Hz
figure, axis = plt.subplots(2, 1)
plt.subplots_adjust(hspace=1)
# Espectro señal filtrada
axis[0].set_title('Espectro de señal filtrada en 10 Hz')
#axis[0].plot(fg,2*np.abs(TR_FILT[:int(ns/2+1)]))
axis[0].plot(f,2*np.abs(TR_FILT_10))
axis[0].set_xlabel('f [Hz]')
axis[0].set_ylabel('Amplitud [m]')

# Espectro de seno de 10 Hz
ns2 = len(s2)
S2 = np.fft.fft(s2)/ns2

axis[1].set_title('Espectro de señal original de 10 Hz')
#axis[1].plot(fg,2*np.abs(S2[:int(ns2/2+1)]))
axis[1].plot(f,2*np.abs(S2))
axis[1].set_xlabel('f [Hz]')
axis[1].set_ylabel('Amplitud [m]')



# RECUPERO SEÑAL DE 5 HZ
# defino frec angular digital
w1_5 = 4/fm              # frec ang dig = frec[Hz]/fm
w2_5 = 6/fm

# aplico filtro pasabanda a la copia que hice anteriormente
tr_filt_5 = tr.copy()
tr_filt_5 = tr_filt_5.filter("bandpass", freqmin=w1_5, freqmax=w2_5, corners=2, zerophase=False)

# grafico la señal filtrada en tiempo
plt.figure(7)
plt.plot(t,tr_filt_5)
plt.title("Recupero 5 Hz",size=20)
plt.xlabel("t [seg]",size=15)
plt.ylabel("Amplitud [m]",size=15)


# hago FFT
ntr_filt = len(tr_filt_5)
TR_FILT_5 = np.fft.fft(tr_filt_5)/ntr_filt

# grafico el espectro de la señal filtrada y comparo con original de 10 Hz
fig, axis = plt.subplots(2, 1)
plt.subplots_adjust(hspace=1)
# Espectro señal filtrada
axis[0].set_title('Espectro de señal filtrada en 5 Hz')
#axis[0].plot(fg,2*np.abs(TR_FILT[:int(ns/2+1)]))
axis[0].plot(f,2*np.abs(TR_FILT_5))
axis[0].set_xlabel('f [Hz]')
axis[0].set_ylabel('Amplitud [m]')

# Espectro de seno de 5 Hz
S1 = np.fft.fft(s1)/ns2

axis[1].set_title('Espectro de señal original de 5 Hz')
#axis[1].plot(fg,2*np.abs(S1[:int(ns2/2+1)]))
axis[1].plot(f,2*np.abs(S1))
axis[1].set_xlabel('f [Hz]')
axis[1].set_ylabel('Amplitud [m]')



# RECUPERO SEÑAL DE 2 HZ
# defino frec angular digital
w1_2 = 1/fm              # frec ang dig = frec[Hz]/fm
w2_2 = 3/fm

# aplico filtro pasabanda a la copia que hice anteriormente
tr_filt_2 = tr.copy()
tr_filt_2 = tr_filt_2.filter("bandpass", freqmin=w1_2, freqmax=w2_2, corners=2, zerophase=False)

# grafico la señal filtrada en tiempo
plt.figure(9)
plt.plot(t,tr_filt_2)
plt.title("Recupero 2 Hz",size=20)
plt.xlabel("t [seg]",size=15)
plt.ylabel("Amplitud [m]",size=15)


# hago FFT
TR_FILT_2 = np.fft.fft(tr_filt_2)/ntr_filt

# grafico el espectro de la señal filtrada y comparo con original de 2 Hz
fig10, axis = plt.subplots(2, 1)
plt.subplots_adjust(hspace=1)
# Espectro señal filtrada
axis[0].set_title('Espectro de señal filtrada en 2 Hz')
#axis[0].plot(fg,2*np.abs(TR_FILT[:int(ns/2+1)]))
axis[0].plot(f,2*np.abs(TR_FILT_2))
axis[0].set_xlabel('f [Hz]')
axis[0].set_ylabel('Amplitud [m]')

# Espectro de seno de 2 Hz
S3 = np.fft.fft(s3)/ns2

axis[1].set_title('Espectro de señal original de 2 Hz')
#axis[1].plot(fg,2*np.abs(S3[:int(ns2/2+1)]))
axis[1].plot(f,2*np.abs(S3))
axis[1].set_xlabel('f [Hz]')
axis[1].set_ylabel('Amplitud [m]')