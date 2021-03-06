# -*- coding: utf-8 -*-
"""
Autora: María Celeste Novak Merquel

Descripción:
    1. abrir la señal 
    2. cortar la parte que quiero analizar --> evento
    3. quitar tendencia al evento
    4. quitar valor medio al punto 3
    5. aplicar filtro pasa-altos (HP) de 1 Hz
    6. aplicar al punto 4 taper hamming de 10% y 15% (graficar ambos)
    7. filtrar con pasabanda de 1Hz. Probar con distintos rangos
    8. calcular fft (graficar) y comparar con fft de la señal sin filtrar
    9. calcular envolvente para la señal obtenida con cada filtro BP
    10. calcular valor máx de cada env y cortar c/ env desde val máx en adelante
    11. suavizar cada envolente quitando detalle
    12. calcular t para máxima amplitud de la envolvente
"""

from obspy.core import read
from obspy import UTCDateTime
import matplotlib.pyplot as plt
import numpy as np
from scipy.fft import rfft, rfftfreq, irfft
from obspy.core import Trace, Stream
import copy
import obspy
import obspy.signal
import datetime
from geopy.distance import geodesic


#%%
#-----------------------------------------------------------------------------
# 1. Leo la señal que quiero analizar
#-----------------------------------------------------------------------------
st = read('C:/Users/Propietario/Desktop/tesis_de_grado/ATENUACION/DATOS/FG16/GI.FG16.00.BHZ.D.2020.064')

tr = st[0]
tr_copy = tr
dt = tr.meta.delta

#%%
#-----------------------------------------------------------------------------
# 2. Corto la parte de la señal que me interesa analizar
#-----------------------------------------------------------------------------
inicio = UTCDateTime('2020-03-04T17:20:10')
fin = UTCDateTime('2020-03-04T17:21:10')
evento = tr_copy.slice(starttime = inicio, endtime = fin)
#evento.plot()

#%%
#-----------------------------------------------------------------------------
# 3. Se quita tendencia al evento
#-----------------------------------------------------------------------------
evento_sin_tendencia = evento.detrend()

#%%
#-----------------------------------------------------------------------------
# 4. Calculo el valor medio y se lo quito al evento
#-----------------------------------------------------------------------------
val_medio_sin_tend = np.mean(evento_sin_tendencia.data)
evento_sin_medio = evento_sin_tendencia - val_medio_sin_tend
evento_sin_medio_copy = copy.copy(evento_sin_medio)

# Gráfico en matplotlib
n = len(evento_sin_medio)
base = datetime.datetime(2020, 3, 4,17,20,10) #año, mes, día
dates = [base + datetime.timedelta(seconds=(dt * i)) for i in range(n)]
N = len(dates)
np.random.seed(19680801)
y = np.cumsum(np.random.randn(N))
#lims = (np.datetime64('2020-03-04 17:20:10'), np.datetime64('2020-03-04 17:21:10'))
lims = (np.datetime64(inicio), np.datetime64(fin))
plt.plot(dates,evento_sin_medio, color='black')
plt.xlim(lims)
titulo = str(inicio) + '-' + str(fin)
plt.title(titulo)
plt.show()

# Agrego los valores a una variable tipo trace
st_edit = Stream(Trace())
st_edit.append(Trace(data=evento_sin_medio_copy))
tr_edit = st_edit[1]

#%%
#-----------------------------------------------------------------------------
# 5. Aplico filtro pasa-alto (HP)
#-----------------------------------------------------------------------------
tr_edit_hp = tr_edit.copy()   # hago una copia para no modificarel evento original

# Defino frec angular digital para el filtro (frec ang dig = frec[Hz]/fm)
dt = tr.meta.delta
fm = 1/dt                   # frec de muestreo
w_hp = 0.1/fm
tr_edit_hp = tr_edit_hp.filter("highpass", freq=w_hp, corners=2, zerophase=False)

ns = len(tr_edit_hp)
dt = tr.meta.delta
f = rfftfreq(ns, dt)
fft_tr_hp = rfft(tr_edit_hp)   # fft de la señal con taper
fft_evento_sin_medio = rfft(evento_sin_medio_copy)   # fft de la señal evento sin taper

# Grafico para comparar
plt.figure(11)
plt.title("Comparación de espectros")
plt.plot(f,np.abs(fft_evento_sin_medio), color="black", label="señal sin HP")
plt.plot(f,np.abs(fft_tr_hp), color="violet", label="señal con HP a 0.1Hz")
plt.legend()

#%%
#-----------------------------------------------------------------------------
# 6. Aplico taper
#-----------------------------------------------------------------------------
# El pre y post evento equivalen a un 18,18% del total de la señal a analizar.
# Es por eso que si aplico un taper del 10%, éste no me va a estar afectando 
# al evento en sí, sino que va a estar afectando solo al pre y post evento 
# y todavía queda parte sin ser modificado.
# Taper del 10%
tr_edit_copy10 = tr_edit_hp.copy()   # hago una copia para no sobreescribir
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

#%%
#-----------------------------------------------------------------------------
# Calculo el espectro de la señal con taper y sin taper y comparo
#-----------------------------------------------------------------------------
# Taper del 10%
ns = len(tr_taper10)
dt = tr.meta.delta
f = rfftfreq(ns, dt)
fft_tr_taper10 = rfft(tr_taper10)   # fft de la señal con taper

# Taper del 15%
fft_tr_taper15 = rfft(tr_taper15)


fig2, (ax1, ax2) = plt.subplots(2)
fig2.suptitle('Comparaciones de espectros')
ax1.plot(f,np.abs(fft_tr_hp), color="black", label="FFT señal con HP 1Hz y sin taper")
ax1.plot(f,np.abs(fft_tr_taper10), color="violet", label="FFT señal con taper del 10%")
#ax1.set_xlim(0,1)
ax1.legend()
ax2.plot(f,np.abs(fft_tr_hp), color="black", label="FFT señal con HP 1Hz y sin taper")
ax2.plot(f,np.abs(fft_tr_taper15), color="violet", label="FFT señal con taper del 15%")
#ax2.set_xlim(0,1)
ax2.legend()

#%%
#-----------------------------------------------------------------------------
# 7. Aplico filtros pasabanda de 1 Hz
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

#%%
#-----------------------------------------------------------------------------
# 8. Calculo FFT del evento previamente trabajado
#-----------------------------------------------------------------------------
f = rfftfreq(ns, dt)
fft_01 = rfft(evento_filtrado_01)   # fft de BP e/ 0.5 - 1.5 Hz
fft_02 = rfft(evento_filtrado_02)   # fft de BP e/ 1.0 - 2.0 Hz
fft_03 = rfft(evento_filtrado_03)   # fft de BP e/ 1.5 - 2.5 Hz
fft_04 = rfft(evento_filtrado_04)   # fft de BP e/ 2.0 - 3.0 Hz
fft_05 = rfft(evento_filtrado_05)   # fft de BP e/ 2.5 - 3.5 Hz
fft_06 = rfft(evento_filtrado_06)   # fft de BP e/ 3.0 - 4.0 Hz
fft_07 = rfft(evento_filtrado_07)   # fft de BP e/ 3.5 - 4.5 Hz
fft_08 = rfft(evento_filtrado_08)   # fft de BP e/ 4.0 - 5.0 Hz

# n es el número de la figura
def graf_fft(n,fft_num,label_fft_num):
    plt.figure(n)
    a = plt.plot(f,np.abs(fft_tr_hp), color='black', label='FFT señal con HP 1Hz y sin taper')
    b = plt.plot(f,np.abs(fft_num), color='violet', label=label_fft_num)
    plt.legend(fontsize=11)
    plt.xlabel("frec [Hz]",size=11)
    plt.yticks(size=11)
    return a, b

graf_fft(3,fft_01,"FFT señal con BP 0.5 - 1.5 Hz")
graf_fft(4,fft_02,"FFT señal con BP 1.0 - 2.0 Hz")
graf_fft(5,fft_03,"FFT señal con BP 1.5 - 2.5 Hz")
graf_fft(6,fft_04,"FFT señal con BP 2.0 - 3.0 Hz")
graf_fft(7,fft_05,"FFT señal con BP 2.5 - 3.5 Hz")
graf_fft(8,fft_06,"FFT señal con BP 3.0 - 4.0 Hz")
graf_fft(9,fft_07,"FFT señal con BP 3.5 - 4.5 Hz")
graf_fft(10,fft_08,"FFT señal con BP 4.0 - 5.0 Hz")

#%%
#-----------------------------------------------------------------------------
# 9. Calculo envolvente
#-----------------------------------------------------------------------------
# función que cambia el tipo de variable de data a Stream
def new_st(data):
    header = evento.stats
    st_nueva = Stream(Trace())
    st_nueva.append(Trace(data=data, header=header))
    return st_nueva

# cálculo del vector tiempo
def vec_tiempo(variable_stream):
    npts = variable_stream[1].stats.npts
    samprate = variable_stream[1].stats.sampling_rate
    t = np.arange(0, npts/samprate, 1/samprate)
    return t

# Grafico de señal y su envolvente
def graf_env(envolvente_num,st_dom_time_num,label):
    a = plt.plot(dates, envolvente_num,'k:', label="envolvente")
    b = plt.plot(dates, st_dom_time_num[1], 'k', label=label)
    plt.xticks(size=7)
    plt.yticks(size=7)
    plt.legend(fontsize=7)
    plt.xlabel("date: 2020-03-04", size=7)
    return a, b


# BP e/ 0.5 - 1.5 Hz
dom_time_01 = irfft(fft_01,n=len(evento))     # cálculo de la transformada inversa (parte real)

# cambio el tipo de variable
st_dom_time_01 = new_st(dom_time_01)
t_01 = vec_tiempo(st_dom_time_01)

# cálculo de la envolvente
envolvente_01 = obspy.signal.filter.envelope(st_dom_time_01[1].data)

# Grafico de señal y su envolvente
plt.subplot(4, 2, 1)
graf_env(envolvente_01,st_dom_time_01,"señal con BP 0.5 - 1.5 Hz")


#----------------------------
# BP e/ 1.0 - 2.0 Hz
dom_time_02 = irfft(fft_02,n=len(evento))     # cálculo de la transformada inversa (parte real)

# cambio el tipo de variable
st_dom_time_02 = new_st(dom_time_02)
t_02 = vec_tiempo(st_dom_time_02) # vector tiempo

# cálculo de la envolvente
envolvente_02 = obspy.signal.filter.envelope(st_dom_time_02[1].data)

# Grafico de señal y su envolvente
plt.subplot(4, 2, 2)
graf_env(envolvente_02,st_dom_time_02,"señal con BP 1 Hz - 2 Hz")


#----------------------------
# BP e/ 1.5 - 2.5 Hz
dom_time_03 = irfft(fft_03,n=len(evento))     # cálculo de la transformada inversa (parte real)

# cambio el tipo de variable
st_dom_time_03 = new_st(dom_time_03)
t_03 = vec_tiempo(st_dom_time_03) # vector tiempo

# cálculo de la envolvente
envolvente_03 = obspy.signal.filter.envelope(st_dom_time_03[1].data)

# Grafico de señal y su envolvente
plt.subplot(4, 2, 3)
graf_env(envolvente_03,st_dom_time_03,"señal con BP 1.5 - 2.5 Hz")


#----------------------------
# BP e/ 2.0 - 3.0 Hz
dom_time_04 = irfft(fft_04,n=len(evento))     # cálculo de la transformada inversa (parte real)

# cambio el tipo de variable
st_dom_time_04 = new_st(dom_time_04)
t_04 = vec_tiempo(st_dom_time_04) # vector tiempo

# cálculo de la envolvente
envolvente_04 = obspy.signal.filter.envelope(st_dom_time_04[1].data)

# Grafico de señal y su envolvente
plt.subplot(4, 2, 4)
graf_env(envolvente_04,st_dom_time_04,"señal con BP 2.0 - 3.0 Hz")


#----------------------------
# BP e/ 2.5 - 3.5 Hz
dom_time_05 = irfft(fft_05,n=len(evento))     # cálculo de la transformada inversa (parte real)

# cambio el tipo de variable
st_dom_time_05 = new_st(dom_time_05)
t_05 = vec_tiempo(st_dom_time_05) # vector tiempo

# cálculo de la envolvente
envolvente_05 = obspy.signal.filter.envelope(st_dom_time_05[1].data)

# Grafico de señal y su envolvente
plt.subplot(4, 2, 5)
graf_env(envolvente_05,st_dom_time_05,"señal con BP 2.5 - 3.5 Hz")


#----------------------------
# BP e/ 3.0 - 4.0 Hz
dom_time_06 = irfft(fft_06,n=len(evento))     # cálculo de la transformada inversa (parte real)

# cambio el tipo de variable
st_dom_time_06 = new_st(dom_time_06)
t_06 = vec_tiempo(st_dom_time_06) # vector tiempo

# cálculo de la envolvente
envolvente_06 = obspy.signal.filter.envelope(st_dom_time_06[1].data)

# Grafico de señal y su envolvente
plt.subplot(4, 2, 6)
graf_env(envolvente_06,st_dom_time_06,"señal con BP 3.0 - 4.0 Hz")


#----------------------------
# BP e/ 3.5 - 4.5 Hz
dom_time_07 = irfft(fft_07,n=len(evento))     # cálculo de la transformada inversa (parte real)

# cambio el tipo de variable
st_dom_time_07 = new_st(dom_time_07)
t_07 = vec_tiempo(st_dom_time_07) # vector tiempo

# cálculo de la envolvente
envolvente_07 = obspy.signal.filter.envelope(st_dom_time_07[1].data)

# Grafico de señal y su envolvente
plt.subplot(4, 2, 7)
graf_env(envolvente_07,st_dom_time_07,"señal con BP 3.5 - 4.5 Hz")


#----------------------------
# BP e/ 4.0 - 5.0 Hz
dom_time_08 = irfft(fft_08,n=len(evento))     # cálculo de la transformada inversa (parte real)

# cambio el tipo de variable
st_dom_time_08 = new_st(dom_time_08)
t_08 = vec_tiempo(st_dom_time_08) # vector tiempo

# cálculo de la envolvente
envolvente_08 = obspy.signal.filter.envelope(st_dom_time_08[1].data)

# Grafico de señal y su envolvente
plt.subplot(4, 2, 8)
graf_env(envolvente_08,st_dom_time_08,"señal con BP 4.0 - 5.0 Hz")

#%%
#-----------------------------------------------------------------------------
# 10. Valor máximo de cada envolvente y corte de la misma
#-----------------------------------------------------------------------------
def posicion_max(envolvente_num):
    val = np.amax(envolvente_num)
    lista = envolvente_num.tolist()
    posicion = lista.index(val)
    return posicion

max_01 = posicion_max(envolvente_01)
max_02 = posicion_max(envolvente_02)
max_03 = posicion_max(envolvente_03)
max_04 = posicion_max(envolvente_04)
max_05 = posicion_max(envolvente_05)
max_06 = posicion_max(envolvente_06)
max_07 = posicion_max(envolvente_07)
max_08 = posicion_max(envolvente_08)

cut_env01 = envolvente_01[max_01:]
cut_env02 = envolvente_02[max_02:]
cut_env03 = envolvente_03[max_03:]
cut_env04 = envolvente_04[max_04:]
cut_env05 = envolvente_05[max_05:]
cut_env06 = envolvente_06[max_06:]
cut_env07 = envolvente_07[max_07:]
cut_env08 = envolvente_08[max_08:]

cut_t01 = dates[max_01:]
cut_t02 = dates[max_02:]
cut_t03 = dates[max_03:]
cut_t04 = dates[max_04:]
cut_t05 = dates[max_05:]
cut_t06 = dates[max_06:]
cut_t07 = dates[max_07:]
cut_t08 = dates[max_08:]

#%%
#-----------------------------------------------------------------------------
# 11. Suavizamiento de la envolvente para quitar detalle
#-----------------------------------------------------------------------------
from scipy.signal import savgol_filter
s_env01 = savgol_filter(cut_env01, 71, 3) # window size 51, polynomial order 3
s_env02 = savgol_filter(cut_env02, 71, 3)
s_env03 = savgol_filter(cut_env03, 71, 3)
s_env04 = savgol_filter(cut_env04, 71, 3)
s_env05 = savgol_filter(cut_env05, 71, 3)
s_env06 = savgol_filter(cut_env06, 71, 3)
s_env07 = savgol_filter(cut_env07, 71, 3)
s_env08 = savgol_filter(cut_env08, 71, 3)

'''
plt.figure(50)
plt.plot(dates, st_dom_time_01[1], 'k', label="señal")
plt.plot(cut_t01,cut_env01, 'k:', label="envolvente sin suavizar")
plt.plot(cut_t01,s_env01, 'violet', label="envolvente suavizada")
lims = (np.datetime64(inicio), np.datetime64(fin))
plt.xlim(lims)
plt.title("Suavizado con Savitzky-Golay filter (SciPy)")
plt.legend()
plt.xlabel("date: 2020-03-04")
'''
'''
op01 = obspy.signal.util.smooth(cut_env01, 31)
# OBS: cuanto mayor es el número, más baja los máximos de la env
plt.figure(51)
plt.plot(t_01, st_dom_time_01[1], 'k', label="señal")
plt.plot(cut_t01, cut_env01, 'k:', label="envolvente sin suavizar")
plt.plot(cut_t01, op01, 'violet', label="envolvente suavizada")
plt.title("Suavizado con central moving average (ObsPy)")
plt.legend()
plt.xlabel("tiempo [seg] ????")
'''

# CONCLUSIÓN: la mejor opción parecería ser savgol_filter()


# Grafico de señal y su envolvente
def graf_env_suav(s_envnum,cut_envnum,cut_tnum,st_dom_time_num,label):
    plt.suptitle("Suavizado con Savitzky-Golay filter (SciPy)")
    a = plt.plot(dates, st_dom_time_num[1], 'k', label=label)
    b = plt.plot(cut_tnum, cut_envnum, 'k:', label="envolvente sin suavizar")
    c = plt.plot(cut_tnum, s_envnum,'violet', label="envolvente suavizada")
    plt.xticks(size=7)
    plt.yticks(size=7)
    plt.legend(fontsize=7)
    plt.xlabel("date: 2020-03-04", size=7)

    return a, b, c

fig = plt.subplots(4,2)

plt.subplot(4, 2, 1)
graf_env_suav(s_env01,cut_env01,cut_t01,st_dom_time_01,"señal con BP 0.5 - 1.5 Hz")

plt.subplot(4, 2, 2)
graf_env_suav(s_env02,cut_env02,cut_t02,st_dom_time_02,"señal con BP 1 Hz - 2 Hz")

plt.subplot(4, 2, 3)
graf_env_suav(s_env03,cut_env03,cut_t03,st_dom_time_03,"señal con BP 1.5 - 2.5 Hz")

plt.subplot(4, 2, 4)
graf_env_suav(s_env04,cut_env04,cut_t04,st_dom_time_04,"señal con BP 2.0 - 3.0 Hz")

plt.subplot(4, 2, 5)
graf_env_suav(s_env05,cut_env05,cut_t05,st_dom_time_05,"señal con BP 2.5 - 3.5 Hz")

plt.subplot(4, 2, 6)
graf_env_suav(s_env06,cut_env06,cut_t06,st_dom_time_06,"señal con BP 3.0 - 4.0 Hz")

plt.subplot(4, 2, 7)
graf_env_suav(s_env07,cut_env07,cut_t07,st_dom_time_07,"señal con BP 3.5 - 4.5 Hz")

plt.subplot(4, 2, 8)
graf_env_suav(s_env08,cut_env08,cut_t08,st_dom_time_08,"señal con BP 4.0 - 5.0 Hz")

#%%
#-----------------------------------------------------------------------------
# 12. Cálculo del tiempo correspondiente al máx valor de la envolvente
#-----------------------------------------------------------------------------
# me devuelve el número de la siguiente forma:
# datetime.datetime(año, mes, día, hora, min, seg, VER)
t_max01 = dates[max_01]
t_max02 = dates[max_02]
t_max03 = dates[max_03]
t_max04 = dates[max_04]
t_max05 = dates[max_05]
t_max06 = dates[max_06]
t_max07 = dates[max_07]
t_max08 = dates[max_08]


#%%
#-----------------------------------------------------------------------------
# Cálculo de distancia en km entre cráter y estación 
#-----------------------------------------------------------------------------
# buscar como esta proyectando para hacer el calculo
FG_8 = (14.43250,-90.93590)
FG_10 = (14.41300,-90.91239)  
FG_12 = (14.43651,-90.83606) 
FG_13 = (14.40677,-90.81859)
FG_16 = (14.455,-90.8508)
crater = (14.474585,-90.880777)

r_fg8 = geodesic(FG_8, crater).kilometers
r_fg10 = geodesic(FG_10, crater).kilometers
r_fg12 = geodesic(FG_12, crater).kilometers
r_fg13 = geodesic(FG_13, crater).kilometers
r_fg16 = geodesic(FG_16, crater).kilometers




#%%
#-----------------------------------------------------------------------------
# Ajuste exponencial
#-----------------------------------------------------------------------------
# uso ScyPy
from scipy.optimize import curve_fit

def func(x, a, b, c):
    return a * np.exp(-b * x) + c

m = len(s_env01)
x = np.arange(0,m*dt,dt)
y = s_env01.copy()

# popt: Valores óptimos para los parámetros de modo que la suma del error al 
# cuadrado de f (xdata, * popt) - ydata se minimice
# pcov: La covarianza estimada de popt. Las diagonales proporcionan la
# varianza de la estimación del parámetro. Para calcular un error de
# desviación estándar en los parámetros, use perr = np.sqrt (np.diag (pcov)).
popt, pcov = curve_fit(func, x, y) 
#plt.plot(x, func(x, *popt))

a, b, c = tuple(popt)

# copio los parámetros de la exponencial en un csv
import pandas as pd
#parametros_ajuste = pd.DataFrame([[a, b, c]], columns=['a', 'b', 'c'])
parametros_ajuste = pd.DataFrame([tuple(popt)],columns=['a', 'b', 'c'])

parametros_ajuste.to_csv('parametros_ajuste.csv')

# leo el csv que hice arriba y printeo
import csv
with open('parametros_ajuste.csv', newline='') as f:
    reader = csv.reader(f)
    for row in reader:
        print(row)



'''
plt.figure(20)
plt.plot(x,s_env01)
plt.plot(x, func(x, *popt))
'''
plt.figure(30)
plt.suptitle("Ajuste exponencial de la envolvente con SciPy")
a = plt.plot(dates, st_dom_time_01[1], 'k', label= "señal con BP 0.5 - 1.5 Hz")
b = plt.plot(cut_t01, cut_env01, 'g:', label="envolvente sin suavizar")
c = plt.plot(cut_t01, s_env01,'violet', label="envolvente suavizada")
d = plt.plot(cut_t01, func(x, *popt), label="ajuste exponencial")
plt.xticks(size=7)
plt.yticks(size=7)
plt.legend(fontsize=7)
plt.xlabel("date: 2020-03-04", size=7)




#%%
#-----------------------------------------------------------------------------
# Gráfico para poner en la tesis. Evento, espectrograma y espectro de ampl
#-----------------------------------------------------------------------------

from obspy import read
from obspy.signal import util
from obspy.imaging.cm import obspy_sequential
def cwt(st, dt, w0, fmin, fmax, nf=100, wl='morlet'):
    """
    Continuous Wavelet Transformation in the Frequency Domain.

    .. seealso:: [Kristekova2006]_, eq. (4)

    :param st: time dependent signal.
    :param dt: time step between two samples in st (in seconds)
    :param w0: parameter for the wavelet, tradeoff between time and frequency
        resolution
    :param fmin: minimum frequency (in Hz)
    :param fmax: maximum frequency (in Hz)
    :param nf: number of logarithmically spaced frequencies between fmin and
        fmax
    :param wl: wavelet to use, for now only 'morlet' is implemented

    :return: time frequency representation of st, type numpy.ndarray of complex
        values, shape = (nf, len(st)).
    """
    npts = len(st) * 2
    tmax = (npts - 1) * dt
    t = np.linspace(0., tmax, npts)
    f = np.logspace(np.log10(fmin), np.log10(fmax), nf)

    cwt = np.zeros((npts // 2, nf), dtype=complex)

    if wl == 'morlet':

        def psi(t):
            return np.pi ** (-.25) * np.exp(1j * w0 * t) * \
                np.exp(-t ** 2 / 2.)

        def scale(f):
            return w0 / (2 * np.pi * f)
    else:
        raise ValueError('wavelet type "' + wl + '" not defined!')

    nfft = util.next_pow_2(npts) * 2
    sf = np.fft.fft(st, n=nfft)

    # Ignore underflows.
    with np.errstate(under="ignore"):
        for n, _f in enumerate(f):
            a = scale(_f)
            # time shift necessary, because wavelet is defined around t = 0
            psih = psi(-1 * (t - t[-1] / 2.) / a).conjugate() / np.abs(a) ** .5
            psihf = np.fft.fft(psih, n=nfft)
            tminin = int(t[-1] / 2. / (t[1] - t[0]))
            cwt[:, n] = np.fft.ifft(psihf * sf)[tminin:tminin + npts // 2] * \
                (t[1] - t[0])

    return cwt.T

def _pcolormesh_same_dim(ax, x, y, v, **kwargs):
    # x, y, v must have the same dimension
    try:
        return ax.pcolormesh(x, y, v, shading='nearest', **kwargs)
    except TypeError:
        # matplotlib versions < 3.3
        return ax.pcolormesh(x, y, v[:-1, :-1], **kwargs)

def plot_tfr(st, dt=0.01, t0=0., fmin=1., fmax=10., nf=100, w0=6, left=0.1,
             bottom=0.1, h_1=0.2, h_2=0.6, w_1=0.2, w_2=0.6, w_cb=0.01,
             d_cb=0.0, show=True, plot_args=['k', 'k'], clim=0.0,
             cmap=obspy_sequential, mode='absolute', fft_zero_pad_fac=0):
    
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    npts = st.shape[-1]
    tmax = (npts - 1) * dt
    t = np.linspace(0., tmax, npts) + t0

    if fft_zero_pad_fac == 0:
        nfft = npts
    else:
        nfft = util.next_pow_2(npts) * fft_zero_pad_fac

    f_lin = np.linspace(0, 0.5 / dt, nfft // 2 + 1)

    if len(st.shape) == 1:
        _w = np.zeros((1, nf, npts), dtype=complex)
        _w[0] = cwt(st, dt, w0, fmin, fmax, nf)
        ntr = 1

        spec = np.zeros((1, nfft // 2 + 1), dtype=complex)
        spec[0] = np.fft.rfft(st, n=nfft) * dt

        st = st.reshape((1, npts))
    else:
        _w = np.zeros((st.shape[0], nf, npts), dtype=complex)
        spec = np.zeros((st.shape[0], nfft // 2 + 1), dtype=complex)

        for i in np.arange(st.shape[0]):
            _w[i] = cwt(st[i], dt, w0, fmin, fmax, nf)
            spec[i] = np.fft.rfft(st[i], n=nfft) * dt

        ntr = st.shape[0]

    if mode == 'absolute':
        _tfr = np.abs(_w)
        spec = np.abs(spec)
    elif mode == 'power':
        _tfr = np.abs(_w) ** 2
        spec = np.abs(spec) ** 2
    else:
        raise ValueError('mode "' + mode + '" not defined!')

    figs = []

    for itr in np.arange(ntr):
        fig = plt.figure()

        # plot signals
        ax_sig = fig.add_axes([left + w_1, bottom, w_2, h_1])
        ax_sig.plot(t, st[itr], plot_args[0])

        # plot TFR
        ax_tfr = fig.add_axes([left + w_1, bottom + h_1, w_2, h_2])

        x, y = np.meshgrid(
            t, np.logspace(np.log10(fmin), np.log10(fmax),
                           _tfr[itr].shape[0]))
        img_tfr = _pcolormesh_same_dim(ax_tfr, x, y, _tfr[itr], cmap=cmap)
        img_tfr.set_rasterized(True)
        ax_tfr.set_yscale("log")
        ax_tfr.set_ylim(fmin, fmax)
        ax_tfr.set_xlim(t[0], t[-1])

        # plot spectrum
        ax_spec = fig.add_axes([left, bottom + h_1, w_1, h_2])
        ax_spec.semilogy(spec[itr], f_lin, plot_args[1])

        # add colorbars
        ax_cb_tfr = fig.add_axes([left + w_1 + w_2 + d_cb + w_cb, bottom +
                                  h_1, w_cb, h_2])
        fig.colorbar(img_tfr, cax=ax_cb_tfr)

        # set limits
        ax_sig.set_ylim(st.min() * 1.1, st.max() * 1.1)
        ax_sig.set_xlim(t[0], t[-1])

        xlim = spec.max() * 1.1

        ax_spec.set_xlim(xlim, 0.)
        ax_spec.set_ylim(fmin, fmax)

        if clim == 0.:
            clim = _tfr.max()

        img_tfr.set_clim(0., clim)

        ax_sig.set_xlabel('Tiempo [seg] después de 2020-03-04T17:20:10.000000Z')
        ax_spec.set_ylabel('Frecuencia [Hz]')

        # remove axis labels
        ax_tfr.xaxis.set_major_formatter(NullFormatter())
        ax_tfr.yaxis.set_major_formatter(NullFormatter())

        figs.append(fig)

    if show:
        plt.show()
    else:
        if ntr == 1:
            return figs[0]
        else:
            return figs



if __name__ == '__main__':
    import doctest
    doctest.testmod(exclude_empty=True)
    
    #log set_cmap
#from obspy.signal.tf_misfit import plot_tfr
#tr = read("https://examples.obspy.org/a02i.2008.240.mseed")[0]
plot_tfr(tr_taper10.data, dt=evento.stats.delta, t0=0, fmin=.01, 
        fmax=50., w0=6, nf=64, fft_zero_pad_fac=4)



