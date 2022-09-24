# -*- coding: utf-8 -*-
"""
Autora: María Celeste Novak Merquel

Descripción:
    1. Se abre el archivo de tags y se busca si ya existe el archivo de salida 
    2. Se localizan en el archivo de entrada los eventos que se van a procesar
 NO 3. Se quitar tendencia al evento
    4. Se calcula y elimina el valor medio de cada evento
    5. Se aplica filtro pasa-alto (HP) de 1 Hz
    6. Se aplicar al dato taper hamming de 10%
    7. Se filtra con pasabandas de ancho 1Hz
    8. Se calcula fft del resultado en 7. Se compara original con filtrado
    9. Cálculo de envolvente para la señal obtenida con cada filtro BP
    10. Cálculo de valor máx de cada envolvente y corte de la misma
 ?? 11. Se suaviza cada envolente quitando detalle
    12. Se calcula el tiempo correspondiente al máx valor de la envolvente
    13. Se calcula de distancia en km entre cráter y estación 
    14. Se calcula el tiempo entre la fuente y la estación FG16
    15. Se calculan los vectores lapsetime
    16. Se corrige por divergencia esférica
    17. Primer ajuste: lineal
    18. Ajuste exponencial a partir de los resultados de 17
    
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
import pandas as pd
import os
from geopy.distance import geodesic
from sklearn.linear_model import LinearRegression
from scipy.optimize import curve_fit
import csv


#%%
#-----------------------------------------------------------------------------
# 1. Se abre el archivo de tags y se busca si ya existe el archivo de salida
#-----------------------------------------------------------------------------

### path del file donde se encuentran los tags
#tags_dir='/media/gaby/Backup Plus/SISMOLOGIA_DATOS_Y_DOCS/Vn_FUEGO/tags'
tags_dir='C:/Users/Propietario/Desktop/tesis_de_grado/ATENUACION/DATOS/FG16'

### tags_dir + 'nombre del archivo que contiene los tags'
#tags_file=tags_dir +'/fecha_tipo_dur_ampl_max_LP_2020_MEJORADO.csv'
tags_file=tags_dir + '/fecha_tipo_dur_ampl_max_LP_2020_MEJORADO_1.csv'

#data_dir='/media/gaby/Backup Plus/SISMOLOGIA_DATOS_Y_DOCS/Vn_FUEGO/FG16'
data_dir='C:/Users/Propietario/Desktop/tesis_de_grado/ATENUACION/DATOS/FG16'

### los eventos a leer son de tipo LP. El tag de un evento comienza con LP y
# finaliza con Noise
EventCode='LP'
EventEnd = 'Noise'
tags=pd.read_csv(tags_file,sep=',',header=0)

"""
Ojo que el csv debe tener en primera fila los titulos de columnas:"fecha,scnl,tipo" 
para que el read_csv  pueda cargar una dataframe con las tres columnas ya tituladas
""" 
dates=tags.fecha
scnl=tags.scnl
events=tags.tipo

### path del directorio donde se guardará un archivo con todos los resultados
dire_salida='C:/Users/Propietario/Desktop/tesis_de_grado/ATENUACION/DATOS/FG16/resultados'
#dire_salida='/media/gaby/Backup Plus/SISMOLOGIA_DATOS_Y_DOCS/Vn_FUEGO/resultados'

existencias=os.listdir(dire_salida)
scnl1 = scnl[0]
red = scnl1[9:11]
archivo_salida = red + '_Q.csv'
path_archivo_salida = dire_salida +'/' + archivo_salida

### si ya existe un archivo de salida creado por el programa, éste lo detectará
# y avisará por terminal que ya existe. En caso contrario, avisará que no 
# existe y lo crea imprimiendo ya los titulos de columnas
if len(existencias)==0: #si el listado de archivos resulta vacio
    existe_archivo=0
    print('No existe archivo',archivo_salida)
    field_names  = ['UTCDateTime','Fc','lapsetime','dist','p_dg','vel','a_lin','b_lin','a_exp','b_exp','Qc_SBS','a1_dif','a2_dif','a3_dif','Qi','d','Qs'] #The name of the columns
    with open(path_archivo_salida, 'w') as f:
        writer = csv.writer(f) #this is the writer object
        writer.writerow(field_names ) # this will list out the names of the columns which are always the first entrries
else:
    for i in range(0,len(existencias),1): # si el listado de archivos no esta vacio
        if existencias[i]==archivo_salida: #ya hay un archivo con ese nombre
            existe_archivo=1
            print('Ya existe archivo...',archivo_salida)
            break
        else:
            existe_archivo=0 #no hay un archivo con ese nombre
            print('No existe archivo',archivo_salida)
            field_names  = ['UTCDateTime','Fc','lapsetime','dist','p_dg','vel','a_lin','b_lin','a_exp','b_exp','Qc_SBS','a1_dif','a2_dif','a3_dif','Qi','d','Qs'] #The name of the columns
            with open(path_archivo_salida, 'w') as f:
                writer = csv.writer(f) #this is the writer object
                writer.writerow(field_names ) # this will list out the names of the columns which are always the first entrries


#%%
#-----------------------------------------------------------------------------
# 2. Se localizan en el archivo de entrada los eventos que se van a procesar
#-----------------------------------------------------------------------------

### días que se analizarán. Deben estar en un archivo txt en donde la primer 
# fila tiene el título y la segunda fila tiene los nombres de los archivos a 
# analizar; los nombres de los archivos tienen que estar separados por: espacio
# coma espacio
selectos_file = data_dir+'/selectos.csv'
#selectos_file = tags_dir+'/selectos.csv'

### se lee el csv
selectosDF = pd.read_csv(selectos_file,sep=',',header=0)

### se guarda la información en una dataframe
dias_selectos = selectosDF.selectos

### a los eventos se les agrega un pre y post evento para que al momento de
# procesar no se pierda información importante del evento en sí.
'''
Revisar pre y post para distintos eventos
'''
pre_evento = 10 #seg
post_evento = 10 #seg

### se procesan todos los días seleccionados
for i in range(0,len(dias_selectos)):   
    # dia_procesar es el i-esimo de la seleccion de dias hecha por Vic y Cel 
    dia_procesar = dias_selectos[i]  #me da la fecha
    dia_procesar = dia_procesar[:10]  #leo los primeros 10 caracteres de i

        
    #print("dia_procesar",dia_procesar)
    ### dates es la fecha del tag
    for j in range(0,len(dates),1):  #range(inicio,fin,paso)
        # dates_procesar es el j-esimo  del listado de tags
        dates_procesar = dates[j]
        dates_procesar = dates_procesar[:10]
       
        
        if dia_procesar == dates_procesar:  #me fijo cuándo aparece la fecha que quiero
            fila_inicio = j  #guardo la posición en la que está el dia_procesar en el archivo con todos los tags
            #print(fila_inicio)
            #print("dates_procesar",dates_procesar)
            #inicio_string = dates[j]
            #fin_string = dates[j+1]
            #inicio = UTCDateTime(inicio_string)
            tipo_inicio = events[j]
            #fin = UTCDateTime(fin_string)
            tipo_fin = events[j+1]
            
            if tipo_inicio == 'LP' and tipo_fin == 'Noise': # confirmo que tiene inicio y fin
                break
    #print('fila inicio',j,'inicio_string',inicio_string)            
    #print('fila fin',j+1,'fin_string',fin_string)
    for k in range(j,len(dates),2):
        ### se recorren todos los tags del dia_procesar
        date_comparar = dates[k]
        date_comparar = date_comparar[:10]
        if date_comparar == dia_procesar:
            inicio_string=dates[k]
            fin_string=dates[k+1]
            inicio=UTCDateTime(inicio_string)
            fin=UTCDateTime(fin_string)
            ### se agrega pre y postevento en segundos
            inicio_pre = inicio - pre_evento
            fin_post = fin + post_evento
            
            ### se generan los vectores vacíos que se van a llenar con los 
            # parámetros del ajuste exponencial para el cálculo de Q y con el 
            # tiempo de inicio del tag para luego guardarlos en el archivo de salida
            if i==0 and existe_archivo==0:
                print('generando vectores')
                a=[]
                b=[]
                c=[]
                tiempos=[]
            
            ### se extrae año y día juliano de inicio_pre y fin_post para 
            # convertirlos en cadenas de caracteres y usarlos para completar
            # el nombre de/de los wave_file/s
            anio=[]
            dia=[]
            hora=[]
            minuto=[]
            segundo=[]
            anio.append(inicio_pre.year)
            dia.append(inicio_pre.julday)
            hora.append(inicio_pre.hour)
            minuto.append(inicio_pre.minute)
            segundo.append(inicio_pre.second)
            
            anio.append(fin_post.year)
            dia.append(fin_post.julday)
            hora.append(fin_post.hour)
            minuto.append(fin_post.minute)
            segundo.append(fin_post.second)
 
            extremos={'anio':anio,'jday':dia,'hour':hora,'minute':minuto,'second':segundo}
            DFextremos=pd.DataFrame(extremos)
            DFextremos['anio']=DFextremos['anio'].apply(str)
            DFextremos['jday']=DFextremos['jday'].apply(str)
            DFextremos['hour']=DFextremos['hour'].apply(str)
            DFextremos['minute']=DFextremos['minute'].apply(str)
            DFextremos['second']=DFextremos['second'].apply(str)
            
            ### se arma una lista con los nombres de archivos de onda de los 
            # días de inicio y fin
            wave_file=[]
            nombre = scnl[k]

            for m in range(0,len(anio),1):
                if len(DFextremos.jday[m])==1:
                    dayofyear_str='00'+DFextremos.jday[m]
                elif len(DFextremos.jday[m])==2:
                    dayofyear_str='0'+DFextremos.jday[m]
                else:
                    dayofyear_str=DFextremos.jday[m]
                
                wave_file_raiz  = data_dir + '/' + nombre[9:11] + '.' + nombre[0:4]+'.'+nombre[12:14]+'.'+nombre[5:8]+'.'+'D.'               
                wave_file.append(wave_file_raiz + DFextremos.anio[m]+'.'+dayofyear_str)

### se controla si los dos nombres de archivos son los mismos y 
# se carga una o dos ondas. Puede pasar que parte de un evento esté
# en el final de un archivo y la otra parque haya quedado en el
# inicio del otro.
# En el segundo caso se hace merge de las dos formas de onda antes
            # de extraer evento
            if wave_file[0]==wave_file[1]:
                wave=read(wave_file[0])
            else:
                wave=read(wave_file[0])
                wave +=read(wave_file[1])
                wave.sort(['starttime'])
                wave.merge(method=1)
            
            ### se extrae evento
            evento=wave.slice(inicio_pre,fin_post)
            #evento.plot()
    
            #continue
            
            
#%%
            '''
#-----------------------------------------------------------------------------
# 3. Se quita tendencia al evento
# se saca esto por ahora porque a trazas cortas les imprime una pendiente
#-----------------------------------------------------------------------------
            evento_sin_tendencia = evento.detrend()
            '''
            
            
#%%
#-----------------------------------------------------------------------------
# 4. Se calcula y elimina el valor medio de cada evento
#-----------------------------------------------------------------------------
            
            traza_evento=evento[0]
            dt = traza_evento.meta.delta # seg
            dato=traza_evento.data
            val_medio_dato=np.mean(dato)
            dato_sin_medio=dato-val_medio_dato
            traza_evento_sin_medio=traza_evento.copy() 
            traza_evento_sin_medio.data=dato_sin_medio
            #traza_evento_sin_medio.plot()
            

#-----------------------------------------------------------------------------
# 5. Se aplica filtro pasa-alto (HP)
#-----------------------------------------------------------------------------
            
### se hace una copia para no modificar el evento original y sobreescribirlo
            tr_edit_hp = traza_evento_sin_medio.copy()   
            
### se defino frec angular digital para el filtro (frec ang dig = frec[Hz]/fm)
            dt = traza_evento_sin_medio.meta.delta
            fm = 1/dt                   # frec de muestreo
            w_hp = 0.1/fm
            tr_edit_hp = tr_edit_hp.filter("highpass", freq=w_hp, corners=2, zerophase=False)
            ns = len(tr_edit_hp)
            f = rfftfreq(ns, dt)
            fft_tr_hp = rfft(tr_edit_hp)   # fft de la señal con taper
            
            '''
            # SI QUIERO COMPARAR ALGUNA TRAZA CON Y SIN FILTRO, DESCOMENTAR:
            fft_evento_sin_medio = rfft(dato)   # fft de la señal evento sin taper
            plt.figure(11)
            plt.title("Comparación de espectros de amplitud")
            plt.plot(f,np.abs(fft_evento_sin_medio), color="black", label="Datos sin filtrar")
            plt.plot(f,np.abs(fft_tr_hp), color="violet", label="Datos filtrados con pasa alto de 0.002 Hz")
            plt.xlabel("frecuencia [Hz]")
            plt.ylabel("amplitud")
            plt.legend()
            '''
            ######## CALCULAR ESTO PARA OTROS EVENTOS

                        
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

            '''
            # Taper del 15%b
            tr_edit_copy15 = traza_evento_sin_medio.copy()   # hago una copia para no sobreescribir
            tr_taper15 = tr_edit_copy15.taper(max_percentage=0.15, type='hamming', max_length=None, side='both')
            
            # SI QUIERO COMPARAR ALGUNA TRAZA CON Y SIN FILTRO, DESCOMENTAR:
            # Grafico para comparar
            fig, (ax1, ax2) = plt.subplots(2)
            fig.suptitle('Comparaciones')
            ax1.plot(dato, color="black", label="señal sin valor medio")
            ax1.plot(tr_taper10, color="violet", label="señal con taper del 10%")
            ax1.legend()
            ax2.plot(dato, color="black", label="señal sin valor medio")
            ax2.plot(tr_taper15, color="violet", label="señal con taper del 15%")
            ax2.legend()
            '''
#%%
#-----------------------------------------------------------------------------
# Calculo el espectro de la señal con taper y sin taper y comparo
#-----------------------------------------------------------------------------
# Taper del 10%
            ns = len(tr_taper10)
            f = rfftfreq(ns, dt)
            fft_tr_taper10 = rfft(tr_taper10)   # fft de la señal con taper
            
            '''
            # Taper del 15%
            fft_tr_taper15 = rfft(tr_taper15)
            
            # SI QUIERO COMPARAR ALGUNA TRAZA CON Y SIN FILTRO, DESCOMENTAR:
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
            '''
            
#%%
#-----------------------------------------------------------------------------
# 7. Aplico filtros pasabanda de 1 Hz
#-----------------------------------------------------------------------------
           
            # 01
            fm = 1/dt                   # frec de muestreo
            #w1 = 0.5/fm
            #w2 = 1.5/fm
            w1=0.5
            w2=1.5
            tr_taper10_copy01 = tr_taper10.copy()   # hago una copia para no sobreescribir
            evento_filtrado_01 = tr_taper10_copy01.filter("bandpass", freqmin=w1, freqmax=w2, corners=2, zerophase=False)
            
            # ----------------------------
            # 02
            # w3 = 1.0/fm
            # w4 = 2.0/fm
            w3 = 1.0
            w4 = 2.0
            tr_taper10_copy02 = tr_taper10.copy()   # hago una copia para no sobreescribir
            evento_filtrado_02 = tr_taper10_copy02.filter("bandpass", freqmin=w3, freqmax=w4, corners=2, zerophase=False)
            
            # ----------------------------
            # 03
            #w5 = 1.5/fm
            #w6 = 2.5/fm
            w5=1.5
            w6=2.5
            tr_taper10_copy03 = tr_taper10.copy()   # hago una copia para no sobreescribir
            evento_filtrado_03 = tr_taper10_copy03.filter("bandpass", freqmin=w5, freqmax=w6, corners=2, zerophase=False)
            
            # ----------------------------
            # 04
            #w7 = 2.0/fm
            #w8 = 3.0/fm
            w7=2.0
            w8=3.0
            tr_taper10_copy04 = tr_taper10.copy()   # hago una copia para no sobreescribir
            evento_filtrado_04 = tr_taper10_copy04.filter("bandpass", freqmin=w7, freqmax=w8, corners=2, zerophase=False)
            
            # ----------------------------
            # 05
            #w9 = 2.5/fm
            #w10 = 3.5/fm
            w9=2.5
            w10=3.5
            tr_taper10_copy05 = tr_taper10.copy()   # hago una copia para no sobreescribir
            evento_filtrado_05 = tr_taper10_copy05.filter("bandpass", freqmin=w9, freqmax=w10, corners=2, zerophase=False)
            
            # ----------------------------
            # 06
            #w11 = 3.0/fm
            #w12 = 4.0/fm
            w11=3.0
            w12=4.0
            tr_taper10_copy06 = tr_taper10.copy()   # hago una copia para no sobreescribir
            evento_filtrado_06 = tr_taper10_copy06.filter("bandpass", freqmin=w11, freqmax=w12, corners=2, zerophase=False)
            
            # ----------------------------
            # 07
            #w13 = 3.5/fm
            #w14 = 4.5/fm
            w13=3.5
            w14=4.5
            tr_taper10_copy07 = tr_taper10.copy()   # hago una copia para no sobreescribir
            evento_filtrado_07 = tr_taper10_copy07.filter("bandpass", freqmin=w13, freqmax=w14, corners=2, zerophase=False)
            
            # ----------------------------
            # 08
            #w15 = 4.0/fm
            #w16 = 5.0/fm
            w15=4.0
            w16=5.0
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
                a = plt.plot(f,np.abs(fft_tr_hp), color='black', label='FFT señal con HP 0.002Hz y sin taper')
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
# 9. Cálculo de envolvente
#-----------------------------------------------------------------------------

# función que cambia el tipo de variable de data a Stream
            def new_st(data):
                header = evento[0].stats
                st_nueva = Stream(Trace())
                st_nueva.append(Trace(data=data, header=header))
                return st_nueva
            
# cálculo del vector tiempo
            def vec_tiempo(variable_stream):
                npts = variable_stream[0].stats.npts
                samprate = variable_stream[0].stats.sampling_rate
                t = np.arange(0, npts/samprate, 1/samprate)
                return t
            
# gráfico de señal y su envolvente
            def graf_env(t_num,envolvente_num,st_dom_time_num,label):
                a = plt.plot(t_num, envolvente_num,'k:', label="envolvente")
                b = plt.plot(t_num, st_dom_time_num[0], 'k', label=label)
                plt.xticks(size=7)
                plt.yticks(size=7)
                plt.legend(fontsize=7)
               # plt.xlabel(inicio_string, size=7)
                return a, b
            
            
            # BP e/ 0.5 - 1.5 Hz
            dom_time_01 = irfft(fft_01,n=len(evento[0]))     # cálculo de la transformada inversa (parte real)
            
            # cambio el tipo de variable
            st_dom_time_01 = (new_st(dom_time_01)).slice(inicio,fin)
            t_01 = vec_tiempo(st_dom_time_01)
            
            # cálculo de la envolvente
            envolvente_01 = obspy.signal.filter.envelope(st_dom_time_01[0].data)
            
            # Grafico de señal y su envolvente
            plt.subplot(4, 2, 1)
            titu=inicio_string+"señal con BP 0.5 - 1.5 Hz"
            graf_env(t_01,envolvente_01,st_dom_time_01,titu)
            
            
            #----------------------------
            # BP e/ 1.0 - 2.0 Hz
            dom_time_02 = irfft(fft_02,n=len(evento[0]))     # cálculo de la transformada inversa (parte real)
            
            # cambio el tipo de variable
            st_dom_time_02 = (new_st(dom_time_02)).slice(inicio,fin)
            t_02 = vec_tiempo(st_dom_time_02) # vector tiempo
            
            # cálculo de la envolvente
            envolvente_02 = obspy.signal.filter.envelope(st_dom_time_02[0].data)
            
            # Grafico de señal y su envolvente
            plt.subplot(4, 2, 2)
            titu=inicio_string+"señal con BP 1 - 2 Hz"
            graf_env(t_02,envolvente_02,st_dom_time_02,titu)
            
            
            #----------------------------
            # BP e/ 1.5 - 2.5 Hz
            dom_time_03 = irfft(fft_03,n=len(evento[0]))     # cálculo de la transformada inversa (parte real)
            
            # cambio el tipo de variable
            st_dom_time_03 = (new_st(dom_time_03)).slice(inicio,fin)
            t_03 = vec_tiempo(st_dom_time_03) # vector tiempo
            
            # cálculo de la envolvente
            envolvente_03 = obspy.signal.filter.envelope(st_dom_time_03[0].data)
            
            # Grafico de señal y su envolvente
            plt.subplot(4, 2, 3)
            titu=inicio_string+"señal con BP 1.5 - 2.5 Hz"
            graf_env(t_03,envolvente_03,st_dom_time_03,titu)
            
                
            #----------------------------
            # BP e/ 2.0 - 3.0 Hz
            dom_time_04 = irfft(fft_04,n=len(evento[0]))     # cálculo de la transformada inversa (parte real)
            
            # cambio el tipo de variable
            st_dom_time_04 = (new_st(dom_time_04)).slice(inicio,fin)
            t_04 = vec_tiempo(st_dom_time_04) # vector tiempo
            
            # cálculo de la envolvente
            envolvente_04 = obspy.signal.filter.envelope(st_dom_time_04[0].data)
            
            # Grafico de señal y su envolvente
            plt.subplot(4, 2, 4)
            titu=inicio_string+"señal con BP 2.0 - 3.0 Hz"
            graf_env(t_04,envolvente_04,st_dom_time_04,titu)
            
            
            #----------------------------
            # BP e/ 2.5 - 3.5 Hz
            dom_time_05 = irfft(fft_05,n=len(evento[0]))     # cálculo de la transformada inversa (parte real)
            
            # cambio el tipo de variable
            st_dom_time_05 = (new_st(dom_time_05)).slice(inicio,fin)
            t_05 = vec_tiempo(st_dom_time_05) # vector tiempo
            
            # cálculo de la envolvente
            envolvente_05 = obspy.signal.filter.envelope(st_dom_time_05[0].data)
            
            # Grafico de señal y su envolvente
            plt.subplot(4, 2, 5)
            titu=inicio_string+"señal con BP 2.5 - 3.5 Hz"
            graf_env(t_05,envolvente_05,st_dom_time_05,titu)
            
            
            #----------------------------
            # BP e/ 3.0 - 4.0 Hz
            dom_time_06 = irfft(fft_06,n=len(evento[0]))     # cálculo de la transformada inversa (parte real)
            
            # cambio el tipo de variable
            st_dom_time_06 = (new_st(dom_time_06)).slice(inicio,fin)
            t_06 = vec_tiempo(st_dom_time_06) # vector tiempo
            
            # cálculo de la envolvente
            envolvente_06 = obspy.signal.filter.envelope(st_dom_time_06[0].data)
            
            # Grafico de señal y su envolvente
            plt.subplot(4, 2, 6)
            titu=inicio_string+"señal con BP 3.0 - 4.0 Hz"
            graf_env(t_06,envolvente_06,st_dom_time_06,titu)
            
            
            #----------------------------
            # BP e/ 3.5 - 4.5 Hz
            dom_time_07 = irfft(fft_07,n=len(evento[0]))     # cálculo de la transformada inversa (parte real)
            
            # cambio el tipo de variable
            st_dom_time_07 = (new_st(dom_time_07)).slice(inicio,fin)
            t_07 = vec_tiempo(st_dom_time_07) # vector tiempo
            
            # cálculo de la envolvente
            envolvente_07 = obspy.signal.filter.envelope(st_dom_time_07[0].data)
            
            # Grafico de señal y su envolvente
            plt.subplot(4, 2, 7)
            titu=inicio_string+"señal con BP 3.5 - 4.5 Hz"
            graf_env(t_07,envolvente_07,st_dom_time_07,titu)
            
            
            #----------------------------
            # BP e/ 4.0 - 5.0 Hz
            dom_time_08 = irfft(fft_08,n=len(evento[0]))     # cálculo de la transformada inversa (parte real)
            
            # cambio el tipo de variable
            st_dom_time_08 = (new_st(dom_time_08)).slice(inicio,fin)
            t_08 = vec_tiempo(st_dom_time_08) # vector tiempo
            
            # cálculo de la envolvente
            envolvente_08 = obspy.signal.filter.envelope(st_dom_time_08[0].data)
            
            # Grafico de señal y su envolvente
            plt.subplot(4, 2, 8)
            titu=inicio_string+"señal con BP 4.0 - 5.0 Hz"
            graf_env(t_08,envolvente_08,st_dom_time_08,titu)
            
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
            
            
# corto la envolvente y me quedo con la parte entre el valor máximo
# y la última muestra antes del post_evento            
            cut_env01 = envolvente_01[max_01:]
            cut_env02 = envolvente_02[max_02:]
            cut_env03 = envolvente_03[max_03:]
            cut_env04 = envolvente_04[max_04:]
            cut_env05 = envolvente_05[max_05:]
            cut_env06 = envolvente_06[max_06:]
            cut_env07 = envolvente_07[max_07:]
            cut_env08 = envolvente_08[max_08:]

# ídem anterior, pero con el vector tiempo
            cut_t01 = t_01[max_01:]
            cut_t02 = t_02[max_02:]
            cut_t03 = t_03[max_03:]
            cut_t04 = t_04[max_04:]
            cut_t05 = t_05[max_05:]
            cut_t06 = t_06[max_06:]
            cut_t07 = t_07[max_07:]
            cut_t08 = t_08[max_08:]
            
#%%
#-----------------------------------------------------------------------------
# 11. Suavizamiento de la envolvente para quitar detalle
#-----------------------------------------------------------------------------
            from scipy.signal import savgol_filter
            s_env01 = savgol_filter(cut_env01, 151, 3) # window size 51, polynomial order 3
            s_env02 = savgol_filter(cut_env02, 151, 3)
            s_env03 = savgol_filter(cut_env03, 151, 3)
            s_env04 = savgol_filter(cut_env04, 151, 3)
            s_env05 = savgol_filter(cut_env05, 151, 3)
            s_env06 = savgol_filter(cut_env06, 151, 3)
            s_env07 = savgol_filter(cut_env07, 151, 3)
            s_env08 = savgol_filter(cut_env08, 151, 3)
            
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
            
            
# Gráfico de señal y su envolvente
            def graf_env_suav(s_envnum,cut_envnum,cut_tnum,st_dom_time_num,t_num,label):
                titu2 = inicio_string + " - Suavizado con Savitzky-Golay filter (SciPy)"
                plt.suptitle(titu2)
                a = plt.plot(t_num, st_dom_time_num[0], 'k', label=label)
                b = plt.plot(cut_tnum, cut_envnum, 'k:', label="envolvente sin suavizar")
                c = plt.plot(cut_tnum, s_envnum,'violet', label="envolvente suavizada")
                plt.xticks(size=7)
                plt.yticks(size=7)
                plt.legend(fontsize=7)
                #plt.xlabel("date: 2020-03-04", size=7)
            
                return a, b, c
            
            fig = plt.subplots(4,2)
            
            plt.subplot(4, 2, 1)
            graf_env_suav(s_env01,cut_env01,cut_t01,st_dom_time_01,t_01,"señal con BP 0.5 - 1.5 Hz")
            
            plt.subplot(4, 2, 2)
            graf_env_suav(s_env02,cut_env02,cut_t02,st_dom_time_02,t_02,"señal con BP 1 Hz - 2 Hz")
            
            plt.subplot(4, 2, 3)
            graf_env_suav(s_env03,cut_env03,cut_t03,st_dom_time_03,t_03,"señal con BP 1.5 - 2.5 Hz")
            
            plt.subplot(4, 2, 4)
            graf_env_suav(s_env04,cut_env04,cut_t04,st_dom_time_04,t_04,"señal con BP 2.0 - 3.0 Hz")
            
            plt.subplot(4, 2, 5)
            graf_env_suav(s_env05,cut_env05,cut_t05,st_dom_time_05,t_05,"señal con BP 2.5 - 3.5 Hz")
            
            plt.subplot(4, 2, 6)
            graf_env_suav(s_env06,cut_env06,cut_t06,st_dom_time_06,t_06,"señal con BP 3.0 - 4.0 Hz")
            
            plt.subplot(4, 2, 7)
            graf_env_suav(s_env07,cut_env07,cut_t07,st_dom_time_07,t_07,"señal con BP 3.5 - 4.5 Hz")
            
            plt.subplot(4, 2, 8)
            graf_env_suav(s_env08,cut_env08,cut_t08,st_dom_time_08,t_08,"señal con BP 4.0 - 5.0 Hz")
            
#%%
#-----------------------------------------------------------------------------
# 12. Cálculo del tiempo correspondiente al máx valor de la envolvente
#-----------------------------------------------------------------------------
# me devuelve el número de la siguiente forma:
# datetime.datetime(año, mes, día, hora, min, seg, VER)
            t_max01 = inicio_pre + t_01[max_01]
            t_max02 = inicio_pre + t_01[max_02]
            t_max03 = inicio_pre + t_01[max_03]
            t_max04 = inicio_pre + t_01[max_04]
            t_max05 = inicio_pre + t_01[max_05]
            t_max06 = inicio_pre + t_01[max_06]
            t_max07 = inicio_pre + t_01[max_07]
            t_max08 = inicio_pre + t_01[max_08]
            
            
#%%
#-----------------------------------------------------------------------------
# 13. Cálculo de distancia en km entre cráter y estación 
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
            r_fg16 = geodesic(FG_16, crater).kilometers # dist e/ cráter y GF16
            
            

#%%
#-----------------------------------------------------------------------------
# 14. Cálculo del tiempo que le lleva a la señal ir desde la fuente hasta la
# estación FG16
#-----------------------------------------------------------------------------

            vel = 3.5  #km/seg
            t_evento = r_fg16/vel #desde cráter a fg16
            

#%%
#-----------------------------------------------------------------------------
# 15. Cálculo de vectores lapsetime
#-----------------------------------------------------------------------------
# Se calcula el tiempo real de cómo se dio el evento. Ya que cuando se 
# extrae el evento para poder procesarlo se le asocia un vector tiempo de dura-
# ción genérico, no teniendo en cuenta el tiempo que tardó el LP en viajar
# de la fuente a la estación, etc.
            delta_max01 = t_max01 - inicio
            lapsetime01 = t_evento + delta_max01 + cut_t01
            
            delta_max02 = t_max02 - inicio
            lapsetime02 = t_evento + delta_max02 + cut_t02
            
            delta_max03 = t_max03 - inicio
            lapsetime03 = t_evento + delta_max03 + cut_t03
            
            delta_max04 = t_max04 - inicio
            lapsetime04 = t_evento + delta_max04 + cut_t04
            
            delta_max05 = t_max05 - inicio
            lapsetime05 = t_evento + delta_max05 + cut_t05
            
            delta_max06 = t_max06 - inicio
            lapsetime06 = t_evento + delta_max06 + cut_t06

            delta_max07 = t_max07 - inicio
            lapsetime07 = t_evento + delta_max07 + cut_t07
            
            delta_max08 = t_max08 - inicio
            lapsetime08 = t_evento + delta_max08 + cut_t08

#%%
#-----------------------------------------------------------------------------
# 16. Corrección por divergencia esférica
#-----------------------------------------------------------------------------

# comentados están las correcciones calculadas con la envolvente suavizada
# ojo porque algunos valores de la env suavizada son negativos!!!!         


            Alpt01 = s_env01 * lapsetime01
            
            Alpt02 = s_env02 * lapsetime02
            
            Alpt03 = s_env03 * lapsetime03
            
            Alpt04 = s_env04 * lapsetime04
            
            Alpt05 = s_env05 * lapsetime05
            
            Alpt06 = s_env06 * lapsetime06
            
            Alpt07 = s_env07 * lapsetime07
            
            Alpt08 = s_env08 * lapsetime08
            
            '''
# a modo de prueba, cálculo de la corrección con la envolvente sin suavizar
            pAlpt01 = cut_env01 * lapsetime01
            
            pAlpt02 = cut_env02 * lapsetime02
            
            pAlpt03 = cut_env03 * lapsetime03
            
            pAlpt04 = cut_env04 * lapsetime04
            
            pAlpt05 = cut_env05 * lapsetime05
            
            pAlpt06 = cut_env06 * lapsetime06
            
            pAlpt07 = cut_env07 * lapsetime07
            
            pAlpt08 = cut_env08 * lapsetime08
            '''

#%%
#-----------------------------------------------------------------------------
# 17. Primer ajuste: lineal
#-----------------------------------------------------------------------------       
# Se hace un ajuste lineal, los datos de entrada serían el tiempo real (no el
# del registo) y el log de las amplitudes.
# El objetivo de este paso es buscar valores para "dárselos" al ajuste
# exponencial de la envolvente, que se hará después, así tiene un rango de 
# valores por donde comenzar a ajustar.

            ##########################################################
            # NUEVO CORTADO
            '''
            CAMBIAR EXPONENTE EN CORREC GEOM
            '''
            # ----------------------------
            # 01
            # se hace el ajute lineal a los datos filtrados por BP 0.5 - 1.5 Hz            
            print("Resultados del ajuste lineal. Frecuencias entre 0.5 - 1.5 Hz")
            
            largo = 150
            
            lpt01_grow = lapsetime01[:largo] ### acepto como mínimo 30 segundo a 50 mps
            s_env01_grow = s_env01[:largo]
            
            Alpt01 = s_env01_grow * lpt01_grow
            # se hace reshape de los vectores para que cumplan la condición que
            # pide el paquete de regresión lineal
            logAlpt01 = np.log(Alpt01).reshape(-1)
            lpt01_grow_reshape = lpt01_grow.reshape(-1,1)
            model01 = LinearRegression().fit(lpt01_grow_reshape, logAlpt01)
            prediccion01 = model01.predict(lpt01_grow_reshape)
            r_sq01 = model01.score(lpt01_grow_reshape, logAlpt01)
            #print("Resultados del ajuste lineal. Frecuencias entre 0.5 - 1.5 Hz")
            #print(f"coefficient of determination: {r_sq01}")
            p2_01_a = model01.intercept_
            #print("intersección: {:.4f}".format(p2_01))
            p1_01_a = model01.coef_
            #print("pendiente: {:.4f}. len(lapsetime01)={:}".format(p1_01[0],len(lapsetime01)))
            
            #largo = largo + 10*50
            #lpt01_grow = lapsetime01[:largo]
            #s_env01_grow = s_env01[:largo]
            
            while len(lpt01_grow) < len(lapsetime01):

                largo = largo + 10*50
                lpt01_grow = lapsetime01[:largo]
                s_env01_grow = s_env01[:largo]                

                Alpt01 = s_env01_grow * lpt01_grow
                logAlpt01 = np.log(Alpt01).reshape(-1)
                lpt01_grow_reshape = lpt01_grow.reshape(-1,1)
                model01 = LinearRegression().fit(lpt01_grow_reshape, logAlpt01)
                prediccion01 = model01.predict(lpt01_grow_reshape)
                r_sq01 = model01.score(lpt01_grow_reshape, logAlpt01)
                #print("Resultados del ajuste lineal. Frecuencias entre 0.5 - 1.5 Hz")
                #print(f"coefficient of determination: {r_sq01}")
                p2_01_b = model01.intercept_
                #print("intersección: {:.4f}".format(p2_01))
                p1_01_b = model01.coef_
                #print("pendiente: {:.4f}. len(lapsetime01)={:}".format(p1_01[0],len(lapsetime01)))
                #print("largo: {}, pendiente: {}".format(largo,p1_01_b[0]))
                
                if p1_01_b[0] > p1_01_a[0]: #>
                    
                    p1_01 = p1_01_a[0]
                    p2_01 = p2_01_a
                    #le resto 10*50 porque significa q me sirve la pendiente anterior
                    #entonces tengo que volver los arrays al lago de ese momento
                    len_lpt01 = len(lpt01_grow) - (10*50)
                    lpt01_grow = lapsetime01[:len_lpt01]
                    print("len_lpt01: {}, pendiente: {}".format(len_lpt01,p1_01))
                    break
                
                elif p1_01_b[0] < p1_01_a[0]: #<
                    #print("pendiente: {:.4f}. len(lapsetime01)={:}".format(p1_01[0],len(lapsetime01)))
                    #largo = largo + 10*50
                    lpt01_grow = lapsetime01[:largo]
                    s_env01_grow = s_env01[:largo]
                    
                    p1_01_a[0] = p1_01_b[0] 

                    
            cut_prediccion01 = prediccion01[:len_lpt01]
            cut_Alpt01 = Alpt01[:len_lpt01]
            cut_cut_env01 = cut_env01[:len_lpt01]
            cut_logAlpt01 = logAlpt01[:len_lpt01]

            plt.figure(1)
            plt.plot(lpt01_grow,cut_prediccion01, label="Alpt = {:.4f}t + {:.4f}".format(p1_01,p2_01))
            #plt.plot(lpt01_grow,cut_Alpt01, label="exponencial corregida por div geom")
            plt.plot(lpt01_grow,cut_logAlpt01, label="logaritmo de la amplitud")
            #plt.plot(lpt01_grow,cut_cut_env01, label="envolvente de la coda")
            plt.xlabel("tiempo [seg]")
            plt.ylabel("amplitud corregida por div geom")
            plt.title("LAPSETIME DETERMINADO CON LOOP")
            plt.legend()    
            


            # ----------------------------
            # 02
            print("Resultados del ajuste lineal. Frecuencias entre 1.0 - 2.0 Hz")
            
            largo = 150
            
            lpt02_grow = lapsetime02[:largo] ### acepto como mínimo 30 segundo a 50 mps
            s_env02_grow = s_env02[:largo]
            
            Alpt02 = s_env02_grow * lpt02_grow

            logAlpt02 = np.log(Alpt02).reshape(-1)
            lpt02_grow_reshape = lpt02_grow.reshape(-1,1)
            model02 = LinearRegression().fit(lpt02_grow_reshape, logAlpt02)
            prediccion02 = model02.predict(lpt02_grow_reshape)
            r_sq02 = model02.score(lpt02_grow_reshape, logAlpt02)
            p2_02_a = model02.intercept_
            #print("intersección: {:.4f}".format(p2_02))
            p1_02_a = model02.coef_
            #print("pendiente: {:.4f}. len(lpt02_grow)={:}".format(p1_02[0],len(lpt02_grow)))
            
            
            while len(lpt02_grow) < len(lapsetime02):

                largo = largo + 10*50
                lpt02_grow = lapsetime02[:largo]
                s_env02_grow = s_env02[:largo]                

                Alpt02 = s_env02_grow * lpt02_grow
                logAlpt02 = np.log(Alpt02).reshape(-1)
                lpt02_grow_reshape = lpt02_grow.reshape(-1,1)
                model02 = LinearRegression().fit(lpt02_grow_reshape, logAlpt02)
                prediccion02 = model02.predict(lpt02_grow_reshape)
                r_sq02 = model02.score(lpt02_grow_reshape, logAlpt02)
                p2_02_b = model02.intercept_
                #print("intersección: {:.4f}".format(p2_02))
                p1_02_b = model02.coef_
                #print("pendiente: {:.4f}. len(lpt02_grow)={:}".format(p1_02[0],len(lpt02_grow)))
                #print("largo: {}, pendiente: {}".format(largo,p1_02_b[0]))
                
                if p1_02_b[0] > p1_02_a[0]:
                    
                    p1_02 = p1_02_a[0]
                    p2_02 = p2_02_a
                    #le resto 10*50 porque significa q me sirve la pendiente anterior
                    #entonces tengo que volver los arrays al lago de ese momento
                    len_lpt02 = len(lpt02_grow) - (10*50)
                    lpt02_grow = lapsetime02[:len_lpt02]
                    print("len_lpt02: {}, pendiente: {}".format(len_lpt02,p1_02))
                    break
                
                elif p1_02_b[0] < p1_02_a[0]:
                    #print("pendiente: {:.4f}. len(len_lpt02)={:}".format(p1_02,len(len_lpt02)))
                    lpt02_grow = lapsetime02[:largo]
                    s_env02_grow = s_env02[:largo]
                    
                    p1_02_a[0] = p1_02_b[0] 

                
            cut_prediccion02 = prediccion02[:len_lpt02]
            cut_Alpt02 = Alpt02[:len_lpt02]
            cut_cut_env02 = cut_env02[:len_lpt02]
            cut_logAlpt02 = logAlpt02[:len_lpt02]

            plt.figure(2)
            plt.plot(lpt02_grow,cut_prediccion02, label="Alpt = {:.4f}t + {:.4f}".format(p1_02,p2_02))
            #plt.plot(lpt02_grow,cut_Alpt02, label="exponencial corregida por div geom")
            plt.plot(lpt02_grow,cut_logAlpt02, label="logaritmo de la amplitud")
            #plt.plot(lpt02_grow,cut_cut_env02, label="envolvente de la coda")
            plt.xlabel("tiempo [seg]")
            plt.ylabel("amplitud corregida por div geom")
            plt.title("LAPSETIME DETERMINADO CON LOOP")
            plt.legend()      
            
            
            
            # ----------------------------
            # 03
            print("Resultados del ajuste lineal. Frecuencias entre 1.5 - 2.5 Hz")
            
            largo = 150
            
            lpt03_grow = lapsetime03[:largo] ### acepto como mínimo 30 segundo a 50 mps
            s_env03_grow = s_env03[:largo]
            
            Alpt03 = s_env03_grow * lpt03_grow

            logAlpt03 = np.log(Alpt03).reshape(-1)
            lpt03_grow_reshape = lpt03_grow.reshape(-1,1)
            model03 = LinearRegression().fit(lpt03_grow_reshape, logAlpt03)
            prediccion03 = model03.predict(lpt03_grow_reshape)
            r_sq03 = model03.score(lpt03_grow_reshape, logAlpt03)
            p2_03_a = model03.intercept_
            p1_03_a = model03.coef_
            
            
            while len(lpt03_grow) < len(lapsetime03):

                largo = largo + 10*50
                lpt03_grow = lapsetime03[:largo]
                s_env03_grow = s_env03[:largo]                

                Alpt03 = s_env03_grow * lpt03_grow
                logAlpt03 = np.log(Alpt03).reshape(-1)
                lpt03_grow_reshape = lpt03_grow.reshape(-1,1)
                model03 = LinearRegression().fit(lpt03_grow_reshape, logAlpt03)
                prediccion03 = model03.predict(lpt03_grow_reshape)
                r_sq03 = model03.score(lpt03_grow_reshape, logAlpt03)
                p2_03_b = model03.intercept_
                #print("intersección: {:.4f}".format(p2_03))
                p1_03_b = model03.coef_
                #print("largo: {}, pendiente: {}".format(largo,p1_03_b[0]))
                
                if p1_03_b[0] > p1_03_a[0]:
                    
                    p1_03 = p1_03_a[0]
                    p2_03 = p2_03_a
                    len_lpt03 = len(lpt03_grow) - (10*50)
                    lpt03_grow = lapsetime03[:len_lpt03]
                    print("len_lpt03: {}, pendiente: {}".format(len_lpt03,p1_03))
                    break
                
                elif p1_03_b[0] < p1_03_a[0]:

                    lpt03_grow = lapsetime03[:largo]
                    s_env03_grow = s_env03[:largo]
                    
                    p1_03_a[0] = p1_03_b[0] 

                
            cut_prediccion03 = prediccion03[:len_lpt03]
            cut_Alpt03 = Alpt03[:len_lpt03]
            cut_cut_env03 = cut_env03[:len_lpt03]
            cut_logAlpt03 = logAlpt03[:len_lpt03]

            plt.figure(3)
            plt.plot(lpt03_grow,cut_prediccion03, label="Alpt = {:.4f}t + {:.4f}".format(p1_03,p2_03))
            #plt.plot(lpt03_grow,cut_Alpt03, label="exponencial corregida por div geom")
            plt.plot(lpt03_grow,cut_logAlpt03, label="logaritmo de la amplitud")
            #plt.plot(lpt03_grow,cut_cut_env03, label="envolvente de la coda")
            plt.xlabel("tiempo [seg]")
            plt.ylabel("amplitud corregida por div geom")
            plt.title("LAPSETIME DETERMINADO CON LOOP")
            plt.legend()   
            
            
            
            # ----------------------------
            # 04
            print("Resultados del ajuste lineal. Frecuencias entre 2.0 - 3.0 Hz")
            
            largo = 150
            
            lpt04_grow = lapsetime04[:largo] ### acepto como mínimo 30 segundo a 50 mps
            s_env04_grow = s_env04[:largo]
            
            Alpt04 = s_env04_grow * lpt04_grow

            logAlpt04 = np.log(Alpt04).reshape(-1)
            lpt04_grow_reshape = lpt04_grow.reshape(-1,1)
            model04 = LinearRegression().fit(lpt04_grow_reshape, logAlpt04)
            prediccion04 = model04.predict(lpt04_grow_reshape)
            r_sq04 = model04.score(lpt04_grow_reshape, logAlpt04)
            p2_04_a = model04.intercept_
            p1_04_a = model04.coef_
            
            
            while len(lpt04_grow) < len(lapsetime04):

                largo = largo + 10*50
                lpt04_grow = lapsetime04[:largo]
                s_env04_grow = s_env04[:largo]                

                Alpt04 = s_env04_grow * lpt04_grow
                logAlpt04 = np.log(Alpt04).reshape(-1)
                lpt04_grow_reshape = lpt04_grow.reshape(-1,1)
                model04 = LinearRegression().fit(lpt04_grow_reshape, logAlpt04)
                prediccion04 = model04.predict(lpt04_grow_reshape)
                r_sq04 = model04.score(lpt04_grow_reshape, logAlpt04)
                p2_04_b = model04.intercept_
                #print("intersección: {:.4f}".format(p2_04))
                p1_04_b = model04.coef_
                #print("largo: {}, pendiente: {}".format(largo,p1_04_b[0]))
                
                if p1_04_b[0] > p1_04_a[0]:
                    
                    p1_04 = p1_04_a[0]
                    p2_04 = p2_04_a
                    len_lpt04 = len(lpt04_grow) - (10*50)
                    lpt04_grow = lapsetime04[:len_lpt04]
                    print("len_lpt04: {}, pendiente: {}".format(len_lpt04,p1_04))
                    break
                
                elif p1_04_b[0] < p1_04_a[0]:

                    lpt04_grow = lapsetime04[:largo]
                    s_env04_grow = s_env04[:largo]
                    
                    p1_04_a[0] = p1_04_b[0] 

                
            cut_prediccion04 = prediccion04[:len_lpt04]
            cut_Alpt04 = Alpt04[:len_lpt04]
            cut_cut_env04 = cut_env04[:len_lpt04]
            cut_logAlpt04 = logAlpt04[:len_lpt04]

            plt.figure(4)
            plt.plot(lpt04_grow,cut_prediccion04, label="Alpt = {:.4f}t + {:.4f}".format(p1_04,p2_04))
            #plt.plot(lpt04_grow,cut_Alpt04, label="exponencial corregida por div geom")
            plt.plot(lpt04_grow,cut_logAlpt04, label="logaritmo de la amplitud")
            #plt.plot(lpt04_grow,cut_cut_env04, label="envolvente de la coda")
            plt.xlabel("tiempo [seg]")
            plt.ylabel("amplitud corregida por div geom")
            plt.title("LAPSETIME DETERMINADO CON LOOP")
            plt.legend()              
            
            
            
            # ----------------------------
            # 05
            print("Resultados del ajuste lineal. Frecuencias entre 2.5 - 3.5 Hz")
            
            largo = 150
            
            lpt05_grow = lapsetime05[:largo] ### acepto como mínimo 30 segundo a 50 mps
            s_env05_grow = s_env05[:largo]
            
            Alpt05 = s_env05_grow * lpt05_grow

            logAlpt05 = np.log(Alpt05).reshape(-1)
            lpt05_grow_reshape = lpt05_grow.reshape(-1,1)
            model05 = LinearRegression().fit(lpt05_grow_reshape, logAlpt05)
            prediccion05 = model05.predict(lpt05_grow_reshape)
            r_sq05 = model05.score(lpt05_grow_reshape, logAlpt05)
            p2_05_a = model05.intercept_
            p1_05_a = model05.coef_
            
            
            while len(lpt05_grow) < len(lapsetime05):

                largo = largo + 10*50
                lpt05_grow = lapsetime05[:largo]
                s_env05_grow = s_env05[:largo]                

                Alpt05 = s_env05_grow * lpt05_grow
                logAlpt05 = np.log(Alpt05).reshape(-1)
                lpt05_grow_reshape = lpt05_grow.reshape(-1,1)
                model05 = LinearRegression().fit(lpt05_grow_reshape, logAlpt05)
                prediccion05 = model05.predict(lpt05_grow_reshape)
                r_sq05 = model05.score(lpt05_grow_reshape, logAlpt05)
                p2_05_b = model05.intercept_
                #print("intersección: {:.4f}".format(p2_05))
                p1_05_b = model05.coef_
                #print("largo: {}, pendiente: {}".format(largo,p1_05_b[0]))
                
                if p1_05_b[0] > p1_05_a[0]:
                    
                    p1_05 = p1_05_a[0]
                    p2_05 = p2_05_a
                    len_lpt05 = len(lpt05_grow) - (10*50)
                    lpt05_grow = lapsetime05[:len_lpt05]
                    print("len_lpt05: {}, pendiente: {}".format(len_lpt05,p1_05))
                    break
                
                elif p1_05_b[0] < p1_05_a[0]:

                    lpt05_grow = lapsetime05[:largo]
                    s_env05_grow = s_env05[:largo]
                    
                    p1_05_a[0] = p1_05_b[0] 

                
            cut_prediccion05 = prediccion05[:len_lpt05]
            cut_Alpt05 = Alpt05[:len_lpt05]
            cut_cut_env05 = cut_env05[:len_lpt05]
            cut_logAlpt05 = logAlpt05[:len_lpt05]

            plt.figure(5)
            plt.plot(lpt05_grow,cut_prediccion05, label="Alpt = {:.4f}t + {:.4f}".format(p1_05,p2_05))
            #plt.plot(lpt05_grow,cut_Alpt05, label="exponencial corregida por div geom")
            plt.plot(lpt05_grow,cut_logAlpt05, label="logaritmo de la amplitud")
            #plt.plot(lpt05_grow,cut_cut_env05, label="envolvente de la coda")
            plt.xlabel("tiempo [seg]")
            plt.ylabel("amplitud corregida por div geom")
            plt.title("LAPSETIME DETERMINADO CON LOOP")
            plt.legend()              
            
            
            
            # ----------------------------
            # 06
            print("Resultados del ajuste lineal. Frecuencias entre 3.0 - 4.0 Hz")
            
            largo = 150
            
            lpt06_grow = lapsetime06[:largo] ### acepto como mínimo 30 segundo a 50 mps
            s_env06_grow = s_env06[:largo]
            
            Alpt06 = s_env06_grow * lpt06_grow

            logAlpt06 = np.log(Alpt06).reshape(-1)
            lpt06_grow_reshape = lpt06_grow.reshape(-1,1)
            model06 = LinearRegression().fit(lpt06_grow_reshape, logAlpt06)
            prediccion06 = model06.predict(lpt06_grow_reshape)
            r_sq06 = model06.score(lpt06_grow_reshape, logAlpt06)
            p2_06_a = model06.intercept_
            p1_06_a = model06.coef_
            
            
            while len(lpt06_grow) < len(lapsetime06):

                largo = largo + 10*50
                lpt06_grow = lapsetime06[:largo]
                s_env06_grow = s_env06[:largo]                

                Alpt06 = s_env06_grow * lpt06_grow
                logAlpt06 = np.log(Alpt06).reshape(-1)
                lpt06_grow_reshape = lpt06_grow.reshape(-1,1)
                model06 = LinearRegression().fit(lpt06_grow_reshape, logAlpt06)
                prediccion06 = model06.predict(lpt06_grow_reshape)
                r_sq06 = model06.score(lpt06_grow_reshape, logAlpt06)
                p2_06_b = model06.intercept_
                #print("intersección: {:.4f}".format(p2_06))
                p1_06_b = model06.coef_
                #print("largo: {}, pendiente: {}".format(largo,p1_06_b[0]))
                
                if p1_06_b[0] > p1_06_a[0]:
                    
                    p1_06 = p1_06_a[0]
                    p2_06 = p2_06_a
                    len_lpt06 = len(lpt06_grow) - (10*50)
                    lpt06_grow = lapsetime06[:len_lpt06]
                    print("len_lpt06: {}, pendiente: {}".format(len_lpt06,p1_06))
                    break
                
                elif p1_06_b[0] < p1_06_a[0]:

                    #largo = largo + 10*50
                    lpt06_grow = lapsetime06[:largo]
                    s_env06_grow = s_env06[:largo]
                    
                    p1_06_a[0] = p1_06_b[0] 

                
            cut_prediccion06 = prediccion06[:len_lpt06]
            cut_Alpt06 = Alpt06[:len_lpt06]
            cut_cut_env06 = cut_env06[:len_lpt06]
            cut_logAlpt06 = logAlpt06[:len_lpt06]

            plt.figure(6)
            plt.plot(lpt06_grow,cut_prediccion06, label="Alpt = {:.4f}t + {:.4f}".format(p1_06,p2_06))
            #plt.plot(lpt06_grow,cut_Alpt06, label="exponencial corregida por div geom")
            plt.plot(lpt06_grow,cut_logAlpt06, label="logaritmo de la amplitud")
            #plt.plot(lpt06_grow,cut_cut_env06, label="envolvente de la coda")
            plt.xlabel("tiempo [seg]")
            plt.ylabel("amplitud corregida por div geom")
            plt.title("LAPSETIME DETERMINADO CON LOOP")
            plt.legend()              
            
            
            
            # ----------------------------
            # 07
            print("Resultados del ajuste lineal. Frecuencias entre 3.5 - 4.5 Hz")
            
            largo = 150
            
            lpt07_grow = lapsetime07[:largo] ### acepto como mínimo 30 segundo a 50 mps
            s_env07_grow = s_env07[:largo]
            
            Alpt07 = s_env07_grow * lpt07_grow

            logAlpt07 = np.log(Alpt07).reshape(-1)
            lpt07_grow_reshape = lpt07_grow.reshape(-1,1)
            model07 = LinearRegression().fit(lpt07_grow_reshape, logAlpt07)
            prediccion07 = model07.predict(lpt07_grow_reshape)
            r_sq07 = model07.score(lpt07_grow_reshape, logAlpt07)
            p2_07_a = model07.intercept_
            p1_07_a = model07.coef_
            
            
            while len(lpt07_grow) < len(lapsetime07):

                largo = largo + 10*50
                lpt07_grow = lapsetime07[:largo]
                s_env07_grow = s_env07[:largo]                

                Alpt07 = s_env07_grow * lpt07_grow
                logAlpt07 = np.log(Alpt07).reshape(-1)
                lpt07_grow_reshape = lpt07_grow.reshape(-1,1)
                model07 = LinearRegression().fit(lpt07_grow_reshape, logAlpt07)
                prediccion07 = model07.predict(lpt07_grow_reshape)
                r_sq07 = model07.score(lpt07_grow_reshape, logAlpt07)
                p2_07_b = model07.intercept_
                #print("intersección: {:.4f}".format(p2_07))
                p1_07_b = model07.coef_
                #print("largo: {}, pendiente: {}".format(largo,p1_07_b[0]))
                
                if p1_07_b[0] > p1_07_a[0]:
                    
                    p1_07 = p1_07_a[0]
                    p2_07 = p2_07_a
                    len_lpt07 = len(lpt07_grow) - (10*50)
                    lpt07_grow = lapsetime07[:len_lpt07]
                    print("len_lpt07: {}, pendiente: {}".format(len_lpt07,p1_07))
                    break
                
                elif p1_07_b[0] < p1_07_a[0]:

                    #largo = largo + 10*50
                    lpt07_grow = lapsetime07[:largo]
                    s_env07_grow = s_env07[:largo]
                    
                    p1_07_a[0] = p1_07_b[0] 

                
            cut_prediccion07 = prediccion07[:len_lpt07]
            cut_Alpt07 = Alpt07[:len_lpt07]
            cut_cut_env07 = cut_env07[:len_lpt07]
            cut_logAlpt07 = logAlpt07[:len_lpt07]

            plt.figure(7)
            plt.plot(lpt07_grow,cut_prediccion07, label="Alpt = {:.4f}t + {:.4f}".format(p1_07,p2_07))
            #plt.plot(lpt07_grow,cut_Alpt07, label="exponencial corregida por div geom")
            plt.plot(lpt07_grow,cut_logAlpt07, label="logaritmo de la amplitud")
            #plt.plot(lpt07_grow,cut_cut_env07, label="envolvente de la coda")
            plt.xlabel("tiempo [seg]")
            plt.ylabel("amplitud corregida por div geom")
            plt.title("LAPSETIME DETERMINADO CON LOOP")
            plt.legend()              
            
            
            
            # ----------------------------
            # 08
            print("Resultados del ajuste lineal. Frecuencias entre 4.0 - 5.0 Hz")
            
            largo = 150
            
            lpt08_grow = lapsetime08[:largo] ### acepto como mínimo 30 segundo a 50 mps
            s_env08_grow = s_env08[:largo]
            
            Alpt08 = s_env08_grow * lpt08_grow

            logAlpt08 = np.log(Alpt08).reshape(-1)
            lpt08_grow_reshape = lpt08_grow.reshape(-1,1)
            model08 = LinearRegression().fit(lpt08_grow_reshape, logAlpt08)
            prediccion08 = model08.predict(lpt08_grow_reshape)
            r_sq08 = model08.score(lpt08_grow_reshape, logAlpt08)
            p2_08_a = model08.intercept_
            p1_08_a = model08.coef_
            
            
            while len(lpt08_grow) < len(lapsetime08):

                largo = largo + 10*50
                lpt08_grow = lapsetime08[:largo]
                s_env08_grow = s_env08[:largo]                

                Alpt08 = s_env08_grow * lpt08_grow
                logAlpt08 = np.log(Alpt08).reshape(-1)
                lpt08_grow_reshape = lpt08_grow.reshape(-1,1)
                model08 = LinearRegression().fit(lpt08_grow_reshape, logAlpt08)
                prediccion08 = model08.predict(lpt08_grow_reshape)
                r_sq08 = model08.score(lpt08_grow_reshape, logAlpt08)
                p2_08_b = model08.intercept_
                #print("intersección: {:.4f}".format(p2_08))
                p1_08_b = model08.coef_
                #print("largo: {}, pendiente: {}".format(largo,p1_08_b[0]))
                
                if p1_08_b[0] > p1_08_a[0]:
                    
                    p1_08 = p1_08_a[0]
                    p2_08 = p2_08_a
                    len_lpt08 = len(lpt08_grow) - (10*50)
                    lpt08_grow = lapsetime08[:len_lpt08]
                    print("len_lpt07: {}, pendiente: {}".format(len_lpt08,p1_08))
                    break
                
                elif p1_08_b[0] < p1_08_a[0]:

                    #largo = largo + 10*50
                    lpt08_grow = lapsetime08[:largo]
                    s_env08_grow = s_env08[:largo]
                    
                    p1_08_a[0] = p1_08_b[0] 

                
            cut_prediccion08 = prediccion08[:len_lpt08]
            cut_Alpt08 = Alpt08[:len_lpt08]
            cut_cut_env08 = cut_env08[:len_lpt08]
            cut_logAlpt08 = logAlpt08[:len_lpt08]

            plt.figure(8)
            plt.plot(lpt08_grow,cut_prediccion08, label="Alpt = {:.4f}t + {:.4f}".format(p1_08,p2_08))
            #plt.plot(lpt08_grow,cut_Alpt08, label="exponencial corregida por div geom")
            plt.plot(lpt08_grow,cut_logAlpt08, label="logaritmo de la amplitud")
            #plt.plot(lpt08_grow,cut_cut_env08, label="envolvente de la coda")
            plt.xlabel("tiempo [seg]")
            plt.ylabel("amplitud corregida por div geom")
            plt.title("LAPSETIME DETERMINADO CON LOOP")
            plt.legend()              
            
            
#%%
#-----------------------------------------------------------------------------
# 18. Segundo ajuste: exponencial y cálculo del Q
#-----------------------------------------------------------------------------
# Se realiza el ajuste exponencial dando como valores de inicio los encontrados en
# el ajuste lineal

            def model_exp(t, a, b):
                return a * np.exp(b*t)

            # ----------------------------
            # 01
            print("--------------------------------")
            print("Frecuencias entre 0.5 - 1.5 Hz")

            if p1_01 < 0:
                # p0 tiene que ser una tupla
                p0_01 = (np.exp(p2_01),p1_01) # hago p1[0] xq es un array y yo necesito el elemento
                # ajuste
                popt_01, pcov_01 = curve_fit(model_exp, lpt01_grow, cut_Alpt01, p0_01)
                Q_01 = -np.pi*1/popt_01[1]
                print("Qc = {}".format(Q_01))
            
            elif p1_01 > 0:
                Q_01 = 0
                print("No se puede calcular Qc")
            
            
            
            # ----------------------------
            # 02
            print("--------------------------------")
            print("Frecuencias entre 1.0 - 2.0 Hz")

            if p1_02 < 0:
                p0_02 = (np.exp(p2_02),p1_02)
                popt_02, pcov_02 = curve_fit(model_exp, lpt02_grow, cut_Alpt02, p0_02)
                Q_02 = -np.pi*1/popt_02[1]
                print("Qc = {:.4f}".format(Q_02))
            
            elif p1_02 > 0:
                Q_02 = 0
                print("No se puede calcular Qc")

            
            # ----------------------------
            # 03
            print("--------------------------------")
            print("Frecuencias entre 1.5 - 2.5 Hz")
            
            if p1_03 < 0:
                p0_03 = (np.exp(p2_03),p1_03)
                popt_03, pcov_03 = curve_fit(model_exp, lpt03_grow, cut_Alpt03, p0_03)
                Q_03 = -np.pi*1/popt_03[1]
                print("Qc = {:.4f}".format(Q_03))
            
            elif p1_03 > 0:
                Q_03 = 0
                print("No se puede calcular Qc")


            # ----------------------------
            # 04
            print("--------------------------------")
            print("Frecuencias entre 2.0 - 3.0 Hz")
            
            if p1_04 < 0:
                p0_04 = (np.exp(p2_04),p1_04)
                popt_04, pcov_04 = curve_fit(model_exp, lpt04_grow, cut_Alpt04, p0_04)
                Q_04 = -np.pi*1/popt_04[1]
                print("Qc= {:.4f}".format(Q_04))
            
            elif p1_04 > 0:
                Q_04 = 0
                print("No se puede calcular Qc")


            # ----------------------------
            # 05
            print("--------------------------------")
            print("Frecuencias entre 2.5 - 3.5 Hz")
            
            if p1_05 < 0:
                p0_05 = (np.exp(p2_05),p1_05)
                popt_05, pcov_05 = curve_fit(model_exp, lpt05_grow, cut_Alpt05, p0_05)
                Q_05 = -np.pi*1/popt_05[1]
                print("Qc = {:.4f}".format(Q_05))
            
            elif p1_05 > 0:
                Q_05 = 0
                print("No se puede calcular Qc")


            # ----------------------------
            # 06
            print("--------------------------------")
            print("Frecuencias entre 3.0 - 4.0 Hz")
            
            if p1_06 < 0:
                p0_06 = (np.exp(p2_06),p1_06)
                popt_06, pcov_06 = curve_fit(model_exp, lpt06_grow, cut_Alpt06, p0_06)
                Q_06 = -np.pi*1/popt_06[1]
                print("Qc = {:.4f}".format(Q_06))
            
            elif p1_06 > 0:
                Q_06 = 0
                print("No se puede calcular Qc")
            
            
            # ----------------------------
            # 07
            print("--------------------------------")
            print("Frecuencias entre 3.5 - 4.5 Hz")
            
            if p1_07 < 0:
                p0_07 = (np.exp(p2_07),p1_07)
                popt_07, pcov_07 = curve_fit(model_exp, lpt07_grow, cut_Alpt07, p0_07)
                Q_07 = -np.pi*1/popt_07[1]
                print("Qc = {:.4f}".format(Q_07))
            
            elif p1_07 > 0:
                Q_07 = 0
                print("No se puede calcular Qc")

            # ----------------------------
            # 08
            print("--------------------------------")
            print("Frecuencias entre 4.0 - 5.0 Hz")
            
            if p1_08 < 0:
                p0_08 = (np.exp(p2_08),p1_08)
                popt_08, pcov_08 = curve_fit(model_exp, lpt08_grow, cut_Alpt08, p0_08)
                Q_08 = -np.pi*1/popt_08[1]
                print("Qc = {:.4f}".format(Q_08))
            
            elif p1_08 > 0:
                Q_08 = 0
                print("No se puede calcular Qc")
           

#%%
#-----------------------------------------------------------------------------
# 19. Difusión
#-----------------------------------------------------------------------------
            def model_dif(t, a1,a2,a3):
                return  a1-a2*t-a3/t
            
            dist = r_fg16
            p = 2 #div geom O.S p=2, O.I p=3


            # ----------------------------
            # 01
            print("--------------------------------")
            print("Frecuencias entre 0.5 - 1.5 Hz")
            
            if p1_01 < 0:
                p0_01_d = (np.exp(p2_01),-p1_01,0)
                popt_01_d, pcov_01_d = curve_fit(model_dif, lpt01_grow, np.log(cut_Alpt01), p0_01_d)
                
                Qi_01 = 2*np.pi*(w1+w2)/(2*popt_01_d[1])
                d_01 = dist**2/(4*popt_01_d[2])
                Qs_01 = np.pi*(w1+w2)*p*d_01/vel**2
                print("Qi={:.4f} y Qs={:.4f}".format(Qi_01,Qs_01))
            
            elif p1_01 > 0:
                Qi_01 = 0
                Qs_01 = 0
                print("No se puede calcular Qi ni Qs")


            # ----------------------------
            # 02
            print("--------------------------------")
            print("Frecuencias entre 1.0 - 2.0 Hz")
            
            if p1_02 < 0:
                p0_02_d = (np.exp(p2_02),-p1_02,0)
                popt_02_d, pcov_02_d = curve_fit(model_dif, lpt02_grow, np.log(cut_Alpt02), p0_02_d)
            
                Qi_02 = 2*np.pi*(w3+w4)/(2*popt_02_d[1])
                d_02 = dist**2/(4*popt_02_d[2])
                Qs_02 = np.pi*(w3+w4)*p*d_02/vel**2
                print("Qi={:.4f} y Qs={:.4f}".format(Qi_02,Qs_02))
            
            elif p1_02 > 0:
                Qi_02 = 0
                Qs_02 = 0
                print("No se puede calcular Qi ni Qs")


            # ----------------------------
            # 03
            print("--------------------------------")
            print("Frecuencias entre 1.5 - 2.5 Hz")
            
            if p1_03 < 0:
                p0_03_d = (np.exp(p2_03),-p1_03,0)
                popt_03_d, pcov_03_d = curve_fit(model_dif, lpt03_grow, np.log(cut_Alpt03), p0_03_d)
                
                Qi_03 = 2*np.pi*(w5+w6)/(2*popt_03_d[1])
                d_03 = dist**2/(4*popt_03_d[2])
                Qs_03 = np.pi*(w5+w6)*p*d_03/vel**2
                print("Qi={:.4f} y Qs={:.4f}".format(Qi_03,Qs_03))
            
            elif p1_03 > 0:
                Qi_03 = 0
                Qs_03 = 0
                print("No se puede calcular Qi ni Qs")


            # ----------------------------
            # 04
            print("--------------------------------")
            print("Frecuencias entre 2.0 - 3.0 Hz")
            
            if p1_04 < 0:
                p0_04_d = (np.exp(p2_04),-p1_04,0)
                popt_04_d, pcov_04_d = curve_fit(model_dif, lpt04_grow, np.log(cut_Alpt04), p0_04_d)
                
                Qi_04 = 2*np.pi*(w7+w8)/(2*popt_04_d[1])
                d_04 = dist**2/(4*popt_04_d[2])
                Qs_04 = np.pi*(w7+w8)*p*d_04/vel**2
                print("Qi={:.4f} y Qs={:.4f}".format(Qi_04,Qs_04))
            
            elif p1_04 > 0:
                Qi_04 = 0
                Qs_04 = 0
                print("No se puede calcular Qi ni Qs")


            # ----------------------------
            # 05
            print("--------------------------------")
            print("Frecuencias entre 2.5 - 3.5 Hz")
            
            if p1_05 < 0:
                p0_05_d = (np.exp(p2_05),-p1_05,0)
                popt_05_d, pcov_05_d = curve_fit(model_dif, lpt05_grow, np.log(cut_Alpt05), p0_05_d)
                
                Qi_05 = 2*np.pi*(w9+w10)/(2*popt_05_d[1])
                d_05 = dist**2/(4*popt_05_d[2])
                Qs_05 = np.pi*(w9+w10)*p*d_05/vel**2
                print("Qi={:.4f} y Qs={:.4f}".format(Qi_05,Qs_05))
            
            elif p1_05 > 0:
                Qi_05 = 0
                Qs_05 = 0
                print("No se puede calcular Qi ni Qs")

            
            # ----------------------------
            # 06
            print("--------------------------------")
            print("Frecuencias entre 3.0 - 4.0 Hz")

            if p1_06 < 0:
                p0_06_d = (np.exp(p2_06),-p1_06,0)
                popt_06_d, pcov_06_d = curve_fit(model_dif, lpt06_grow, np.log(cut_Alpt06), p0_06_d)
                
                Qi_06 = 2*np.pi*(w11+w12)/(2*popt_06_d[1])
                d_06 = dist**2/(4*popt_06_d[2])
                Qs_06 = np.pi*(w11+w12)*p*d_06/vel**2
                print("Qi={:.4f} y Qs={:.4f}".format(Qi_06,Qs_06))
            
            elif p1_06 > 0:
                Qi_06 = 0
                Qs_06 = 0
                print("No se puede calcular Qi ni Qs")


            # ----------------------------
            # 07
            print("--------------------------------")
            print("Frecuencias entre 3.5 - 4.5 Hz")
            
            if p1_07 < 0:
                p0_07_d = (np.exp(p2_07),-p1_07,0)
                popt_07_d, pcov_07_d = curve_fit(model_dif, lpt07_grow, np.log(cut_Alpt07), p0_07_d)
                
                Qi_07 = 2*np.pi*(w13+w14)/(2*popt_07_d[1])
                d_07 = dist**2/(4*popt_07_d[2])
                Qs_07 = np.pi*(w13+w14)*p*d_07/vel**2
                print("Qi={:.4f} y Qs={:.4f}".format(Qi_07,Qs_07))
            
            elif p1_07 > 0:
                Qi_07 = 0
                Qs_07 = 0
                print("No se puede calcular Qi ni Qs")


            # ----------------------------
            # 08
            print("--------------------------------")
            print("Frecuencias entre 4.0 - 5.0 Hz")
            
            if p1_08 < 0:
                p0_08_d = (np.exp(p2_08),-p1_08,0)
                popt_08_d, pcov_08_d = curve_fit(model_dif, lpt08_grow, np.log(cut_Alpt08), p0_08_d)
                
                Qi_08 = 2*np.pi*(w15+w16)/(2*popt_08_d[1])
                d_08 = dist**2/(4*popt_08_d[2])
                Qs_08 = np.pi*(w15+w16)*p*d_08/vel**2
                print("Qi={:.4f} y Qs={:.4f}".format(Qi_08,Qs_08))
            
            elif p1_08 > 0:
                Qi_08 = 0
                Qs_08 = 0
                print("No se puede calcular Qi ni Qs")

#%%



