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



#-----------------------------------------------------------------------------
# 1. Se abre el archivo de tags y se busca si ya existe el archivo de salida
#-----------------------------------------------------------------------------

### path del file donde se encuentran los tags
#tags_dir='/media/gaby/Backup Plus/SISMOLOGIA_DATOS_Y_DOCS/Vn_FUEGO/tags'
tags_dir='C:/Users/Propietario/Desktop/tesis_de_grado/ATENUACION/DATOS/FG16'

### tags_dir + 'nombre del archivo que contiene los tags'
#tags_file=tags_dir +'/fecha_tipo_dur_ampl_max_LP_2020_MEJORADO_1.csv'
tags_file=tags_dir + '/buenos_y_malos_fecha_tipo_dur_ampl_max_LP_2020_MEJORADO_1.csv'

#data_dir='/media/gaby/Backup Plus/SISMOLOGIA_DATOS_Y_DOCS/Vn_FUEGO/FG16'
data_dir='C:/Users/Propietario/Desktop/tesis_de_grado/ATENUACION/DATOS/FG16'

### los eventos a leer son de tipo LP. El tag de un evento comienza con LP y
# finaliza con Noise
EventCode='LP_B'
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
            
            if tipo_inicio == 'LP_B' and tipo_fin == 'Noise': # confirmo que tiene inicio y fin
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

#-----------------------------------------------------------------------------
# Calculo el espectro de la señal con taper y sin taper y comparo
#-----------------------------------------------------------------------------
# Taper del 10%
            ns = len(tr_taper10)
            f = rfftfreq(ns, dt)
            fft_tr_taper10 = rfft(tr_taper10)   # fft de la señal con taper
           
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# 01 FILTRO PASA BANDA - Fc = 1Hz - TODOS LOS CÁLCULOS
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------

### FIltro Pasa Banda a la señal           
            # 01
            fm = 1/dt                   # frec de muestreo
            #w1 = 0.5/fm
            #w2 = 1.5/fm
            w1=0.5
            w2=1.5
            tr_taper10_copy01 = tr_taper10.copy()   # hago una copia para no sobreescribir
            evento_filtrado_01 = tr_taper10_copy01.filter("bandpass", freqmin=w1, freqmax=w2, corners=2, zerophase=False)


#------------------------------------------------
### Calculo FFT del evento previamente trabajado

            f = rfftfreq(ns, dt)
            fft_01 = rfft(evento_filtrado_01)   # fft de BP e/ 0.5 - 1.5 Hz
      
        
#------------------------------------------------
### Cálculo de envolvente

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
            

            dom_time_01 = irfft(fft_01,n=len(evento[0]))     # cálculo de la transformada inversa (parte real)
            
            # cambio el tipo de variable
            st_dom_time_01 = (new_st(dom_time_01)).slice(inicio,fin)
            t_01 = vec_tiempo(st_dom_time_01)
            
            # cálculo de la envolvente
            envolvente_01 = obspy.signal.filter.envelope(st_dom_time_01[0].data)
            
            
#------------------------------------------------
### Valor máximo de cada envolvente y corte de la misma

            def posicion_max(envolvente_num):
                val = np.amax(envolvente_num)
                lista = envolvente_num.tolist()
                posicion = lista.index(val)
                return posicion
            
            max_01 = posicion_max(envolvente_01)            
            
# corto la envolvente y me quedo con la parte entre el valor máximo
# y la última muestra antes del post_evento            
            cut_env01 = envolvente_01[max_01:]

# ídem anterior, pero con el vector tiempo
            cut_t01 = t_01[max_01:]
            
            
#------------------------------------------------            
### Cálculo del tiempo correspondiente al máx valor de la envolvente

# me devuelve el número de la siguiente forma:
# datetime.datetime(año, mes, día, hora, min, seg, VER)
            t_max01 = inicio_pre + t_01[max_01]
            

#------------------------------------------------    
### Cálculo de distancia en km entre cráter y estación 

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
            
                       
#------------------------------------------------               
### Cálculo del tiempo que le lleva a la señal ir desde la fuente hasta la
# estación FG16

            vel = 3.5  #km/seg
            t_evento = r_fg16/vel #desde cráter a fg16
            
                       
#------------------------------------------------               
### Cálculo de vectores lapsetime

# Se calcula el tiempo real de cómo se dio el evento. Ya que cuando se 
# extrae el evento para poder procesarlo se le asocia un vector tiempo de dura-
# ción genérico, no teniendo en cuenta el tiempo que tardó el LP en viajar
# de la fuente a la estación, etc.
            delta_max01 = t_max01 - inicio
            lapsetime01 = t_evento + delta_max01 + cut_t01


###### CÁLCULOS PRIORIZANDO BUENOS AJUSTES:
#------------------------------------------------               
### Primer ajuste: lineal
#-----------------------------------------------------------------------------       
# Se hace un ajuste lineal, los datos de entrada serían el tiempo real (no el
# del registo) y el log de las amplitudes.
# El objetivo de este paso es buscar valores para "dárselos" al ajuste
# exponencial de la envolvente, que se hará después, así tiene un rango de 
# valores por donde comenzar a ajustar.

            print(inicio)
            # ----------------------------
            # 01
            # se hace el ajute lineal a los datos filtrados por BP 0.5 - 1.5 Hz            
            print("Resultados del ajuste lineal. Frecuencias entre 0.5 - 1.5 Hz")
            
            largo = 1500
            
            lpt01_grow = lapsetime01[:largo] ### acepto como mínimo 30 segundo a 50 mps
            cut_env01_grow = cut_env01[:largo]
            
            Alpt01 = cut_env01_grow * lpt01_grow
            # se hace reshape de los vectores para que cumplan la condición que
            # pide el paquete de regresión lineal
            logAlpt01 = np.log(Alpt01).reshape(-1)
            lpt01_grow_reshape = lpt01_grow.reshape(-1,1)
            model01 = LinearRegression().fit(lpt01_grow_reshape, logAlpt01)
            prediccion01 = model01.predict(lpt01_grow_reshape)
            r_sq01 = model01.score(lpt01_grow_reshape, logAlpt01)
            p2_01_a = model01.intercept_
            p1_01_a = model01.coef_

            
            while len(lpt01_grow) < len(lapsetime01):

                largo = largo + 10*50
                lpt01_grow = lapsetime01[:largo]
                cut_env01_grow = cut_env01[:largo]                

                Alpt01 = cut_env01_grow * lpt01_grow
                logAlpt01 = np.log(Alpt01).reshape(-1)
                lpt01_grow_reshape = lpt01_grow.reshape(-1,1)
                model01 = LinearRegression().fit(lpt01_grow_reshape, logAlpt01)
                prediccion01 = model01.predict(lpt01_grow_reshape)
                r_sq01 = model01.score(lpt01_grow_reshape, logAlpt01)
                p2_01_b = model01.intercept_
                p1_01_b = model01.coef_
        
                
                if p1_01_b[0] > p1_01_a[0]: #>                    
                    p1_01 = p1_01_a[0]
                    p2_01 = p2_01_a
                    #le resto 10*50 porque significa q me sirve la pendiente anterior
                    #entonces tengo que volver los arrays al lago de ese momento
                    #len_lpt01 = len(lpt01_grow) - (10*50)
                    len_lpt01 = len(lapsetime01[:(largo-(50*10))])
                    lpt01_grow = lapsetime01[:len_lpt01]
                    print("len_lpt01: {}, pendiente: {:.4f}".format(len_lpt01,p1_01))
                    break
                
                elif p1_01_b[0] < p1_01_a[0]: #<
                    lpt01_grow = lapsetime01[:largo]
                    cut_env01_grow = cut_env01[:largo]                    
                    p1_01_a[0] = p1_01_b[0] 
            
            if len(lpt01_grow) == len(lapsetime01):
                print("estoy aca")
                p1_01 = p1_01_a[0]
                p2_01 = p2_01_a
                len_lpt01 = len(lpt01_grow)   
                p1_01_b = [0] # variable creada para que no salte error en el siguiente if
                p1_01_b = np.array(p1_01_b)
                print("len_lpt01: {}, pendiente: {:.4f}".format(len_lpt01,p1_01))

            if p1_01_b[0] == p1_01_a[0]:
                p1_01 = p1_01_a[0]
                p2_01 = p2_01_a
                len_lpt01 = len(lpt01_grow)   
                print("len_lpt01: {}, pendiente: {:.4f}".format(len_lpt01,p1_01))
    
            cut_prediccion01 = prediccion01[:len_lpt01]
            cut_Alpt01 = Alpt01[:len_lpt01]
            cut_cut_env01 = cut_env01[:len_lpt01]
            cut_logAlpt01 = logAlpt01[:len_lpt01]


#------------------------------------------------ 
### Segundo ajuste: exponencial y cálculo del Q

# Se realiza el ajuste exponencial dando como valores de inicio los encontrados en
# el ajuste lineal

            def model_exp(t, a, b):
                return a * np.exp(b*t)

            print(inicio)
            print("--------------------------------")
            print("Frecuencias entre 0.5 - 1.5 Hz")

            if p1_01 < 0:
                # p0 tiene que ser una tupla
                p0_01 = (np.exp(p2_01),p1_01) # hago p1[0] xq es un array y yo necesito el elemento
                # ajuste
                popt_01, pcov_01 = curve_fit(model_exp, lpt01_grow, cut_Alpt01, p0_01,maxfev=15000)
                Q_01 = -np.pi*1/popt_01[1]
                print("Qc = {}".format(Q_01))
            
            elif p1_01 > 0:
                Q_01 = 0
                popt_01 = np.zeros(shape=2)
                print("No se puede calcular Qc")

#------------------------------------------------             
### Difusión
### gráfico difusión

            def ajuste_dif(a1,a2,a3,t):
                return np.exp(a1 + a2*t + a3*(1/t))


            def model_dif(t, a1,a2,a3):
                return  a1 + a2*t + a3/t
            
            dist = r_fg16
            p = 2 #div geom O.S p=2, O.I p=3

            print(inicio)
            print("--------------------------------")
            print("Frecuencias entre 0.5 - 1.5 Hz")
            
            if p1_01 < 0:
                lapsetime01_entero = t_evento + delta_max01 + t_01 # lapsetime desde el inicio del tag
                
                len_diferencia_lpt_01 = len(lapsetime01_entero)-len(lapsetime01) # veo cuántas muestas hay entre el lapsetime desde max amplt y el lapsetime desde el inicio
                ## la difetencia anterior está correcta xq lo verifique mirando grafico plt.plot(t_02,envolvente_02)
                
                cortar_lpt_final_01 = len_lpt01 + len_diferencia_lpt_01 # porque len_lpt esta considerado desde amplt max hasta el mejor ajuste
                lapsetime01_entero = lapsetime01_entero[:cortar_lpt_final_01]
                Alpt01_entero = envolvente_01[:cortar_lpt_final_01] * lapsetime01_entero
                
                popt_01_d, pcov_01_d = curve_fit(model_dif, lapsetime01_entero, np.log(Alpt01_entero))

                f01 = (w1+w2)/2 # frec central en Hz
                Qi_01 = (-2)*np.pi*f01/(popt_01_d[1]) # temporal
                
                d_01 = (-1)*dist**2/(4*popt_01_d[2]) # coef de difusividad

                Qs_01 = (2)*np.pi*f01*p*d_01/vel**2 # espacial
                
                print("Qi={:.4f} y Qs={:.4f}".format(Qi_01,Qs_01))
                
                
            elif p1_01 > 0:
                Qi_01 = 0
                Qs_01 = 0
                f01 = (w1+w2)/2
                popt_01_d = np.zeros(shape=3)
                d_01 = 0
                print("No se puede calcular Qi ni Qs")
            

#------------------------------------------------  
### Grabado de datos
            #Fc = f01
            field_names = ['UTCDateTime','Fc','lapsetime','dist','p_dg','vel','a_lin','b_lin','a_exp','b_exp','Qc_SBS','a1_dif','a2_dif','a3_dif','Qi','d','Qs']

            #lpt=lapsetime01[-1:] # asi en cada lapse time, nos quedamos con el ultimo elemento como float para imprimirlo
            #lptfin=lpt[0]
            dict = {'UTCDateTime':inicio,'Fc':f01,'lapsetime':len_lpt01,'dist':dist,'p_dg':p,'vel':vel,'a_lin':np.exp(p2_01),'b_lin':-p1_01,'a_exp':popt_01[0],'b_exp':popt_01[1],'Qc_SBS': Q_01,'a1_dif':popt_01_d[0],'a2_dif':popt_01_d[1],'a3_dif':popt_01_d[2],'Qi':Qi_01,'d':d_01,'Qs':Qs_01}
            
            with open(path_archivo_salida, 'a') as csv_file:
                dict_object = csv.DictWriter(csv_file, fieldnames=field_names) 
              
                dict_object.writerow(dict)

                

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# 02 FILTRO PASA BANDA - Fc = 1.5Hz - TODOS LOS CÁLCULOS
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
                
### FIltro Pasa Banda a la señal           
            w3 = 1.0
            w4 = 2.0
            tr_taper10_copy02 = tr_taper10.copy()   # hago una copia para no sobreescribir
            evento_filtrado_02 = tr_taper10_copy02.filter("bandpass", freqmin=w3, freqmax=w4, corners=2, zerophase=False)


#------------------------------------------------
### Calculo FFT del evento previamente trabajado

            fft_02 = rfft(evento_filtrado_02)   # fft de BP e/ 0.5 - 1.5 Hz
      
              
#------------------------------------------------
### Cálculo de envolvente

            dom_time_02 = irfft(fft_02,n=len(evento[0]))     # cálculo de la transformada inversa (parte real)
            st_dom_time_02 = (new_st(dom_time_02)).slice(inicio,fin)
            t_02 = vec_tiempo(st_dom_time_02) # vector tiempo
            envolvente_02 = obspy.signal.filter.envelope(st_dom_time_02[0].data)
            
                        
#------------------------------------------------
### Valor máximo de cada envolvente y corte de la misma
            
            max_02 = posicion_max(envolvente_02)            
            
# corto la envolvente y me quedo con la parte entre el valor máximo
# y la última muestra antes del post_evento            
            cut_env02 = envolvente_02[max_02:]

# ídem anterior, pero con el vector tiempo
            cut_t02 = t_02[max_02:]
            
                       
#------------------------------------------------            
### Cálculo del tiempo correspondiente al máx valor de la envolvente

# me devuelve el número de la siguiente forma:
# datetime.datetime(año, mes, día, hora, min, seg, VER)
            t_max02 = inicio_pre + t_02[max_02]
            

#------------------------------------------------               
### Cálculo de vectores lapsetime

# Se calcula el tiempo real de cómo se dio el evento. Ya que cuando se 
# extrae el evento para poder procesarlo se le asocia un vector tiempo de dura-
# ción genérico, no teniendo en cuenta el tiempo que tardó el LP en viajar
# de la fuente a la estación, etc.
            delta_max02 = t_max02 - inicio
            lapsetime02 = t_evento + delta_max02 + cut_t02


###### CÁLCULOS PRIORIZANDO BUENOS AJUSTES:
#------------------------------------------------            
### Primer ajuste: lineal
       
# Se hace un ajuste lineal, los datos de entrada serían el tiempo real (no el
# del registo) y el log de las amplitudes.
# El objetivo de este paso es buscar valores para "dárselos" al ajuste
# exponencial de la envolvente, que se hará después, así tiene un rango de 
# valores por donde comenzar a ajustar.

            print("Resultados del ajuste lineal. Frecuencias entre 1.0 - 2.0 Hz")
            
            largo = 1500
            
            lpt02_grow = lapsetime02[:largo] ### acepto como mínimo 30 segundo a 50 mps
            cut_env02_grow = cut_env02[:largo]
            
            Alpt02 = cut_env02_grow * lpt02_grow

            logAlpt02 = np.log(Alpt02).reshape(-1)
            lpt02_grow_reshape = lpt02_grow.reshape(-1,1)
            model02 = LinearRegression().fit(lpt02_grow_reshape, logAlpt02)
            prediccion02 = model02.predict(lpt02_grow_reshape)
            r_sq02 = model02.score(lpt02_grow_reshape, logAlpt02)
            p2_02_a = model02.intercept_
            p1_02_a = model02.coef_
                        
            while len(lpt02_grow) < len(lapsetime02):
                
                largo = largo + 10*50
                lpt02_grow = lapsetime02[:largo]
                cut_env02_grow = cut_env02[:largo]                
                    
                Alpt02 = cut_env02_grow * lpt02_grow
                logAlpt02 = np.log(Alpt02).reshape(-1)
                lpt02_grow_reshape = lpt02_grow.reshape(-1,1)
                model02 = LinearRegression().fit(lpt02_grow_reshape, logAlpt02)
                prediccion02 = model02.predict(lpt02_grow_reshape)
                r_sq02 = model02.score(lpt02_grow_reshape, logAlpt02)
                p2_02_b = model02.intercept_
                p1_02_b = model02.coef_

                if p1_02_b[0] > p1_02_a[0]:
                    
                    p1_02 = p1_02_a[0]
                    p2_02 = p2_02_a
                    #le resto 10*50 porque significa q me sirve la pendiente anterior
                    #entonces tengo que volver los arrays al lago de ese momento
                    len_lpt02 = len(lapsetime02[:(largo-(50*10))])
                    lpt02_grow = lapsetime02[:len_lpt02]
                    print("len_lpt02: {}, pendiente: {:.4f}".format(len_lpt02,p1_02))
                    break
                
                elif p1_02_b[0] < p1_02_a[0]:
                    lpt02_grow = lapsetime02[:largo]
                    cut_env02_grow = cut_env02[:largo]
                    p1_02_a[0] = p1_02_b[0] 
                    
            if len(lpt02_grow) == len(lapsetime02):
                print("estoy aca")
                p1_02 = p1_02_a[0]
                p2_02 = p2_02_a
                len_lpt02 = len(lpt02_grow)   
                p1_02_b = [0] # variable creada para que no salte error en el siguiente if
                p1_02_b = np.array(p1_02_b)                
                print("len_lpt02: {}, pendiente: {:.4f}".format(len_lpt02,p1_02))

            if p1_02_b[0] == p1_02_a[0]:
                p1_02 = p1_02_a[0]
                p2_02 = p2_02_a
                len_lpt02 = len(lpt02_grow)   
                print("len_lpt02: {}, pendiente: {:.4f}".format(len_lpt02,p1_02))

            cut_prediccion02 = prediccion02[:len_lpt02]
            cut_Alpt02 = Alpt02[:len_lpt02]
            cut_cut_env02 = cut_env02[:len_lpt02]
            cut_logAlpt02 = logAlpt02[:len_lpt02]


#------------------------------------------------ 
### Segundo ajuste: exponencial y cálculo del Q

# Se realiza el ajuste exponencial dando como valores de inicio los encontrados en
# el ajuste lineal

            print("--------------------------------")
            print("Frecuencias entre 1.0 - 2.0 Hz")

            if p1_02 < 0:
                p0_02 = (np.exp(p2_02),p1_02)
                popt_02, pcov_02 = curve_fit(model_exp, lpt02_grow, cut_Alpt02, p0_02)
                Q_02 = -np.pi*1/popt_02[1]
                print("Qc = {:.4f}".format(Q_02))
            
            elif p1_02 > 0:
                Q_02 = 0
                popt_02 = np.zeros(shape=2)
                print("No se puede calcular Qc")


#------------------------------------------------             
### Difusión

            print("--------------------------------")
            print("Frecuencias entre 1.0 - 2.0 Hz")
            
            if p1_02 < 0:
                
                lapsetime02_entero = t_evento + delta_max02 + t_02 # lapsetime desde el inicio del tag
                
                len_diferencia_lpt_02 = len(lapsetime02_entero)-len(lapsetime02) # veo cuántas muestas hay entre el lapsetime desde max amplt y el lapsetime desde el inicio
                ## la difetencia anterior está correcta xq lo verifique mirando grafico plt.plot(t_02,envolvente_02)
                
                cortar_lpt_final_02 = len_lpt02 + len_diferencia_lpt_02 # porque len_lpt esta considerado desde amplt max hasta el mejor ajuste
                lapsetime02_entero = lapsetime02_entero[:cortar_lpt_final_02]
                Alpt02_entero = envolvente_02[:cortar_lpt_final_02] * lapsetime02_entero
                
                popt_02_d, pcov_02_d = curve_fit(model_dif, lapsetime02_entero, np.log(Alpt02_entero))
                
                f02 = (w3+w4)/2
                Qi_02 = (-2)*np.pi*f02/(popt_02_d[1]) # temporal
                
                d_02 = (-1)*dist**2/(4*popt_02_d[2]) # coef de difusividad
                
                Qs_02 = (2)*np.pi*f02*p*d_02/vel**2 # espacial
                
                print("Qi={:.4f} y Qs={:.4f}".format(Qi_02,Qs_02))


            elif p1_02 > 0:
                Qi_02 = 0
                Qs_02 = 0
                f02 = (w3+w4)/2
                popt_02_d = np.zeros(shape=3)
                d_02 = 0
                print("No se puede calcular Qi ni Qs")



#------------------------------------------------  
### Grabado de datos

            field_names = ['UTCDateTime','Fc','lapsetime','dist','p_dg','vel','a_lin','b_lin','a_exp','b_exp','Qc_SBS','a1_dif','a2_dif','a3_dif','Qi','d','Qs']

            dict = {'UTCDateTime':inicio,'Fc':f02,'lapsetime':len_lpt02,'dist':dist,'p_dg':p,'vel':vel,'a_lin':np.exp(p2_02),'b_lin':-p1_02,'a_exp':popt_02[0],'b_exp':popt_02[1],'Qc_SBS': Q_02,'a1_dif':popt_02_d[0],'a2_dif':popt_02_d[1],'a3_dif':popt_02_d[2],'Qi':Qi_02,'d':d_02,'Qs':Qs_02}
            
            with open(path_archivo_salida, 'a') as csv_file:
                dict_object = csv.DictWriter(csv_file, fieldnames=field_names) 
              
                dict_object.writerow(dict)


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# 03 FILTRO PASA BANDA - Fc = 2Hz - TODOS LOS CÁLCULOS
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
                
### FIltro Pasa Banda a la señal           

            w5=1.5
            w6=2.5
            tr_taper10_copy03 = tr_taper10.copy()   # hago una copia para no sobreescribir
            evento_filtrado_03 = tr_taper10_copy03.filter("bandpass", freqmin=w5, freqmax=w6, corners=2, zerophase=False)


#------------------------------------------------
### Calculo FFT del evento previamente trabajado

            fft_03 = rfft(evento_filtrado_03)   # fft de BP e/ 0.5 - 1.5 Hz
      
              
#------------------------------------------------
### Cálculo de envolvente

            dom_time_03 = irfft(fft_03,n=len(evento[0]))     # cálculo de la transformada inversa (parte real)

            st_dom_time_03 = (new_st(dom_time_03)).slice(inicio,fin)
            t_03 = vec_tiempo(st_dom_time_03) # vector tiempo
            
            envolvente_03 = obspy.signal.filter.envelope(st_dom_time_03[0].data)
            
                        
#------------------------------------------------
### Valor máximo de cada envolvente y corte de la misma y del tiempo
            
            max_03 = posicion_max(envolvente_03)            
                     
            cut_env03 = envolvente_03[max_03:]

            cut_t03 = t_03[max_03:]
            
                       
#------------------------------------------------            
### Cálculo del tiempo correspondiente al máx valor de la envolvente

            t_max03 = inicio_pre + t_03[max_03]
            

#------------------------------------------------               
### Cálculo de vectores lapsetime

            delta_max03 = t_max03 - inicio
            lapsetime03 = t_evento + delta_max03 + cut_t03


###### CÁLCULOS PRIORIZANDO BUENOS AJUSTES:
#------------------------------------------------            
### Primer ajuste: lineal

            print("Resultados del ajuste lineal. Frecuencias entre 1.5 - 2.5 Hz")
            
            largo = 1500
            
            lpt03_grow = lapsetime03[:largo] ### acepto como mínimo 30 segundo a 50 mps
            cut_env03_grow = cut_env03[:largo]
            
            Alpt03 = cut_env03_grow * lpt03_grow

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
                cut_env03_grow = cut_env03[:largo]                

                Alpt03 = cut_env03_grow * lpt03_grow
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
                    #len_lpt03 = len(lpt03_grow) - (10*50)
                    len_lpt03 = len(lapsetime03[:(largo-(50*10))])
                    lpt03_grow = lapsetime03[:len_lpt03]                    
                    print("len_lpt03: {}, pendiente: {:.4f}".format(len_lpt03,p1_03))
                    break
                
                elif p1_03_b[0] < p1_03_a[0]:
                    lpt03_grow = lapsetime03[:largo]
                    cut_env03_grow = cut_env03[:largo]                    
                    p1_03_a[0] = p1_03_b[0] 
                    
            if len(lpt03_grow) == len(lapsetime03):
                print("estoy aca")
                p1_03 = p1_03_a[0]
                p2_03 = p2_03_a
                len_lpt03 = len(lpt03_grow)   
                p1_03_b = [0] # variable creada para que no salte error en el siguiente if
                p1_03_b = np.array(p1_03_b)                
                print("len_lpt03: {}, pendiente: {:.4f}".format(len_lpt03,p1_03))
            
            if p1_03_b[0] == p1_03_a[0]:
                p1_03 = p1_03_a[0]
                p2_03 = p2_03_a
                len_lpt03 = len(lpt03_grow)   
                print("len_lpt03: {}, pendiente: {:.4f}".format(len_lpt03,p1_03))
                
            cut_prediccion03 = prediccion03[:len_lpt03]
            cut_Alpt03 = Alpt03[:len_lpt03]
            cut_cut_env03 = cut_env03[:len_lpt03]
            cut_logAlpt03 = logAlpt03[:len_lpt03]


#------------------------------------------------ 
### Segundo ajuste: exponencial y cálculo del Q

# Se realiza el ajuste exponencial dando como valores de inicio los encontrados en
# el ajuste lineal

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
                popt_03 = np.zeros(shape=2)
                print("No se puede calcular Qc")


#------------------------------------------------             
### Difusión

            print("--------------------------------")
            print("Frecuencias entre 1.5 - 2.5 Hz")
            
            if p1_03 < 0:
                lapsetime03_entero = t_evento + delta_max03 + t_03
                len_diferencia_lpt_03 = len(lapsetime03_entero)-len(lapsetime03)
                cortar_lpt_final_03 = len_lpt03 + len_diferencia_lpt_03
                lapsetime03_entero = lapsetime03_entero[:cortar_lpt_final_03]
                Alpt03_entero = envolvente_03[:cortar_lpt_final_03] * lapsetime03_entero
                
                popt_03_d, pcov_03_d = curve_fit(model_dif, lapsetime03_entero, np.log(Alpt03_entero))
                
                f03 = (w5+w6)/2
                Qi_03 = (-2)*np.pi*f03/(popt_03_d[1]) # temporal
                
                d_03 = (-1)*dist**2/(4*popt_03_d[2]) # coef de difusividad
                
                Qs_03 = (2)*np.pi*f03*p*d_03/vel**2 # espacial

                print("Qi={:.4f} y Qs={:.4f}".format(Qi_03,Qs_03))
                
            
            elif p1_03 > 0:
                Qi_03 = 0
                Qs_03 = 0
                f03 = (w5+w6)/2
                popt_03_d = np.zeros(shape=3)
                d_03 = 0
                print("No se puede calcular Qi ni Qs")



#------------------------------------------------  
### Grabado de datos

            field_names = ['UTCDateTime','Fc','lapsetime','dist','p_dg','vel','a_lin','b_lin','a_exp','b_exp','Qc_SBS','a1_dif','a2_dif','a3_dif','Qi','d','Qs']

            #lpt=lapsetime01[-1:] # asi en cada lapse time, nos quedamos con el ultimo elemento como float para imprimirlo
            #lptfin=lpt[0]
            dict = {'UTCDateTime':inicio,'Fc':f03,'lapsetime':len_lpt03,'dist':dist,'p_dg':p,'vel':vel,'a_lin':np.exp(p2_03),'b_lin':-p1_03,'a_exp':popt_03[0],'b_exp':popt_03[1],'Qc_SBS': Q_03,'a1_dif':popt_03_d[0],'a2_dif':popt_03_d[1],'a3_dif':popt_03_d[2],'Qi':Qi_03,'d':d_03,'Qs':Qs_03}
            
            with open(path_archivo_salida, 'a') as csv_file:
                dict_object = csv.DictWriter(csv_file, fieldnames=field_names) 
              
                dict_object.writerow(dict)
                
                
               
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# 04 FILTRO PASA BANDA - Fc = 2.5Hz - TODOS LOS CÁLCULOS
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
                
### FIltro Pasa Banda a la señal           

            w7=2.0
            w8=3.0
            tr_taper10_copy04 = tr_taper10.copy()   # hago una copia para no sobreescribir
            evento_filtrado_04 = tr_taper10_copy04.filter("bandpass", freqmin=w7, freqmax=w8, corners=2, zerophase=False)


#------------------------------------------------
### Calculo FFT del evento previamente trabajado

            fft_04 = rfft(evento_filtrado_04)
      
              
#------------------------------------------------
### Cálculo de envolvente

            dom_time_04 = irfft(fft_04,n=len(evento[0]))     # cálculo de la transformada inversa (parte real)

            st_dom_time_04 = (new_st(dom_time_04)).slice(inicio,fin)
            t_04 = vec_tiempo(st_dom_time_04) # vector tiempo
            
            envolvente_04 = obspy.signal.filter.envelope(st_dom_time_04[0].data)
            
                        
#------------------------------------------------
### Valor máximo de cada envolvente y corte de la misma y del tiempo
            
            max_04 = posicion_max(envolvente_04)            
                     
            cut_env04 = envolvente_04[max_04:]

            cut_t04 = t_04[max_04:]
            
                       
#------------------------------------------------            
### Cálculo del tiempo correspondiente al máx valor de la envolvente

            t_max04 = inicio_pre + t_04[max_04]
            

#------------------------------------------------               
### Cálculo de vectores lapsetime

            delta_max04 = t_max04 - inicio
            lapsetime04 = t_evento + delta_max04 + cut_t04


###### CÁLCULOS PRIORIZANDO BUENOS AJUSTES:
#------------------------------------------------            
### Primer ajuste: lineal

            print("Resultados del ajuste lineal. Frecuencias entre 2.0 - 3.0 Hz")
            
            largo = 1500
            
            lpt04_grow = lapsetime04[:largo] ### acepto como mínimo 30 segundo a 50 mps
            cut_env04_grow = cut_env04[:largo]
            
            Alpt04 = cut_env04_grow * lpt04_grow

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
                cut_env04_grow = cut_env04[:largo]                

                Alpt04 = cut_env04_grow * lpt04_grow
                logAlpt04 = np.log(Alpt04).reshape(-1)
                lpt04_grow_reshape = lpt04_grow.reshape(-1,1)
                model04 = LinearRegression().fit(lpt04_grow_reshape, logAlpt04)
                prediccion04 = model04.predict(lpt04_grow_reshape)
                r_sq04 = model04.score(lpt04_grow_reshape, logAlpt04)
                p2_04_b = model04.intercept_
                p1_04_b = model04.coef_
                
                if p1_04_b[0] > p1_04_a[0]:                    
                    p1_04 = p1_04_a[0]
                    p2_04 = p2_04_a
                    len_lpt04 = len(lapsetime04[:(largo-(50*10))])
                    lpt04_grow = lapsetime04[:len_lpt04]
                    print("len_lpt04: {}, pendiente: {:.4f}".format(len_lpt04,p1_04))
                    break
                
                elif p1_04_b[0] < p1_04_a[0]:
                    lpt04_grow = lapsetime04[:largo]
                    cut_env04_grow = cut_env04[:largo]
                    p1_04_a[0] = p1_04_b[0]
                
            if len(lpt04_grow) == len(lapsetime04):
                print("estoy aca")
                p1_04 = p1_04_a[0]
                p2_04 = p2_04_a
                len_lpt04 = len(lpt04_grow)   
                p1_04_b = [0] # variable creada para que no salte error en el siguiente if
                p1_04_b = np.array(p1_04_b)                
                print("len_lpt04: {}, pendiente: {:.4f}".format(len_lpt04,p1_04))

            if p1_04_b[0] == p1_04_a[0]:
                p1_04 = p1_04_a[0]
                p2_04 = p2_04_a
                len_lpt04 = len(lpt04_grow)
                print("len_lpt04: {}, pendiente: {:.4f}".format(len_lpt04,p1_04))
  
            cut_prediccion04 = prediccion04[:len_lpt04]
            cut_Alpt04 = Alpt04[:len_lpt04]
            cut_cut_env04 = cut_env04[:len_lpt04]
            cut_logAlpt04 = logAlpt04[:len_lpt04]


#------------------------------------------------ 
### Segundo ajuste: exponencial y cálculo del Q

            print("--------------------------------")
            print("Frecuencias entre 2.0 - 3.0 Hz")
            
            if p1_04 < 0:
                p0_04 = (np.exp(p2_04),p1_04)
                popt_04, pcov_04 = curve_fit(model_exp, lpt04_grow, cut_Alpt04, p0_04)
                Q_04 = -np.pi*1/popt_04[1]
                print("Qc= {:.4f}".format(Q_04))
            
            elif p1_04 > 0:
                Q_04 = 0
                popt_04 = np.zeros(shape=2)
                print("No se puede calcular Qc")


#------------------------------------------------             
### Difusión

            print("--------------------------------")
            print("Frecuencias entre 2.0 - 3.0 Hz")
            
            if p1_04 < 0:
                lapsetime04_entero = t_evento + delta_max04 + t_04
                len_diferencia_lpt_04 = len(lapsetime04_entero)-len(lapsetime04)
                cortar_lpt_final_04 = len_lpt04 + len_diferencia_lpt_04
                lapsetime04_entero = lapsetime04_entero[:cortar_lpt_final_04]
                Alpt04_entero = envolvente_04[:cortar_lpt_final_04] * lapsetime04_entero
                
                popt_04_d, pcov_04_d = curve_fit(model_dif, lapsetime04_entero, np.log(Alpt04_entero))
                
                f04 = (w7+w8)/2
                Qi_04 = (-2)*np.pi*f04/(popt_04_d[1]) # temporal
                
                d_04 = (-1)*dist**2/(4*popt_04_d[2]) # coef de difusividad
                
                Qs_04 = (2)*np.pi*f04*p*d_04/vel**2 # espacial
                
                print("Qi={:.4f} y Qs={:.4f}".format(Qi_04,Qs_04))

            
            elif p1_04 > 0:
                Qi_04 = 0
                Qs_04 = 0
                f04 = (w7+w8)/2
                popt_04_d = np.zeros(shape=3)
                d_04 = 0
                print("No se puede calcular Qi ni Qs")

#------------------------------------------------  
### Grabado de datos

            field_names = ['UTCDateTime','Fc','lapsetime','dist','p_dg','vel','a_lin','b_lin','a_exp','b_exp','Qc_SBS','a1_dif','a2_dif','a3_dif','Qi','d','Qs']

            #lpt=lapsetime01[-1:] # asi en cada lapse time, nos quedamos con el ultimo elemento como float para imprimirlo
            #lptfin=lpt[0]
            dict = {'UTCDateTime':inicio,'Fc':f04,'lapsetime':len_lpt04,'dist':dist,'p_dg':p,'vel':vel,'a_lin':np.exp(p2_04),'b_lin':-p1_04,'a_exp':popt_04[0],'b_exp':popt_04[1],'Qc_SBS': Q_04,'a1_dif':popt_04_d[0],'a2_dif':popt_04_d[1],'a3_dif':popt_04_d[2],'Qi':Qi_04,'d':d_04,'Qs':Qs_04}
            
            with open(path_archivo_salida, 'a') as csv_file:
                dict_object = csv.DictWriter(csv_file, fieldnames=field_names) 
              
                dict_object.writerow(dict)
                


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# 05 FILTRO PASA BANDA - Fc = 3Hz - TODOS LOS CÁLCULOS
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
                
### FIltro Pasa Banda a la señal           

            w9=2.5
            w10=3.5
            tr_taper10_copy05 = tr_taper10.copy()   # hago una copia para no sobreescribir
            evento_filtrado_05 = tr_taper10_copy05.filter("bandpass", freqmin=w9, freqmax=w10, corners=2, zerophase=False)


#------------------------------------------------
### Calculo FFT del evento previamente trabajado

            fft_05 = rfft(evento_filtrado_05)
      
              
#------------------------------------------------
### Cálculo de envolvente

            dom_time_05 = irfft(fft_05,n=len(evento[0]))     # cálculo de la transformada inversa (parte real)

            st_dom_time_05 = (new_st(dom_time_05)).slice(inicio,fin)
            t_05 = vec_tiempo(st_dom_time_05) # vector tiempo
            
            envolvente_05 = obspy.signal.filter.envelope(st_dom_time_05[0].data)
            
                        
#------------------------------------------------
### Valor máximo de cada envolvente y corte de la misma y del tiempo
            
            max_05 = posicion_max(envolvente_05)            
                     
            cut_env05 = envolvente_05[max_05:]

            cut_t05 = t_05[max_05:]
            
                       
#------------------------------------------------            
### Cálculo del tiempo correspondiente al máx valor de la envolvente

            t_max05 = inicio_pre + t_05[max_05]
            

#------------------------------------------------               
### Cálculo de vectores lapsetime

            delta_max05 = t_max05 - inicio
            lapsetime05 = t_evento + delta_max05 + cut_t05


###### CÁLCULOS PRIORIZANDO BUENOS AJUSTES:
#------------------------------------------------            
### Primer ajuste: lineal

            print("Resultados del ajuste lineal. Frecuencias entre 2.5 - 3.5 Hz")
            
            largo = 1500
            
            lpt05_grow = lapsetime05[:largo] ### acepto como mínimo 30 segundo a 50 mps
            cut_env05_grow = cut_env05[:largo]
            
            Alpt05 = cut_env05_grow * lpt05_grow

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
                cut_env05_grow = cut_env05[:largo]                

                Alpt05 = cut_env05_grow * lpt05_grow
                logAlpt05 = np.log(Alpt05).reshape(-1)
                lpt05_grow_reshape = lpt05_grow.reshape(-1,1)
                model05 = LinearRegression().fit(lpt05_grow_reshape, logAlpt05)
                prediccion05 = model05.predict(lpt05_grow_reshape)
                r_sq05 = model05.score(lpt05_grow_reshape, logAlpt05)
                p2_05_b = model05.intercept_
                p1_05_b = model05.coef_

                if p1_05_b[0] > p1_05_a[0]:
                    p1_05 = p1_05_a[0]
                    p2_05 = p2_05_a
                    len_lpt05 = len(lapsetime05[:(largo-(50*10))])
                    lpt05_grow = lapsetime05[:len_lpt05]
                    print("len_lpt05: {}, pendiente: {:.4f}".format(len_lpt05,p1_05))
                    break
                
                elif p1_05_b[0] < p1_05_a[0]:
                    lpt05_grow = lapsetime05[:largo]
                    cut_env05_grow = cut_env05[:largo]
                    p1_05_a[0] = p1_05_b[0]
                
            if len(lpt05_grow) == len(lapsetime05):
                print("estoy aca")
                p1_05 = p1_05_a[0]
                p2_05 = p2_05_a
                len_lpt05 = len(lpt05_grow)   
                p1_05_b = [0] # variable creada para que no salte error en el siguiente if
                p1_05_b = np.array(p1_05_b)                
                print("len_lpt05: {}, pendiente: {:.4f}".format(len_lpt05,p1_05))

            if p1_05_b[0] == p1_05_a[0]:
                p1_05 = p1_05_a[0]
                p2_05 = p2_05_a
                len_lpt05 = len(lpt05_grow)   
                print("len_lpt05: {}, pendiente: {:.4f}".format(len_lpt05,p1_05))
                             
            cut_prediccion05 = prediccion05[:len_lpt05]
            cut_Alpt05 = Alpt05[:len_lpt05]
            cut_cut_env05 = cut_env05[:len_lpt05]
            cut_logAlpt05 = logAlpt05[:len_lpt05]

#------------------------------------------------ 
### Segundo ajuste: exponencial y cálculo del Q

            print("--------------------------------")
            print("Frecuencias entre 2.5 - 3.5 Hz")
            
            if p1_05 < 0:
                p0_05 = (np.exp(p2_05),p1_05)
                popt_05, pcov_05 = curve_fit(model_exp, lpt05_grow, cut_Alpt05, p0_05)
                Q_05 = -np.pi*1/popt_05[1]
                print("Qc = {:.4f}".format(Q_05))
            
            elif p1_05 > 0:
                Q_05 = 0
                popt_05 = np.zeros(shape=2)
                print("No se puede calcular Qc")


#------------------------------------------------             
### Difusión

            print("--------------------------------")
            print("Frecuencias entre 2.5 - 3.5 Hz")
            
            if p1_05 < 0:
                lapsetime05_entero = t_evento + delta_max05 + t_05
                len_diferencia_lpt_05 = len(lapsetime05_entero)-len(lapsetime05)
                cortar_lpt_final_05 = len_lpt05 + len_diferencia_lpt_05
                lapsetime05_entero = lapsetime05_entero[:cortar_lpt_final_05]
                Alpt05_entero = envolvente_05[:cortar_lpt_final_05] * lapsetime05_entero
                
                popt_05_d, pcov_05_d = curve_fit(model_dif, lapsetime05_entero, np.log(Alpt05_entero))
                               
                f05 = (w9+w10)/2
                Qi_05 = (-2)*np.pi*f05/(popt_05_d[1]) # temporal
                
                d_05 = (-1)*dist**2/(4*popt_05_d[2]) # coef de difusividad
                
                Qs_05 = (2)*np.pi*f05*p*d_05/vel**2 # espacial
                
                print("Qi={:.4f} y Qs={:.4f}".format(Qi_05,Qs_05))
                

            elif p1_05 > 0:
                Qi_05 = 0
                Qs_05 = 0
                f05 = (w9+w10)/2
                popt_05_d = np.zeros(shape=3)
                d_05 = 0
                print("No se puede calcular Qi ni Qs")

#------------------------------------------------  
### Grabado de datos

            field_names = ['UTCDateTime','Fc','lapsetime','dist','p_dg','vel','a_lin','b_lin','a_exp','b_exp','Qc_SBS','a1_dif','a2_dif','a3_dif','Qi','d','Qs']

            #lpt=lapsetime01[-1:] # asi en cada lapse time, nos quedamos con el ultimo elemento como float para imprimirlo
            #lptfin=lpt[0]
            dict = {'UTCDateTime':inicio,'Fc':f05,'lapsetime':len_lpt05,'dist':dist,'p_dg':p,'vel':vel,'a_lin':np.exp(p2_05),'b_lin':-p1_05,'a_exp':popt_05[0],'b_exp':popt_05[1],'Qc_SBS': Q_05,'a1_dif':popt_05_d[0],'a2_dif':popt_05_d[1],'a3_dif':popt_05_d[2],'Qi':Qi_05,'d':d_05,'Qs':Qs_05}
            
            with open(path_archivo_salida, 'a') as csv_file:
                dict_object = csv.DictWriter(csv_file, fieldnames=field_names) 
              
                dict_object.writerow(dict)
                


#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# 06 FILTRO PASA BANDA - Fc = 3.5Hz - TODOS LOS CÁLCULOS
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------                
### FIltro Pasa Banda a la señal           

            w11=3.0
            w12=4.0
            tr_taper10_copy06 = tr_taper10.copy()   # hago una copia para no sobreescribir
            evento_filtrado_06 = tr_taper10_copy06.filter("bandpass", freqmin=w11, freqmax=w12, corners=2, zerophase=False)


#------------------------------------------------
### Calculo FFT del evento previamente trabajado

            fft_06 = rfft(evento_filtrado_06)
      
              
#------------------------------------------------
### Cálculo de envolvente

            dom_time_06 = irfft(fft_06,n=len(evento[0]))     # cálculo de la transformada inversa (parte real)

            st_dom_time_06 = (new_st(dom_time_06)).slice(inicio,fin)
            t_06 = vec_tiempo(st_dom_time_06) # vector tiempo
            
            envolvente_06 = obspy.signal.filter.envelope(st_dom_time_06[0].data)
            
                        
#------------------------------------------------
### Valor máximo de cada envolvente y corte de la misma y del tiempo
            
            max_06 = posicion_max(envolvente_06)            
                     
            cut_env06 = envolvente_06[max_06:]

            cut_t06 = t_06[max_06:]
            
                       
#------------------------------------------------            
### Cálculo del tiempo correspondiente al máx valor de la envolvente

            t_max06 = inicio_pre + t_06[max_06]
            

#------------------------------------------------               
### Cálculo de vectores lapsetime

            delta_max06 = t_max06 - inicio
            lapsetime06 = t_evento + delta_max06 + cut_t06


###### CÁLCULOS PRIORIZANDO BUENOS AJUSTES:
#------------------------------------------------            
### Primer ajuste: lineal

            print("Resultados del ajuste lineal. Frecuencias entre 3.0 - 4.0 Hz")
            
            largo = 1500
            
            lpt06_grow = lapsetime06[:largo] ### acepto como mínimo 30 segundo a 50 mps
            cut_env06_grow = cut_env06[:largo]
            
            Alpt06 = cut_env06_grow * lpt06_grow

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
                cut_env06_grow = cut_env06[:largo]                

                Alpt06 = cut_env06_grow * lpt06_grow
                logAlpt06 = np.log(Alpt06).reshape(-1)
                lpt06_grow_reshape = lpt06_grow.reshape(-1,1)
                model06 = LinearRegression().fit(lpt06_grow_reshape, logAlpt06)
                prediccion06 = model06.predict(lpt06_grow_reshape)
                r_sq06 = model06.score(lpt06_grow_reshape, logAlpt06)
                p2_06_b = model06.intercept_
                p1_06_b = model06.coef_
                
                if p1_06_b[0] > p1_06_a[0]:                    
                    p1_06 = p1_06_a[0]
                    p2_06 = p2_06_a
                    len_lpt06 = len(lapsetime06[:(largo-(50*10))])
                    lpt06_grow = lapsetime06[:len_lpt06]
                    print("len_lpt06: {}, pendiente: {:.4f}".format(len_lpt06,p1_06))
                    break
                
                elif p1_06_b[0] < p1_06_a[0]:
                    lpt06_grow = lapsetime06[:largo]
                    cut_env06_grow = cut_env06[:largo]
                    p1_06_a[0] = p1_06_b[0] 
  
            if len(lpt06_grow) == len(lapsetime06):
                print("estoy aca")
                p1_06 = p1_06_a[0]
                p2_06 = p2_06_a
                len_lpt06 = len(lpt06_grow)   
                p1_06_b = [0] # variable creada para que no salte error en el siguiente if
                p1_06_b = np.array(p1_06_b)                
                print("len_lpt06: {}, pendiente: {:.4f}".format(len_lpt06,p1_06))

            if p1_06_b[0] == p1_06_a[0]:
                p1_06 = p1_06_a[0]
                p2_06 = p2_06_a
                len_lpt06 = len(lpt06_grow)   
                print("len_lpt06: {}, pendiente: {:.4f}".format(len_lpt06,p1_06))
               
            cut_prediccion06 = prediccion06[:len_lpt06]
            cut_Alpt06 = Alpt06[:len_lpt06]
            cut_cut_env06 = cut_env06[:len_lpt06]
            cut_logAlpt06 = logAlpt06[:len_lpt06]
            
            
#------------------------------------------------ 
### Segundo ajuste: exponencial y cálculo del Q

            print("--------------------------------")
            print("Frecuencias entre 3.0 - 4.0 Hz")
            
            if p1_06 < 0:
                p0_06 = (np.exp(p2_06),p1_06)
                popt_06, pcov_06 = curve_fit(model_exp, lpt06_grow, cut_Alpt06, p0_06)
                Q_06 = -np.pi*1/popt_06[1]
                print("Qc = {:.4f}".format(Q_06))
            
            elif p1_06 > 0:
                Q_06 = 0
                popt_06 = np.zeros(shape=2)
                print("No se puede calcular Qc")
            

#------------------------------------------------             
### Difusión

            print("--------------------------------")
            print("Frecuencias entre 3.0 - 4.0 Hz")

            if p1_06 < 0:
                lapsetime06_entero = t_evento + delta_max06 + t_06
                len_diferencia_lpt_06 = len(lapsetime06_entero)-len(lapsetime06)
                cortar_lpt_final_06 = len_lpt06 + len_diferencia_lpt_06
                lapsetime06_entero = lapsetime06_entero[:cortar_lpt_final_06]
                Alpt06_entero = envolvente_06[:cortar_lpt_final_06] * lapsetime06_entero
                
                popt_06_d, pcov_06_d = curve_fit(model_dif, lapsetime06_entero, np.log(Alpt06_entero))
                
                f06 = (w11+w12)/2
                Qi_06 = (-2)*np.pi*f06/(popt_06_d[1]) # temporal
                
                d_06 = (-1)*dist**2/(4*popt_06_d[2]) # coef de difusividad
                
                Qs_06 = (2)*np.pi*f06*p*d_06/vel**2 # espacial
                
                print("Qi={:.4f} y Qs={:.4f}".format(Qi_06,Qs_06))
                

            elif p1_06 > 0:
                Qi_06 = 0
                Qs_06 = 0
                f06 = (w11+w12)/2
                popt_06_d = np.zeros(shape=3)
                d_06 = 0
                print("No se puede calcular Qi ni Qs")


#------------------------------------------------  
### Grabado de datos

            field_names = ['UTCDateTime','Fc','lapsetime','dist','p_dg','vel','a_lin','b_lin','a_exp','b_exp','Qc_SBS','a1_dif','a2_dif','a3_dif','Qi','d','Qs']

            #lpt=lapsetime01[-1:] # asi en cada lapse time, nos quedamos con el ultimo elemento como float para imprimirlo
            #lptfin=lpt[0]
            dict = {'UTCDateTime':inicio,'Fc':f06,'lapsetime':len_lpt06,'dist':dist,'p_dg':p,'vel':vel,'a_lin':np.exp(p2_06),'b_lin':-p1_06,'a_exp':popt_06[0],'b_exp':popt_06[1],'Qc_SBS': Q_06,'a1_dif':popt_06_d[0],'a2_dif':popt_06_d[1],'a3_dif':popt_06_d[2],'Qi':Qi_06,'d':d_06,'Qs':Qs_06}
            
            with open(path_archivo_salida, 'a') as csv_file:
                dict_object = csv.DictWriter(csv_file, fieldnames=field_names) 
              
                dict_object.writerow(dict)
                
                
              
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# 07 FILTRO PASA BANDA - Fc = 4Hz - TODOS LOS CÁLCULOS
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------                
### FIltro Pasa Banda a la señal           

            w13=3.5
            w14=4.5
            tr_taper10_copy07 = tr_taper10.copy()   # hago una copia para no sobreescribir
            evento_filtrado_07 = tr_taper10_copy07.filter("bandpass", freqmin=w13, freqmax=w14, corners=2, zerophase=False)


#------------------------------------------------
### Calculo FFT del evento previamente trabajado

            fft_07 = rfft(evento_filtrado_07)
      
              
#------------------------------------------------
### Cálculo de envolvente

            dom_time_07 = irfft(fft_07,n=len(evento[0]))     # cálculo de la transformada inversa (parte real)

            st_dom_time_07 = (new_st(dom_time_07)).slice(inicio,fin)
            t_07 = vec_tiempo(st_dom_time_07) # vector tiempo
            
            envolvente_07 = obspy.signal.filter.envelope(st_dom_time_07[0].data)
            
                        
#------------------------------------------------
### Valor máximo de cada envolvente y corte de la misma y del tiempo
            
            max_07 = posicion_max(envolvente_07)            
                     
            cut_env07 = envolvente_07[max_07:]

            cut_t07 = t_07[max_07:]
            
                       
#------------------------------------------------            
### Cálculo del tiempo correspondiente al máx valor de la envolvente

            t_max07 = inicio_pre + t_07[max_07]
            

#------------------------------------------------               
### Cálculo de vectores lapsetime

            delta_max07 = t_max07 - inicio
            lapsetime07 = t_evento + delta_max07 + cut_t07


###### CÁLCULOS PRIORIZANDO BUENOS AJUSTES:
#------------------------------------------------            
### Primer ajuste: lineal

            print("Resultados del ajuste lineal. Frecuencias entre 3.5 - 4.5 Hz")
            
            largo = 1500
            
            lpt07_grow = lapsetime07[:largo] ### acepto como mínimo 30 segundo a 50 mps
            cut_env07_grow = cut_env07[:largo]     
            
            Alpt07 = cut_env07_grow * lpt07_grow
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
                cut_env07_grow = cut_env07[:largo]  
                
                Alpt07 = cut_env07_grow * lpt07_grow
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
                    len_lpt07 = len(lapsetime07[:(largo-(50*10))])
                    lpt07_grow = lapsetime07[:len_lpt07]
                    print("len_lpt07: {}, pendiente: {:.4f}".format(len_lpt07,p1_07))
                    break
                
                elif p1_07_b[0] < p1_07_a[0]:
                    lpt07_grow = lapsetime07[:largo]
                    cut_env07_grow = cut_env07[:largo]
                    p1_07_a[0] = p1_07_b[0] 
                    
            if len(lpt07_grow) == len(lapsetime07):
                print("estoy aca")
                p1_07 = p1_07_a[0]
                p2_07 = p2_07_a
                len_lpt07 = len(lpt07_grow)   
                p1_07_b = [0] # variable creada para que no salte error en el siguiente if
                p1_07_b = np.array(p1_06_b)                
                print("len_lpt07: {}, pendiente: {:.4f}".format(len_lpt07,p1_07))
                
            if p1_07_b[0] == p1_07_a[0]:
                p1_07 = p1_07_a[0]
                p2_07 = p2_07_a
                len_lpt07 = len(lpt07_grow)   
                print("len_lpt07: {}, pendiente: {:.4f}".format(len_lpt07,p1_07))
                
            cut_prediccion07 = prediccion07[:len_lpt07]
            cut_Alpt07 = Alpt07[:len_lpt07]
            cut_cut_env07 = cut_env07[:len_lpt07]
            cut_logAlpt07 = logAlpt07[:len_lpt07]
            
            
#------------------------------------------------ 
### Segundo ajuste: exponencial y cálculo del Q

            print("--------------------------------")
            print("Frecuencias entre 3.5 - 4.5 Hz")
            
            if p1_07 < 0:
                p0_07 = (np.exp(p2_07),p1_07)
                popt_07, pcov_07 = curve_fit(model_exp, lpt07_grow, cut_Alpt07, p0_07)
                Q_07 = -np.pi*1/popt_07[1]
                print("Qc = {:.4f}".format(Q_07))
            
            elif p1_07 > 0:
                Q_07 = 0
                popt_07 = np.zeros(shape=2)
                print("No se puede calcular Qc")
            

#------------------------------------------------             
### Difusión

            print("--------------------------------")
            print("Frecuencias entre 3.5 - 4.5 Hz")
            
            if p1_07 < 0:
                lapsetime07_entero = t_evento + delta_max07 + t_07
                len_diferencia_lpt_07 = len(lapsetime07_entero)-len(lapsetime07)
                cortar_lpt_final_07 = len_lpt07 + len_diferencia_lpt_07
                lapsetime07_entero = lapsetime07_entero[:cortar_lpt_final_07]
                Alpt07_entero = envolvente_07[:cortar_lpt_final_07] * lapsetime07_entero
                
                popt_07_d, pcov_07_d = curve_fit(model_dif, lapsetime07_entero, np.log(Alpt07_entero))
                
                f07 = (w13+w14)/2
                Qi_07 = (-2)*np.pi*f07/(popt_07_d[1]) # temporal
                
                d_07 = (-1)*dist**2/(4*popt_07_d[2]) # coef de difusividad
                
                Qs_07 = (2)*np.pi*f07*p*d_07/vel**2 # espacial
                
                print("Qi={:.4f} y Qs={:.4f}".format(Qi_07,Qs_07))
                

            elif p1_07 > 0:
                Qi_07 = 0
                Qs_07 = 0
                f07 = (w13+w14)/2
                popt_07_d = np.zeros(shape=3)
                d_07 = 0
                print("No se puede calcular Qi ni Qs")


#------------------------------------------------  
### Grabado de datos

            field_names = ['UTCDateTime','Fc','lapsetime','dist','p_dg','vel','a_lin','b_lin','a_exp','b_exp','Qc_SBS','a1_dif','a2_dif','a3_dif','Qi','d','Qs']

            #lpt=lapsetime01[-1:] # asi en cada lapse time, nos quedamos con el ultimo elemento como float para imprimirlo
            #lptfin=lpt[0]
            dict = {'UTCDateTime':inicio,'Fc':f07,'lapsetime':len_lpt07,'dist':dist,'p_dg':p,'vel':vel,'a_lin':np.exp(p2_07),'b_lin':-p1_07,'a_exp':popt_07[0],'b_exp':popt_07[1],'Qc_SBS': Q_07,'a1_dif':popt_07_d[0],'a2_dif':popt_07_d[1],'a3_dif':popt_07_d[2],'Qi':Qi_07,'d':d_07,'Qs':Qs_07}
            
            with open(path_archivo_salida, 'a') as csv_file:
                dict_object = csv.DictWriter(csv_file, fieldnames=field_names) 
              
                dict_object.writerow(dict)
                
                
               
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# 08 FILTRO PASA BANDA - Fc = 4.5Hz - TODOS LOS CÁLCULOS
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------                
### FIltro Pasa Banda a la señal           

            w15=4.0
            w16=5.0
            tr_taper10_copy08 = tr_taper10.copy()   # hago una copia para no sobreescribir
            evento_filtrado_08 = tr_taper10_copy08.filter("bandpass", freqmin=w15, freqmax=w16, corners=2, zerophase=False)


#------------------------------------------------
### Calculo FFT del evento previamente trabajado

            fft_08 = rfft(evento_filtrado_08)
      
              
#------------------------------------------------
### Cálculo de envolvente

            dom_time_08 = irfft(fft_08,n=len(evento[0]))     # cálculo de la transformada inversa (parte real)

            st_dom_time_08 = (new_st(dom_time_08)).slice(inicio,fin)
            t_08 = vec_tiempo(st_dom_time_08) # vector tiempo
            
            envolvente_08 = obspy.signal.filter.envelope(st_dom_time_08[0].data)
            
                        
#------------------------------------------------
### Valor máximo de cada envolvente y corte de la misma y del tiempo
            
            max_08 = posicion_max(envolvente_08)            
                     
            cut_env08 = envolvente_08[max_08:]

            cut_t08 = t_08[max_08:]
            
                       
#------------------------------------------------            
### Cálculo del tiempo correspondiente al máx valor de la envolvente

            t_max08 = inicio_pre + t_08[max_08]
            

#------------------------------------------------               
### Cálculo de vectores lapsetime

            delta_max08 = t_max08 - inicio
            lapsetime08 = t_evento + delta_max08 + cut_t08


###### CÁLCULOS PRIORIZANDO BUENOS AJUSTES:
#------------------------------------------------            
### Primer ajuste: lineal

            print("Resultados del ajuste lineal. Frecuencias entre 4.0 - 5.0 Hz")
            
            largo = 1500
            
            lpt08_grow = lapsetime08[:largo] ### acepto como mínimo 30 segundo a 50 mps
            cut_env08_grow = cut_env08[:largo]         
            
            Alpt08 = cut_env08_grow * lpt08_grow
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
                cut_env08_grow = cut_env08[:largo] 
                
                Alpt08 = cut_env08_grow * lpt08_grow
                logAlpt08 = np.log(Alpt08).reshape(-1)
                lpt08_grow_reshape = lpt08_grow.reshape(-1,1)
                model08 = LinearRegression().fit(lpt08_grow_reshape, logAlpt08)
                prediccion08 = model08.predict(lpt08_grow_reshape)
                r_sq08 = model08.score(lpt08_grow_reshape, logAlpt08)
                p2_08_b = model08.intercept_
                p1_08_b = model08.coef_
                
                if p1_08_b[0] > p1_08_a[0]:                    
                    p1_08 = p1_08_a[0]
                    p2_08 = p2_08_a
                    len_lpt08 = len(lapsetime08[:(largo-(50*10))])
                    lpt08_grow = lapsetime08[:len_lpt08]
                    print("len_lpt08: {}, pendiente: {:.4f}".format(len_lpt08,p1_08))
                    break
                
                elif p1_08_b[0] < p1_08_a[0]:
                    #largo = largo + 10*50
                    lpt08_grow = lapsetime08[:largo]
                    cut_env08_grow = cut_env08[:largo]
                    p1_08_a[0] = p1_08_b[0] 

            if len(lpt08_grow) == len(lapsetime08):
                print("estoy aca")
                p1_08 = p1_08_a[0]
                p2_08 = p2_08_a
                len_lpt08 = len(lpt08_grow)   
                p1_08_b = [0] # variable creada para que no salte error en el siguiente if
                p1_08_b = np.array(p1_08_b)                
                print("len_lpt08: {}, pendiente: {:.4f}".format(len_lpt08,p1_08))

            if p1_08_b[0] == p1_08_a[0]:
                p1_08 = p1_08_a[0]
                p2_08 = p2_08_a
                len_lpt08 = len(lpt08_grow)   
                print("len_lpt08: {}, pendiente: {:.4f}".format(len_lpt08,p1_08))
                
            cut_prediccion08 = prediccion08[:len_lpt08]
            cut_Alpt08 = Alpt08[:len_lpt08]
            cut_cut_env08 = cut_env08[:len_lpt08]
            cut_logAlpt08 = logAlpt08[:len_lpt08]
            
            
#------------------------------------------------ 
### Segundo ajuste: exponencial y cálculo del Q

            print("--------------------------------")
            print("Frecuencias entre 4.0 - 5.0 Hz")
            
            if p1_08 < 0:
                p0_08 = (np.exp(p2_08),p1_08)
                popt_08, pcov_08 = curve_fit(model_exp, lpt08_grow, cut_Alpt08, p0_08)
                Q_08 = -np.pi*1/popt_08[1]
                print("Qc = {:.4f}".format(Q_08))
            
            elif p1_08 > 0:
                Q_08 = 0
                popt_08 = np.zeros(shape=2)
                print("No se puede calcular Qc")
            

#------------------------------------------------             
### Difusión

            print("--------------------------------")
            print("Frecuencias entre 4.0 - 5.0 Hz")
            
            if p1_08 < 0:
                lapsetime08_entero = t_evento + delta_max08 + t_08
                len_diferencia_lpt_08 = len(lapsetime08_entero)-len(lapsetime08)
                cortar_lpt_final_08 = len_lpt08 + len_diferencia_lpt_08
                lapsetime08_entero = lapsetime08_entero[:cortar_lpt_final_08]
                Alpt08_entero = envolvente_08[:cortar_lpt_final_08] * lapsetime08_entero
                
                popt_08_d, pcov_08_d = curve_fit(model_dif, lapsetime08_entero, np.log(Alpt08_entero))
                
                f08 = (w15+w16)/2
                Qi_08 = (-2)*np.pi*f08/(popt_08_d[1]) # temporal
                
                d_08 = (-1)*dist**2/(4*popt_08_d[2]) # coef de difusividad
                
                Qs_08 = (2)*np.pi*f08*p*d_08/vel**2 # espacial
                
                print("Qi={:.4f} y Qs={:.4f}".format(Qi_08,Qs_08))
                

            elif p1_08 > 0:
                Qi_08 = 0
                Qs_08 = 0
                f08 = (w15+w16)/2
                popt_08_d = np.zeros(shape=3)
                d_08 = 0
                print("No se puede calcular Qi ni Qs")


#------------------------------------------------  
### Grabado de datos

            field_names = ['UTCDateTime','Fc','lapsetime','dist','p_dg','vel','a_lin','b_lin','a_exp','b_exp','Qc_SBS','a1_dif','a2_dif','a3_dif','Qi','d','Qs']

            #lpt=lapsetime01[-1:] # asi en cada lapse time, nos quedamos con el ultimo elemento como float para imprimirlo
            #lptfin=lpt[0]
            dict = {'UTCDateTime':inicio,'Fc':f08,'lapsetime':len_lpt08,'dist':dist,'p_dg':p,'vel':vel,'a_lin':np.exp(p2_08),'b_lin':-p1_08,'a_exp':popt_08[0],'b_exp':popt_08[1],'Qc_SBS': Q_08,'a1_dif':popt_08_d[0],'a2_dif':popt_08_d[1],'a3_dif':popt_08_d[2],'Qi':Qi_08,'d':d_08,'Qs':Qs_08}
            
            with open(path_archivo_salida, 'a') as csv_file:
                dict_object = csv.DictWriter(csv_file, fieldnames=field_names) 
              
                dict_object.writerow(dict)
                
#%%
#------------------------------------------------
### Figuras

            '''fecha = inicio_string[0:4] + '_' + inicio_string[5:7] + '_' + inicio_string[8:10] + '_' + inicio_string[11:13] + '_' + inicio_string[14:16] + '_' + inicio_string[17:19]
            path_fig_dif = 'C:/Users/Propietario/Desktop/tesis_de_grado/ATENUACION/DATOS/FG16/resultados/figuras_dif'

# envolvente y ajuste dif
            name_dif01 = '{}_ajuste_dif_{}_Hz.png'.format(fecha,f01)
            ruta01 = os.path.join(path_fig_dif,name_dif01)
            if p1_01 < 0:
                plt.figure()
                plt.title(inicio)
                plt.plot(lapsetime01_entero,Alpt01_entero,color='black',label='envolvente fc={}Hz'.format(f01))                
                plt.plot(lapsetime01_entero,ajuste_dif(popt_01_d[0],popt_01_d[1],popt_01_d[2],lapsetime01_entero),color='violet',label="ajuste difusión")
                plt.xlabel('tiempo de viaje [seg]')
                plt.legend()
                plt.savefig(ruta01)
                plt.close()

            name_dif02 = '{}_ajuste_dif_{}_Hz.png'.format(fecha,f02)
            ruta02 = os.path.join(path_fig_dif,name_dif02)                
            if p1_02 < 0:
                plt.figure()
                plt.title(inicio)
                plt.plot(lapsetime02_entero,Alpt02_entero,color='black',label='envolvente fc={}Hz'.format(f02))                
                plt.plot(lapsetime02_entero,ajuste_dif(popt_02_d[0],popt_02_d[1],popt_02_d[2],lapsetime02_entero),color='violet',label="ajuste difusión")
                plt.xlabel('tiempo de viaje [seg]')
                plt.legend()
                plt.savefig(ruta02)
                plt.close()

            name_dif03 = '{}_ajuste_dif_{}_Hz.png'.format(fecha,f03)
            ruta03 = os.path.join(path_fig_dif,name_dif03)                
            if p1_03 < 0:
                plt.figure()
                plt.title(inicio)
                plt.plot(lapsetime03_entero,Alpt03_entero,color='black',label='envolvente fc={}Hz'.format(f03))                
                plt.plot(lapsetime03_entero,ajuste_dif(popt_03_d[0],popt_03_d[1],popt_03_d[2],lapsetime03_entero),color='violet',label="ajuste difusión")
                plt.xlabel('tiempo de viaje [seg]')
                plt.legend()
                plt.savefig(ruta03)
                plt.close()

            name_dif04 = '{}_ajuste_dif_{}_Hz.png'.format(fecha,f04)
            ruta04 = os.path.join(path_fig_dif,name_dif04)                                
            if p1_04 < 0:
                plt.figure()
                plt.title(inicio)
                plt.plot(lapsetime04_entero,Alpt04_entero,color='black',label='envolvente fc={}Hz'.format(f04))                
                plt.plot(lapsetime04_entero,ajuste_dif(popt_04_d[0],popt_04_d[1],popt_04_d[2],lapsetime04_entero),color='violet',label="ajuste difusión")
                plt.xlabel('tiempo de viaje [seg]')
                plt.legend()
                plt.savefig(ruta04)
                plt.close()

            name_dif05 = '{}_ajuste_dif_{}_Hz.png'.format(fecha,f05)
            ruta05 = os.path.join(path_fig_dif,name_dif05)                    
            if p1_05 < 0:
                plt.figure()
                plt.title(inicio)
                plt.plot(lapsetime05_entero,Alpt05_entero,color='black',label='envolvente fc={}Hz'.format(f05))                
                plt.plot(lapsetime05_entero,ajuste_dif(popt_05_d[0],popt_05_d[1],popt_05_d[2],lapsetime05_entero),color='violet',label="ajuste difusión")
                plt.xlabel('tiempo de viaje [seg]')
                plt.legend()
                plt.savefig(ruta05)
                plt.close()

            name_dif06 = '{}_ajuste_dif_{}_Hz.png'.format(fecha,f06)
            ruta06 = os.path.join(path_fig_dif,name_dif06)                    
            if p1_06 < 0:
                plt.figure()
                plt.title(inicio)
                plt.plot(lapsetime06_entero,Alpt06_entero,color='black',label='envolvente fc={}Hz'.format(f06))                
                plt.plot(lapsetime06_entero,ajuste_dif(popt_06_d[0],popt_06_d[1],popt_06_d[2],lapsetime06_entero),color='violet',label="ajuste difusión")
                plt.xlabel('tiempo de viaje [seg]')
                plt.legend()
                plt.savefig(ruta06)
                plt.close()

            name_dif07 = '{}_ajuste_dif_{}_Hz.png'.format(fecha,f07)
            ruta07 = os.path.join(path_fig_dif,name_dif07)                    
            if p1_07 < 0:
                plt.figure()
                plt.title(inicio)
                plt.plot(lapsetime07_entero,Alpt07_entero,color='black',label='envolvente fc={}Hz'.format(f07))                
                plt.plot(lapsetime07_entero,ajuste_dif(popt_07_d[0],popt_07_d[1],popt_07_d[2],lapsetime07_entero),color='violet',label="ajuste difusión")
                plt.xlabel('tiempo de viaje [seg]')
                plt.legend()
                plt.savefig(ruta07)
                plt.close()

            name_dif08 = '{}_ajuste_dif_{}_Hz.png'.format(fecha,f08)
            ruta08 = os.path.join(path_fig_dif,name_dif08)
            if p1_08 < 0:
                plt.figure()
                plt.title(inicio)
                plt.plot(lapsetime08_entero,Alpt08_entero,color='black',label='envolvente fc={}Hz'.format(f08))                
                plt.plot(lapsetime08_entero,ajuste_dif(popt_08_d[0],popt_08_d[1],popt_08_d[2],lapsetime08_entero),color='violet',label="ajuste difusión")
                plt.xlabel('tiempo de viaje [seg]')
                plt.legend()
                plt.savefig(ruta08)
                plt.close()


# envolvente y ajuste qc

            path_fig_exp = 'C:/Users/Propietario/Desktop/tesis_de_grado/ATENUACION/DATOS/FG16/resultados/figuras_exp'
            
            name_exp01 = '{}_ajuste_exp_{}_Hz.png'.format(fecha,f01)
            ruta01 = os.path.join(path_fig_exp,name_exp01)            
            lapsetime01_entero = t_evento + t_01
            if p1_01 < 0:
                plt.figure()
                plt.title(inicio)
                plt.plot(lapsetime01,Alpt01,color='black',label='envolvente corregida fc={}Hz'.format(f01))                
                plt.plot(lpt01_grow,model_exp(lpt01_grow,popt_01[0],popt_01[1]),color='violet',label="ajuste no lineal")
                plt.plot(lapsetime01_entero,st_dom_time_01[0])
                plt.xlabel('tiempo de viaje [seg]')
                plt.legend()
                #plt.show()
                plt.savefig(ruta01)
                plt.close()

            name_exp02 = '{}_ajuste_exp_{}_Hz.png'.format(fecha,f02)
            ruta02 = os.path.join(path_fig_exp,name_exp02)            
            lapsetime02_entero = t_evento + t_02
            if p1_02 < 0:
                plt.figure()
                plt.title(inicio)
                plt.plot(lapsetime02,Alpt02,color='black',label='envolvente corregida fc={}Hz'.format(f02))                
                plt.plot(lpt02_grow,model_exp(lpt02_grow,popt_02[0],popt_02[1]),color='violet',label="ajuste no lineal")
                plt.plot(lapsetime02_entero,st_dom_time_02[0])
                plt.xlabel('tiempo de viaje [seg]')
                plt.legend()
                #plt.show()
                plt.savefig(ruta02)
                plt.close()

            name_exp03 = '{}_ajuste_exp_{}_Hz.png'.format(fecha,f03)
            ruta03 = os.path.join(path_fig_exp,name_exp03) 
            lapsetime03_entero = t_evento + t_03
            if p1_03 < 0:
                plt.figure()
                plt.title(inicio)
                plt.plot(lapsetime03,Alpt03,color='black',label='envolvente corregida fc={}Hz'.format(f03))                
                plt.plot(lpt03_grow,model_exp(lpt03_grow,popt_03[0],popt_03[1]),color='violet',label="ajuste no lineal")
                plt.plot(lapsetime03_entero,st_dom_time_03[0])
                plt.xlabel('tiempo de viaje [seg]')
                plt.legend()
                #plt.show()
                plt.savefig(ruta03)
                plt.close()

            name_exp04 = '{}_ajuste_exp_{}_Hz.png'.format(fecha,f04)
            ruta04 = os.path.join(path_fig_exp,name_exp04)                
            lapsetime04_entero = t_evento + t_04
            if p1_04 < 0:
                plt.figure()
                plt.title(inicio)
                plt.plot(lapsetime04,Alpt04,color='black',label='envolvente corregida fc={}Hz'.format(f04))                
                plt.plot(lpt04_grow,model_exp(lpt04_grow,popt_04[0],popt_04[1]),color='violet',label="ajuste no lineal")
                plt.plot(lapsetime04_entero,st_dom_time_04[0])
                plt.xlabel('tiempo de viaje [seg]')
                plt.legend()
                #plt.show()
                plt.savefig(ruta04)
                plt.close()

            name_exp05 = '{}_ajuste_exp_{}_Hz.png'.format(fecha,f05)
            ruta05 = os.path.join(path_fig_exp,name_exp05)                 
            lapsetime05_entero = t_evento + t_05
            if p1_05 < 0:
                plt.figure()
                plt.title(inicio)
                plt.plot(lapsetime05,Alpt05,color='black',label='envolvente corregida fc={}Hz'.format(f05))                
                plt.plot(lpt05_grow,model_exp(lpt05_grow,popt_05[0],popt_05[1]),color='violet',label="ajuste no lineal")
                plt.plot(lapsetime05_entero,st_dom_time_05[0])
                plt.xlabel('tiempo de viaje [seg]')
                plt.legend()
                #plt.show()
                plt.savefig(ruta05)
                plt.close()

            name_exp06 = '{}_ajuste_exp_{}_Hz.png'.format(fecha,f06)
            ruta06 = os.path.join(path_fig_exp,name_exp06)
            lapsetime06_entero = t_evento + t_06
            if p1_06 < 0:
                plt.figure()
                plt.title(inicio)
                plt.plot(lapsetime06,Alpt06,color='black',label='envolvente corregida fc={}Hz'.format(f06))                
                plt.plot(lpt06_grow,model_exp(lpt06_grow,popt_06[0],popt_06[1]),color='violet',label="ajuste no lineal")
                plt.plot(lapsetime06_entero,st_dom_time_06[0])
                plt.xlabel('tiempo de viaje [seg]')
                plt.legend()
                #plt.show()
                plt.savefig(ruta06)
                plt.close()

            name_exp07 = '{}_ajuste_exp_{}_Hz.png'.format(fecha,f07)
            ruta07 = os.path.join(path_fig_exp,name_exp07)
            lapsetime07_entero = t_evento + t_07
            if p1_07 < 0:
                plt.figure()
                plt.title(inicio)
                plt.plot(lapsetime07,Alpt07,color='black',label='envolvente corregida fc={}Hz'.format(f07))                
                plt.plot(lpt07_grow,model_exp(lpt07_grow,popt_07[0],popt_07[1]),color='violet',label="ajuste no lineal")
                plt.plot(lapsetime07_entero,st_dom_time_07[0])
                plt.xlabel('tiempo de viaje [seg]')
                plt.legend()
                #plt.show()
                plt.savefig(ruta07)
                plt.close()

            name_exp08 = '{}_ajuste_exp_{}_Hz.png'.format(fecha,f08)
            ruta08 = os.path.join(path_fig_exp,name_exp08)
            lapsetime08_entero = t_evento + t_08
            if p1_08 < 0:
                plt.figure()
                plt.title(inicio)
                plt.plot(lapsetime08,Alpt08,color='black',label='envolvente corregida fc={}Hz'.format(f08))                
                plt.plot(lpt08_grow,model_exp(lpt08_grow,popt_08[0],popt_08[1]),color='violet',label="ajuste no lineal")
                plt.plot(lapsetime08_entero,st_dom_time_08[0])
                plt.xlabel('tiempo de viaje [seg]')
                plt.legend()
                #plt.show()
                plt.savefig(ruta08)
                plt.close()'''
