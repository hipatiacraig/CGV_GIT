#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 2021 07 13

@author: gaby
"""

# -*- coding: utf-8 -*-


import pandas as pd
#from obspy.signal.freqattributes import bandwidth
from obspy.core import read
from obspy import UTCDateTime
import os
#import statistics
"""
Aqui cambiar el nombre del archivo
"""
tags_dir='/media/gaby/Backup Plus/SISMOLOGIA DATOS Y DOCS/Vn_FUEGO/tags/'
tags_file=tags_dir + 'fecha_tipo_dur_ampl_max_LP_2020.csv'
EventCode='FLP'
tags=pd.read_csv(tags_file,sep=',',header=0)
"""
Ojo que el csv debe tener en primera fila los titulos de columnas:"Fecha,SCNL,evento" 
para que el read_csv  pueda cargar una dataframe con las tres columnas ya tituladas
""" 
dates=tags.fecha
events=tags.tipo
dur=tags.dur
ampl=tags.ampl_max



"""
# Si se va a entrar una unica onda... aqui cambiar el nombre del archivo de onda 
"""
#wave_file='/mnt/hgfs/Compartida_Ubuntu/DATOS/Trazas_Fuego/GI.FG16.00.BHZ.D.2020.214'


datos_dir='/media/gaby/Backup Plus/SISMOLOGIA DATOS Y DOCS/Vn_FUEGO/FG16'
'''
armo el listado de archivos de datos en la base
'''
lista_nombres=os.listdir(datos_dir)




"""
En este caso va a generar a partir de los tags el nombre de la forma de onda
 supuesta una unica estacion (FG16). La raiz del nobre de archivo debe cambiarse de 
 acuerdo a la estructura de carpetas en la maquina donde esten los datos. P. ej:
 wave_file='/mnt/hgfs/Compartida_Ubuntu/DATOS/Trazas_Fuego/GI.FG16.00.BHZ.D.'
"""


"""
Extraigo anio y dia juliano de dates para convertirlos en cadenas de caracteres y
usarlos para completar el nombre del wave_file
"""
anio=[]
dia=[]
hora=[]
minuto=[]
segundo=[]

for i in range(0,len(dates),1):
    inicio=UTCDateTime(dates[i])
    anio.append(inicio.year)
    dia.append(inicio.julday)
    hora.append(inicio.hour)
    minuto.append(inicio.minute)
    segundo.append(inicio.second)
 
newdates={'anio':anio,'jday':dia,'hour':hora,'minute':minuto,'second':segundo}
DFdates=pd.DataFrame(newdates)
DFdates['anio']=DFdates['anio'].apply(str)
DFdates['jday']=DFdates['jday'].apply(str)
DFdates['hour']=DFdates['hour'].apply(str)
DFdates['minute']=DFdates['minute'].apply(str)
DFdates['second']=DFdates['second'].apply(str)




for i in range(0,len(dates),1):
    inicio=UTCDateTime(dates[i])-5
    fin=inicio+dur[i]+5
    wave_file='/media/gaby/Backup Plus/SISMOLOGIA DATOS Y DOCS/Vn_FUEGO/FG16/GI.FG16.00.BHZ.D.'
    if len(DFdates.jday[i])==1:
        dayofyear_str='00'+DFdates.jday[i]
    elif len(DFdates.jday[i])==2:
        dayofyear_str='0'+DFdates.jday[i]
    else:
        dayofyear_str=DFdates.jday[i]
    wave_file +=DFdates.anio[i]+'.'+dayofyear_str
    wave=read(wave_file)
    evento=wave.slice(inicio,fin)
    #evento.plot()
    
    #continue
       
    traza_evento=evento[0]
    #bandw_signal = bandwidth(traza_evento)
    dato=traza_evento.data
    #bandw_signal = bandwidth(data=traza_evento.data, fs=traza_evento.meta.sampling_rate, smoothie=5, fk=10)
      
    newdato={'COPZ':dato}
    DFdato=pd.DataFrame(newdato)
    DFdato['COPZ']=DFdato['COPZ'].apply(str)
    ascii_file =tags_dir+DFdates.anio[i]+DFdates.jday[i]+DFdates.hour[i]+DFdates.minute[i]+'.FG16Z'+EventCode
    DFdato.to_csv(ascii_file,sep=' ',index=False)
   
#    input('Press any key to continue...')