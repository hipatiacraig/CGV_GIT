#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 24 21:49:21 2021

@author: gaby
"""

# -*- coding: utf-8 -*-


import pandas as pd
from obspy.core import read
from obspy import UTCDateTime
"""
Aqui cambiar el nombre del archivo
"""
tags_file='/mnt/hgfs/Compartida_Ubuntu/DATOS/Trazas_Fuego/TagsFuego.csv'
tags=pd.read_csv(tags_file,sep=',',header=0)
"""
Ojo que el csv debe tener en primera fila los titulos de columnas:"Fecha,SCNL,evento" 
para que el read_csv  pueda cargar una dataframe con las tres columnas ya tituladas
""" 
dates=tags.Fecha
tipos=tags.Evento


"""
# Si se va a entrar una unica onda... aqui cambiar el nombre del archivo de onda 
"""
#wave_file='/mnt/hgfs/Compartida_Ubuntu/DATOS/Trazas_Fuego/GI.FG16.00.BHZ.D.2020.214'


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


"""
# crea un vector con el rango de tags/2 para recorrer inicio y fin de cada evento
# asume que los tags estan ordenados con una linea para inicio con tipo de evento y
# la siguiente con el final etiquetado como noise
"""

for i in range(len(dates)):
    
    if tipos[i] != 'Noise':
        wave_file='/mnt/hgfs/Compartida_Ubuntu/DATOS/Trazas_Fuego/GI.FG16.00.BHZ.D'
        ascii_file='/mnt/hgfs/Compartida_Ubuntu/DATOS/Trazas_Fuego/LAKIY/'
        inicio=UTCDateTime(dates[i])

        if tipos[i+1] =='Noise':
            fin=UTCDateTime(dates[i+1])
            if tipos[i]=='Explosion':
                EventCode='EX'
            elif tipos[i]=='LP':
                EventCode='LP'
        else:
            if tipos[i] =='Explosion':
                fin=inicio+80
                EventCode='EX'
            elif tipos[i] =='LP':
                fin=inicio+100
                EventCode='LP'
    else:
        continue
                
    wave_file +='.'+DFdates.anio[i]+'.'+DFdates.jday[i]
    wave=read(wave_file)
    evento=wave.slice(inicio,fin)
    evento.plot()
    traza_evento=evento[0]
    dato=traza_evento.data
    newdato={'COPZ':dato}
    DFdato=pd.DataFrame(newdato)
    DFdato['COPZ']=DFdato['COPZ'].apply(str)
    ascii_file +=DFdates.anio[i]+DFdates.jday[i]+DFdates.hour[i]+DFdates.minute[i]+'.FG16Z'+EventCode
    DFdato.to_csv(ascii_file,sep=' ',index=False)
   
#    input('Press any key to continue...')