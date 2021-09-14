# -*- coding: utf-8 -*-
"""
Autora: María Celeste Novak Merquel
"""

import pandas as pd
import numpy as np
from obspy.core import read
from obspy import UTCDateTime
import os
#import matplotlib.pyplot as plt


# Leo los tags
# Primer fila del archivo: fecha,evento,duracion,amplitud --> indica el orden
# de las columnas. Las comas delimitan las columnas
#tags_file = '/media/gaby/Backup Plus/SISMOLOGIA DATOS Y DOCS/Vn_FUEGO/tags/fecha_tipo_dur_ampl_max_LP_2020.csv'
tags_dir='/media/gaby/Backup Plus/SISMOLOGIA DATOS Y DOCS/Vn_FUEGO/tags/'
tags_file=tags_dir + 'fecha_tipo_dur_ampl_max_LP_2020.csv'
#tags_file=tags_dir + 'swarm.csv'
#EventCode='FIS'
EventCode='FLP'

# Genero una variable tipo DataFrame con los datos de tags_file, le indico que 
# las comas separan columnas y que en la fila=0 está el header
tags = pd.read_csv(tags_file,sep=',',header=0)

'''
 Extraigo vectores de datos y los pongo en variables

para LP y TR
'''
dates = tags.fecha
tipos = tags.tipo
duracion = tags.dur

''' 
para Infrasonido
'''
#dates=tags.fecha
# tomo una duracion promedio de los LP
#duracion = 60
'''
para ruido
'''
#duracion = 5

# Leo un archivo de señal
#wave_file = 'C:/Users/Propietario/Desktop/tesis_de_grado/ATENUACION/DATOS/FG16/GI.FG16.00.BHE.D.2020.062'

'''
A partir de los tags, se genera el nombre de la forma de onda a partir de una
única estación (FG16 o FG12). La raiz del nombre de archivo debe cambiarse de 
 acuerdo a la estructura de carpetas en la maquina donde esten los datos. P. ej:
 wave_file='/mnt/hgfs/Compartida_Ubuntu/DATOS/Trazas_Fuego/GI.FG16.00.BHZ.D.'
 Las EX fueron identificadas en FG12 por lo tanto se seleccionan en FG12, lo mismo que los LPs para comparar
'''
#datos_dir='/media/gaby/Backup Plus/SISMOLOGIA DATOS Y DOCS/Vn_FUEGO/FG16'
datos_dir='/media/gaby/Backup Plus/SISMOLOGIA DATOS Y DOCS/Vn_FUEGO/FG12'
'''
armo el listado de archivos de datos en la base
'''
lista_nombres=os.listdir(datos_dir)
'''
Se recorre de inicio a fin cada evento y se lo distingue según su tipo
De dates extraigo año y día juliano para convertirlos en cadenas de caracteres
 y usarlo para completar el nombre del wave_file
'''


for i in range(0,len(dates),1):
    inicio = UTCDateTime(dates[i])
    anio=inicio.year
    dia=inicio.julday
    hora=inicio.hour
    minuto=inicio.minute
    segundo=inicio.second

# Se genera un diccionario con toda la información    
    dates_array=np.array([anio, dia, hora, minuto ,segundo])
    newdates=pd.DataFrame(dates_array,columns=['str'])
    newdates['str'] = newdates['str'].apply(str)    


    #wave_file = os.path.join(datos_dir,'GI.FG16.00.BHZ.D.')
    wave_file = os.path.join(datos_dir,'GI.FG12.00.BHZ.D.')
    preinicio=inicio-5
   #para ruido
    #preinicio=inicio-15
    diaj = newdates.str[1]
    if len(diaj)==1:
        dayofyear='00'+diaj
    elif len(diaj)==2:
        dayofyear='0'+diaj
    else:
        dayofyear=diaj
        
    horastr = newdates.str[2]
    if len(horastr)==1:
        horastr2='0'+horastr
    elif len(horastr)==2:
        horastr2=horastr

    minstr = newdates.str[3]
    if len(minstr)==1:
        minstr2='0'+minstr
    elif len(minstr)==2:
        minstr2=minstr
        
        
    wave_file +=newdates.str[0]+'.'+dayofyear
   # print(wave_file)
    for h in range(len(lista_nombres)):
       # print(h)
        lista_nombres2 = os.path.join(datos_dir,lista_nombres[h])
       # print(lista_nombres2)
        if wave_file != lista_nombres2:
            continue
        else:
            wave=read(wave_file)
            # para los LP y TR que se leen con duracion cada evento
            evento=wave.slice(preinicio,inicio+duracion[i]+5)
            # para IS que NO hay duracion provista para cada evento
            #evento=wave.slice(preinicio,inicio+duracion+5)
            #para ruido
            #evento=wave.slice(preinicio,preinicio+duracion+5)
            print('encontrado'+lista_nombres2)
            traza_evento=evento[0]
            dato=traza_evento.data
            #evento.plot()
                
            newdato={'COPZ':dato}
            DFdato=pd.DataFrame(newdato)
            DFdato['COPZ']=DFdato['COPZ'].apply(str)
            # los LP,TR y RD se etiquetaron en FG16
            #ascii_file =tags_dir+newdates.str[0]+dayofyear+horastr2+minstr2+'.FG16Z'+EventCode
            ascii_file =tags_dir+newdates.str[0]+dayofyear+horastr2+minstr2+'.FG12Z'+EventCode
            # los infrasonidos se etiquetaron en FG12 sensor de infrasonido pero los extraemos en sismometro
            #ascii_file =tags_dir+newdates.str[0]+dayofyear+horastr2+minstr2+'.FG12Z'+EventCode
            DFdato.to_csv(ascii_file,sep=' ',index=False)            
            #input('presione una tecla...')
            break
    #print('ya cargue archivo y plotee')
        

