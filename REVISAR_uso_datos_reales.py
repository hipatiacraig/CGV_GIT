# -*- coding: utf-8 -*-
"""
Autora: María Celeste Novak Merquel
"""

import pandas as pd
from obspy.core import read
from obspy import UTCDateTime

# Leo los tags
# Primer fila del archivo: fecha,evento,duracion,amplitud --> indica el orden
# de las columnas. Las comas delimitan las columnas
tags_file = 'C:/Users/Propietario/Desktop/tesis_de_grado/ATENUACION/DATOS/catalog/fecha_tipo_dur_ampl_max_LP_2020.csv'

# Genero una variable tipo DataFrame con los datos de tags_file, le indico que 
# las comas separan columnas y que en la fila=0 está el header
tags = pd.read_csv(tags_file,sep=',',header=0)


# Extraigo vectores de datos y los pongo en variables
dates = tags.fecha
tipos = tags.evento
duracion = tags.duracion

# Leo un archivo de señal
#wave_file = 'C:/Users/Propietario/Desktop/tesis_de_grado/ATENUACION/DATOS/FG16/GI.FG16.00.BHE.D.2020.062'


# A partir de los tags, se genera el nombre de la forma de onda a partir de una
# única estación (FG16). La raiz del nombre de archivo debe cambiarse de 
# acuerdo a la estructura de carpetas en la maquina donde esten los datos. P. ej:
# wave_file='/mnt/hgfs/Compartida_Ubuntu/DATOS/Trazas_Fuego/GI.FG16.00.BHZ.D.'

# De dates extraigo año y día juliano para convertirlos en cadenas de caracteres
# y usarlo para completar el nombre del wave_file

# Genero vectores vacíos que después en el for van a ser completados
anio = []
dia = []
hora = []
minuto = []
segundo = []

for i in range(0,len(dates),1):
    inicio = UTCDateTime(dates[i])
    anio.append(inicio.year)
    dia.append(inicio.julday)
    hora.append(inicio.hour)
    minuto.append(inicio.minute)
    segundo.append(inicio.second)

# Se genera un diccionario con toda la información    
newdates = {'anio':anio,'jday':dia,'hour':hora,'minute':minuto,'second':segundo}
DFdates = pd.DataFrame(newdates)
DFdates['anio'] = DFdates['anio'].apply(str)
DFdates['jday'] = DFdates['jday'].apply(str)
DFdates['hour'] = DFdates['hour'].apply(str)
DFdates['minute'] = DFdates['minute'].apply(str)
DFdates['second'] = DFdates['second'].apply(str)

nombre_csv = 'C:/Users/Propietario/Desktop/tesis_de_grado/ATENUACION/DATOS/FG16/lista_nombre_archivos.csv'
nombre = pd.read_csv(nombre_csv,header=0)
lista_nombres = nombre.nombre


# Se recorre de inicio a fin cada evento y se lo distingue según su tipo
for i in range(6854,len(dates)):
    wave_file = 'C:/Users/Propietario/Desktop/tesis_de_grado/ATENUACION/DATOS/FG16/GI.FG16.00.BHZ.D'
    inicio = UTCDateTime(dates[i])
    diaj = DFdates.jday[i]
    
    if len(diaj)==1:
        
        dayofyear='00'+diaj
    elif len(diaj)==2:
        dayofyear='0'+diaj
    else:
        dayofyear=diaj

    wave_file +='.'+DFdates.anio[i]+'.'+dayofyear

    
    for h in range(len(lista_nombres)):
        lista_nombres2 = 'C:/Users/Propietario/Desktop/tesis_de_grado/ATENUACION/DATOS/FG16/' + lista_nombres[h]
        if wave_file == lista_nombres2:
            wave=read(wave_file)
            evento=wave.slice(inicio,inicio+duracion[i])
            evento.plot()
            traza_evento=evento[0]
            dato=traza_evento.data
            print(lista_nombres2)
        else:
            continue
        

