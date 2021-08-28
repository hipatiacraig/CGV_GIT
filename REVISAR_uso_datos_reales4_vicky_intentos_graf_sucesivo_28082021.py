# -*- coding: utf-8 -*-
"""
Autora: María Celeste Novak Merquel
"""

import pandas as pd
import numpy as np
from obspy.core import read
from obspy import UTCDateTime
import os
import matplotlib.pyplot as plt
from matplotlib.pyplot import ion, draw, show

# Leo los tags
# Primer fila del archivo: fecha,evento,duracion,amplitud --> indica el orden
# de las columnas. Las comas delimitan las columnas
tags_file = 'd:\Victoria\Desktop\PHD\SCRIPTS\自分の\TRASH\catalogo_LP_2020-Copy.csv'

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

dir_base='d:\Victoria\Desktop\PHD\FUEGO_INSIVUMEH\DATOS_05052021\FG16'
'''
armo el listado de archivos de datos en la base
'''
lista_nombres=os.listdir(dir_base)

# De dates extraigo año y día juliano para convertirlos en cadenas de caracteres
# y usarlo para completar el nombre del wave_file



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

# Se recorre de inicio a fin cada evento y se lo distingue según su tipo
    #for i in range(6854,len(dates)):
    wave_file = os.path.join(dir_base,'GI.FG16.00.BHZ.D.')
    inicio=inicio-5
    diaj = newdates.str[1]
    if len(diaj)==1:
        dayofyear='00'+diaj
    elif len(diaj)==2:
        dayofyear='0'+diaj
    else:
        dayofyear=diaj
    wave_file +=newdates.str[0]+'.'+dayofyear
   # print(wave_file)
    for h in range(len(lista_nombres)):
       # print(h)
        lista_nombres2 = os.path.join(dir_base,lista_nombres[h])
       # print(lista_nombres2)
        if wave_file != lista_nombres2:
            continue
        else:
            wave=read(wave_file)
            evento=wave.slice(inicio,inicio+ 5 +duracion[i])
            print('encontrado'+lista_nombres2)
            traza_evento=evento[0]
            dato=traza_evento.data
            plt.ion() # hace que el gráfico sea interactivo en teoría
            evento.plot()
            # para que el grafico se muestre de inmediato veo que hay veces que agregar un plt.draw(),
            #no parece hacernos falta
            
            #En teoría aclar block = True o block = False en show() me permite manejar 
            #si el código continua corriendo o no dsp que me mostró el gráfico o si espera que yo lo cierre
            #de forma manual para continuar...no lo logro hacer andar
            plt.show(block=True)
            answer = input('Cargar otro slice? ')
            if answer == 'y':
                print('Joya!')
                plt.close("all")
            else:
                print('ok, listo...')
                break
      


