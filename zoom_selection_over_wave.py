# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 14:58:49 2021

@author: Victoria
"""
import numpy as np
import matplotlib.pyplot as plt
from obspy.core import read 

#leo el archivo miniSEED desde mi carpeta de datos
st = read('d:/Victoria/Desktop/PHD/FUEGO_INSIVUMEH/DATOS_05052021/FG8/GI.FG8.00.BHE.D.2020.061')

#convierto mi dato en traza
tr = st[0]

#imprimo metadatos del dato, como la fecha y hora de inicio y fin del dato
print(tr.stats)


#armo figura con señal y etiquetas tipo matplotlib(porque queremos la magia del zoom y esas cosas bellas que tiene),
#usando el "plot" de obspy con sus virtudes que incluyen que entienda
#que tr.times son los tiempos (osea que con pocos comandos puedo traer atributos del metadato de la señal)
#ante dudas ir a https://docs.obspy.org/tutorial/code_snippets/waveform_plotting_tutorial.html

#fig,ax0 = plt.subplots(sharex=True)
fig = plt.figure()
ax0 = fig.add_subplot(1,1,1)
#plt.plot tmb funciona
ax0.plot(tr.times("matplotlib"),tr.data, label='signal')
ax0.xaxis_date()
fig.autofmt_xdate()



def onclick(event):
    print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          ('double' if event.dblclick else 'single', event.button,
           event.x, event.y, event.xdata, event.ydata))
  
    
cid = fig.canvas.mpl_connect('button_press_event', onclick)



#falta asociar el xdata con xdate (osea, tengo el dato para x pero no logro relacionarlo con el date 
#tal que en vez de un numero me devuelva fecha y hora de mi pick)

#se me ocurre que se puede difereciar 1 click (que se hace al hacer zoom por ejemplo) de 2 click para dato seleccionado
#tomar selección de dato "if" arranca la linea de texto con "double" en vez de con "single"...es una opción