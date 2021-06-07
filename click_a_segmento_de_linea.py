# -*- coding: utf-8 -*-
"""
Created on Mon May 31 09:51:54 2021

@author: Victoria
"""

"""
Atributos de evento
este programa se supone (porque no me estaríasaliendo) que tiene que graficar un segmento 
de línea cada vez que hago click.


Todos los eventos de Matplotlib heredan de la clase base matplotlib.backend_bases.Event, 
que almacena los atributos:

name (el nombre del evento)
canva(la instancia de FigureCanvas que genera el evento)
guiEvent(el evento GUI que desencadenó el evento Matplotlib)

Los eventos más comunes que son los eventos de click / soltar de teclas y los eventos de click / soltar 
y movimiento del mouse. 
Las clases KeyEvent y MouseEvent que manejan estos eventos se derivan de LocationEvent, que tiene los siguientes atributos:

x, y(posición x e y del mouse en píxeles desde la izquierda y la parte inferior del canva)
inaxes(la instancia de axisas sobre la que se encuentra el mouse, si existe; más None)
xdata, ydata(posición x e y del mouse en coordenadas de datos, si el mouse está sobre un eje) 

"""

"""
al final para obtener grafico interactivo de matplotlib
https://stackoverflow.com/questions/54840619/how-to-plot-in-a-new-window-in-spyder-idepython
"""
from matplotlib import pyplot as plt

class LineBuilder:
    def __init__(self, line):
        self.line = line
        self.xs = list(line.get_xdata())
        self.ys = list(line.get_ydata())
        self.cid = line.figure.canvas.mpl_connect('button_press_event', self)

    def __call__(self, event):
        print('click', event)
        if event.inaxes!=self.line.axes: return
        self.xs.append(event.xdata)
        self.ys.append(event.ydata)
        self.line.set_data(self.xs, self.ys)
        self.line.figure.canvas.draw()



fig, ax = plt.subplots()
ax.set_title('click to build line segments')
line, = ax.plot([0], [0])  # empty line
linebuilder = LineBuilder(line)

plt.show()
