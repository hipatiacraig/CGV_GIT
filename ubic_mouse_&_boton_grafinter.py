# -*- coding: utf-8 -*-
"""
Created on Mon May 31 09:36:56 2021

@author: Victoria
"""

"""
Fuente: https://matplotlib.org/stable/users/event_handling.html
Para recibir eventos, hay que escribir una función tipo callback
y luego conectar mi función con el "event manager", que es parte de FigureCanvasBase. 
(https://matplotlib.org/stable/api/backend_bases_api.html#matplotlib.backend_bases.FigureCanvasBase)
Este ejemplo muestra, en teoría, la ubicación del click de mouse y qué btón fue apretado
"""
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots()
ax.plot(np.random.rand(10))

def onclick(event):
    print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          ('double' if event.dblclick else 'single', event.button,
           event.x, event.y, event.xdata, event.ydata))

cid = fig.canvas.mpl_connect('button_press_event', onclick)
