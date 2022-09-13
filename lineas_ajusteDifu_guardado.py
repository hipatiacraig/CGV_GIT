#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 11:47:43 2022

@author: gaby
"""
import csv
def model_dif(t, a1,a2,a3):
    return  a1-a2*t-a3/t

p0_01_03 = (np.exp(p2_01),-p1_01[0],-1000) # 
p0_01_04 = (np.exp(p2_01),-p1_01[0],0) # 
popt_01_03, pcov_01_03 = curve_fit(model_dif, lapsetime01, np.log(pAlpt01),p0_01_03)
popt_01_04, pcov_01_04 = curve_fit(model_dif, lapseleto, np.log(pAlpt0leto),p0_01_04)
Qi=2*np.pi*(w1+w2)/(2*popt_01_03[1])
dist=r_fg16
d=dist**2/(4*popt_01_03[2])
p=2 #div geom O.S p=2, O.I p=3
Qs=np.pi*(w1+w2)*p*d/vel**2

gaby=100
while len(lapsetime01)>30*50:
    gaby=gaby-50
    print(gaby)
    if gaby < 0:
        break
    elif gaby>0:
        print(gaby)
        lapsetime01=lapsetime01[:len(lapsetime01)-500]
print(len(lapsetime01))       

#######################################################################################
# ojo que parte de esto ya esta en el programa
dire_salida='/media/gaby/Backup Plus/SISMOLOGIA_DATOS_Y_DOCS/Vn_FUEGO/resultados'

existencias=os.listdir(dire_salida)
scnl1 = scnl[0]
red = scnl1[9:11]
archivo_salida = red + '_Q.csv'
path_archivo_salida = dire_salida +'/' + archivo_salida

### si ya existe un archivo de salida (.csv) creado por el programa, éste lo detectará
# y avisará por terminal que ya existe. En caso contrario, avisa que no 
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

##### Se puede definir las frec medias si sirve para simplificar, pero eso segun te parezca
Fc=(w1+w2)/2
Fc=np.append(Fc,(w3+w4)/2)
Fc=np.append(Fc,(w5+w6)/2)
Fc=np.append(Fc,(w7+w8)/2)
Fc=np.append(Fc,(w9+w10)/2)
Fc=np.append(Fc,(w11+w12)/2)
Fc=np.append(Fc,(w13+w14)/2)
Fc=np.append(Fc,(w15+w16)/2)
Fc
Out[90]: array([ 1. ,  1.5,  2. ,  2.5,  3. ,  3.5,  4. ,  4.5])

## ########## Una vez calculados todos los valores para una traza
############# Definimos un diccionario para cada frec e imprimimos una linea por vez
### es decit que la siguiente linea hay que repetirla para todas las frecuencias
field_names = ['UTCDateTime','Fc','lapsetime','dist','p_dg','vel','a_lin','b_lin','a_exp','b_exp','Qc_SBS','a1_dif','a2_dif','a3_dif','Qi','d','Qs']
##  
# Dictionary
lpt=lapsetime01[-1:] # asi en cada lapse time, nos quedamos con el ultimo elemento como float para imprimirlo
lptfin=lpt[0]
dict = {'UTCDateTime':inicio,'Fc':Fc[0],'lapsetime':lptfin,'dist':dist,'p_dg':p,'vel':vel,'a_lin':np.exp(p2_01),'b_lin':-p1_01[0],'a_exp':popt_01[0],'b_exp':popt_01[1],'Qc_SBS': Q_01,'a1_dif':popt_01_03[0],'a2_dif':popt_01_03[1],'a3_dif':popt_01_03[2],'Qi':Qi,'d':d,'Qs':Qs}

with open(path_archivo_salida, 'a') as csv_file:
    dict_object = csv.DictWriter(csv_file, fieldnames=field_names) 
  
    dict_object.writerow(dict)






ETRPYdaydf['tiempos']=ETRPYdaydf['tiempos'].apply(str)
ETRPYdaydf['tsallis']=ETRPYdaydf['tsallis'].apply(str)
ETRPYdaydf['renyi']=ETRPYdaydf['renyi'].apply(str)
#ETRPYdaydf['complexity']=ETRPYdaydf['complexity'].apply(str)
ETRPYdaydf['permutation']=ETRPYdaydf['permutation'].apply(str)
ETRPYdaydf['entropy']=ETRPYdaydf['entropy'].apply(str)

ETRPYdaydf.to_csv(path_oldoutfileETday,sep=' ',index=False)   
