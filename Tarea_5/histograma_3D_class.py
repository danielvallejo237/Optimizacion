# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 13:27:12 2015
@author: odin
"""

import argparse
import numpy as np
import cv2

#Variables Globales
#------------------#
num_bins = 8
num_clases = 2
status_clase = 1
status_lb = 0
#------------------#

#Funcion para dibujar un circulo en el punto
def dibujar_punto(x, y):
    global status_clase
    color = [(255,255,255),(0,0,255),(255,0,0),(0,255,0),(0,255,255)]
    row = matriz_clases_inicial.shape[0]
    col = matriz_clases_inicial.shape[1]
    
    
    #Validar indices
    x = x if x>=0 else 0
    x = x if x< col else col-1
    y = y if y>=0 else 0
    y = y if y< row else row-1
     
    #Guardar clase en la vecindad
    for i in range(-7,8):
        for j in range(-7,8):
            ind_x = x + i
            ind_y = y + j
            
            if ind_x < 0 or ind_x >= col:
                continue
            if ind_y < 0 or ind_y >= row:
                continue
            matriz_clases_inicial[ind_y,ind_x] = status_clase
            

    cv2.circle(img_1 , (x,y), 10, color[status_clase], -1)
    cv2.imshow('Imagen_1', img_1)



#Funcion para recuperar los clic
def onMouse(event, x, y, flags, param):   

    global puntos_1
    global num_clases
    global status_clase
    global status_lb
       

    #Clic IZQ-DOWN del mouse
    if event == cv2.EVENT_LBUTTONDOWN:
        status_lb = 1
        dibujar_punto(x,y)
        
    #Clic IZQ-UP del mouse
    if event == cv2.EVENT_LBUTTONUP:
        status_lb = 0

    #Clic DER-UP del mouse
    if event == cv2.EVENT_RBUTTONDOWN:
        status_clase = status_clase + 1 if (status_clase + 1)<= num_clases else 1 
        print( 'Activa la clase: ', status_clase)
    
    #Movimiento del mouse
    if event == cv2.EVENT_MOUSEMOVE:
        if status_lb==1:
            dibujar_punto(x,y)


#Funcion para calcular los histograms de cada clase 
def calcular_histogramas_OPTI(img, num_bins, clases, num_clases):
    
    #Dimensiones de la imagen
    row = img.shape[0]
    col = img.shape[1]
    
    #Inicializar histogramas a cero
    hist =np.zeros((num_clases,num_bins,num_bins,num_bins,))
    
    #Recorrer cada pixel de la imagen y acumular en la casilla del histograma correspondiente
    for i in range(row):
        for j in range(col):
            
            #Verificar si este pixel fue "seleccionado".
            if clases[i,j]==0:
                continue
            
            #Calcular la posicion de este pixel en los bins del histograma
            clase_id = clases[i,j] - 1
            bin_x = int(float(img[i,j,0])/256.0*num_bins)
            bin_y = int(float(img[i,j,1])/256.0*num_bins)
            bin_z = int(float(img[i,j,2])/256.0*num_bins)
            
            #Acumular en el histograma
            hist[clase_id, bin_x, bin_y, bin_z] = hist[clase_id, bin_x, bin_y, bin_z] + 1

    return hist



#Funcion para guardar histograma (se guarda como un vector)
def guardar_histograma(hist,FileName):
    
    #Abrir archivo
    w = open(FileName, "w")

    #Escribir dimensiones
    w.write('{} {} {}\n'.format(hist.shape[0], hist.shape[1], hist.shape[2]))

    #Escribir entradas de Histograma_3D. Se guardara
    for i in range(hist.shape[0]):
        for j in range(hist.shape[1]):
            for k in range(hist.shape[2]):
                w.write("%d\n" % hist[i,j,k])
        
    #Cerrrar archivo
    w.close()


#Funcion main
if __name__ == '__main__':  
  
    #PARAMETROS DE ENTRADA: Path de la imagen a segmentar y numero de bins
    #EJEMPLO EJECUCION: python histograma_3D_class.py --image rose.png --bins 3
  
    ap = argparse.ArgumentParser()
    ap.add_argument('--image', required = True, help = 'Path de la imagen 1')
    ap.add_argument('--bins', required = True, type = int, help = 'Numero de bins')
    args = vars(ap.parse_args())
    
    #Cargar imagenes y convertirlas a escala de grises
    imagen_1 = cv2.imread(args['image'])
    img_1 = imagen_1.copy()
    row = img_1.shape[0]
    col = img_1.shape[1]
    
    #Determinar el numero de bins
    num_bins = args['bins']
   
    #Alojar memoria para definir las clases
    matriz_clases_inicial = np.zeros((row,col), dtype='int32')
     
    #Nombrar ventana
    cv2.namedWindow('Imagen_1',cv2.WINDOW_NORMAL)

    #Establecer eventos para el mouse
    cv2.setMouseCallback('Imagen_1', onMouse)
    
    #Mostrar imagenes para que el usuario selecciones los 8-pares de puntos
    print( 'Activa la clase: ', status_clase)
    cv2.imshow('Imagen_1', img_1)    
    #Para cerrar la imagen basta con oprimir la tecla 'Esc'
    cv2.waitKey(0)
    
    #Guardar clases dadas por el usuaario
    cv2.imwrite('Strokes.png',img_1)
    
    #Calcular el histograma por cada clase
    H = calcular_histogramas_OPTI(imagen_1, num_bins, matriz_clases_inicial, num_clases)
    
    #Guardar cada histograma
    for i in range(num_clases):
        guardar_histograma(H[i,:,:],'H_'+str(i)+'.txt')

        