#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 13:24:07 2020

@author: Soma
"""

#Import Libraries
import numpy as np
import gudhi as gd
import glob
import csv 
import pathlib

#Saving the path of this file location
#Need to keep track of where to get the input data (npy files)
#And Need to keep track of where to output data (csv files)
path = pathlib.Path('__file__').parent.absolute()

print(path)



#The python script reads in the cubical image chunks that
#center on the tumor region of interest,
#Computes the cubical complex persistence homology
#ANd outputs the Results

#Only run this code after you have generated the cubical image chunks from 
#the
    
#Creating function to read in an image (array) and name
#And Spitting out a csv cubical complex
def cube_hom(name, array):
    cubical = gd.CubicalComplex(dimensions = array.shape,
                        top_dimensional_cells =  array.flatten())
    
    #Calculate the persistent homology of the filtered cubical complex
    phom = cubical.persistence()

    #Persistence homology per dimension
    phom_0 = cubical.persistence_intervals_in_dimension(0)
    phom_1 = cubical.persistence_intervals_in_dimension(1)
    phom_2 = cubical.persistence_intervals_in_dimension(2)

    #Creating a third column listing dim number
    if len(phom_0) != 0:
        phom_0_form = (np.append(phom_0,np.zeros([len(phom_0),1]),1) +
              np.array([0, 0, 0]))
    else:
        phom_0_form = phom_0
    
    if len(phom_1) != 0:
        phom_1_form = (np.append(phom_1,np.zeros([len(phom_1),1]),1) +
              np.array([0, 0, 1]))
    else:
        phom_1_form = phom_1

    if len(phom_2) != 0:
        phom_2_form = (np.append(phom_2,np.zeros([len(phom_2),1]),1) +
              np.array([0, 0, 2]))
    else :
        phom_2_form = phom_2
    
    file_name = str(path) + '/python_hom/' + name + ".csv"

    # opening the csv file in 'w+' mode 
    file = open(file_name, 'w+', newline ='') 

    # writing the data into the file 
    with file:     
        write = csv.writer(file) 
        write.writerows(phom_0_form) 
        write.writerows(phom_1_form)
        write.writerows(phom_2_form)
  


####Iterating over numpy vars dictionary####
#Reading in all the npy files and storing in dictionary
#Modify string to represent patient ID as dict key
numpy_vars = {}
for np_name in glob.glob(str(path) + '/Numpy_Arrays/*.np[yz]'):
    numpy_vars[np_name.split("/")[-1]] = np.load(np_name)
    
#ITerating function defined above over dictionary items   
for key, value in numpy_vars.items():
    cube_hom(key, value)



 
    