import pandas as pd
from numpy import nan
import numpy as np
from Funciones import fold_change, Promedio

#IMPORTAMOS LOS GENES A FILTRAR

Genes = pd.read_csv('C:/Users/cogge/Desktop/Investigación/Genes/genes.csv', header = None)
Genes2 = pd.read_csv('C:/Users/cogge/Desktop/Investigación/Genes/genes-metabolismo.csv', header = None)
Genes3 = pd.read_excel('C:/Users/cogge/Desktop/Investigación/Genes/genes2Nuevo.xlsx',header = None)

#BASE DE DATOS

Data = pd.read_csv('C:/Users/cogge/Desktop/Investigación/GSE185557/GSE185557_count_matrix.csv/GSE185557_count_matrix.csv')
Data.index = list(Data['ID_REF'])
Data.drop(columns = ['ID_REF'],inplace = True)

Traduccion = pd.read_csv('C:/Users/cogge/Desktop/Investigación/GSE185557/Traduccion_Genes.csv')

#TOMANDO EN CUENTA QUE DE LOS GENES QUE NO TENIAMOS TRADUCCION ERAN 1551 Y DE ELLOS MENOS DEL 5%
#ENCONTRAMOS TRADUCCION ENTONCES LOS GEN_REF SE ASIGNAN SU RESPECTIVO ID_REF

for _2 in list(_ for _ in range (len(Traduccion.index)) if type(Traduccion.iloc[_,1]) == float):
    Traduccion.iloc[_2,1] = Traduccion.iloc[_2,0]

DB = Data.copy()
DB.index = list(Traduccion['GEN_REF'])

#NORMALIZAMOS RESPECTO AL GEN GAPDH

for GAPDH, Columna in zip(list(DB.loc['GAPDH']), list(DB.columns)):
    DB[Columna] = DB[Columna]/GAPDH
    
#SEPARAMOS LOS GRUPOS

#Cord blood COVID-19

C_B_COVID_19 = DB.loc[:,list(C_B_COVID for C_B_COVID in list(DB.columns) if C_B_COVID.startswith('CB-SA'))]

#Maternal blood COVID-19

MB_COVID_19 = DB.loc[:,list(MB_COVID for MB_COVID in list(DB.columns) if MB_COVID.startswith('MB-SA'))]

#Cord blood Control

C_B_Control = DB.loc[:,list(C_B_C for C_B_C in list(DB.columns) if C_B_C.startswith('CB-C'))]

#Maternal blood Control

MB_Control = DB.loc[:,list(MB_C for MB_C in list(DB.columns) if MB_C.startswith('MB-C'))]

#SACAMOS LOS PROMEDIOS

Promedio(C_B_COVID_19,'CB - COVID')
Promedio(MB_COVID_19,'MB - COVID')
Promedio(C_B_Control,'CB - Control')
Promedio(MB_Control,'MB - Control')

#fold-change

C_B_COVID_19['fold-change CB'] = fold_change(C_B_Control['Promedio CB - Control'], C_B_COVID_19['Promedio CB - COVID'])
MB_COVID_19['fold-change MB'] = fold_change(MB_Control['Promedio MB - Control'],MB_COVID_19['Promedio MB - COVID'])

#UNIMOS LA INFORMACION OBTENIDA

DB = pd.concat([DB,C_B_COVID_19.loc[:,['Promedio CB - COVID','fold-change CB']]],axis = 1)
DB = pd.concat([DB,MB_COVID_19.loc[:,['Promedio MB - COVID','fold-change MB']]],axis = 1)

#EXPORTAMOS LA INFORMACION

#Arreglo para los genes que no se encuentran en la base de datos pero los necesitamos buscar

Arreglo0 = pd.DataFrame({g1:[nan] for g1 in list(DB.columns)},
                        index = ['ISG15', 'IL17A', 'HLA-F-AS1', 'TLR8-AS1', 'SOCS2-AS1', 'SOCS7','GLUD2','CR','DHNTP'])

DB = pd.concat([DB,Arreglo0])

DB.to_csv('Resultados_Analisis_GSE185557.csv')

DB.loc[list(Genes[0]),['fold-change CB','fold-change MB']].to_csv('fold_change_genes_1_GSE185557.csv')
DB.loc[list(Genes2[0]),['fold-change CB','fold-change MB']].to_csv('fold_change_genes_2_GSE185557.csv')
DB.loc[list(Genes3[0]),['fold-change CB','fold-change MB']].to_csv('fold_change_genes_3_GSE185557.csv')