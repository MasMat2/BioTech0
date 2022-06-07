import pandas as pd
from numpy import nan
from Funciones import fold_change, Promedio
from Analisis_GSE179277 import gene_name_codigo, gene_name_descriptivo

#IMPORTAMOS LOS GENES A FILTRAR

Genes = pd.read_csv('C:/Users/cogge/Desktop/Investigación/Genes/genes.csv', header = None)
Genes2 = pd.read_csv('C:/Users/cogge/Desktop/Investigación/Genes/genes-metabolismo.csv', header = None)
Genes3 = pd.read_excel('C:/Users/cogge/Desktop/Investigación/Genes/genes2Nuevo.xlsx',header = None)

#IMPORTAMOS LAS BASES DE DATOS

NBD = 'GSM5585685,GSM5585686,GSM5585687,GSM5585688,GSM5585689,GSM5585690,GSM5585691,GSM5585692,GSM5585693,GSM5585694,GSM5585695,GSM5585696'.split(',')

Data0 = pd.DataFrame()

for _ in NBD:
    D = pd.read_excel('C:/Users/cogge/Desktop/Investigación/GSE184390/GSE184390_RAW/'+_+'.xlsx')
    D.index = list(map(lambda string : string[:string.rfind('.')],list(D['gene_id'])))
    Data0 = pd.concat([Data0,D.loc[:,'FPKM'].to_frame()],axis = 1)
    Data0.rename(columns = {'FPKM':_},inplace = True)

filas2 = pd.read_excel('C:/Users/cogge/Desktop/Investigación/Genes_descripcion_codigo.xlsx', header = None)
filas2.columns = ['Codigo','Descripcion']
filas2.index = list(filas2['Codigo'])
filas2.drop(columns = ['Codigo'],inplace = True)

#CREAMOS UNA COPIA DE LA BASE DE DATOS PARA PODER TRABJAR CON ELLA

Data = Data0.copy()

#NORMALIZAMOS

for GAPDH,Columna in zip(list(Data.loc['ENSG00000111640']),list(Data.columns)):
    Data[Columna] = Data[Columna]/GAPDH

#SEPARAMOS LOS GRUPOS

MOI1 = Data.iloc[:,0:4].copy()

MOI10 = Data.iloc[:,4:8].copy()

Controles = Data.iloc[:,-4:].copy()

#SACAMOS LOS PROMEDIOS

Promedio(MOI1,'MOI1')
Promedio(MOI10,'MOI10')
Promedio(Controles,'controles')

#fold-change

MOI1['fold-change MOI1'] = fold_change(Controles['Promedio controles'],MOI1['Promedio MOI1'])
MOI10['fold-change MOI10'] = fold_change(Controles['Promedio controles'],MOI10['Promedio MOI10'])

#UNIMOS LA INFORMACION OBTENIDA

Data = pd.concat([Data,MOI1.loc[:,['Promedio MOI1','fold-change MOI1']]],axis = 1)
Data = pd.concat([Data,MOI10.loc[:,['Promedio MOI10','fold-change MOI10']]],axis = 1)
Data = pd.concat([Data,filas2], axis = 1)

#EXPORTAMOS LA INFORMACION

Genes = list(Genes[0])
Genes2 = list(Genes2[0])
Genes3 = list(Genes3[0])

#VAMOS A ELIMINAR AQUELLOS QUE NO ENCONTRAMOS SU NOMBRE DESCRITO

Data['Descripcion'] = Data['Descripcion'].replace(nan,'_')
Data = Data[Data['Descripcion'] != '_'].copy()
Data.index = list(Data['Descripcion'])
Data.drop(columns = ['Descripcion'],inplace = True)

Data.loc[Genes,['fold-change MOI1','fold-change MOI10']].to_csv('fold_change_genes_1_GSE184390.csv')
Data.loc[Genes2,['fold-change MOI1','fold-change MOI10']].to_csv('fold_change_genes_2_GSE184390.csv')

#Arreglo para Genes 3

Arreglo0 = pd.DataFrame({_:[nan] for _ in list(Data.columns)},index = ['CR', 'DHNTP'])

Data = pd.concat([Data,Arreglo0])

Data.loc[Genes3,['fold-change MOI1','fold-change MOI10']].to_csv('fold_change_genes_3_GSE184390.csv')

Data.to_csv('Resultados_Analisis_GSE184390.csv')