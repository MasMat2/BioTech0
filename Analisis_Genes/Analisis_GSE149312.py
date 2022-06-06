import pandas as pd
import statistics as sts
from numpy import nan
from Funciones import fold_change

#IMPORTAMOS LOS GENES A FILTRAR

Genes = pd.read_csv('C:/Users/cogge/Desktop/Investigación/Genes/genes.csv', header = None)
Genes2 = pd.read_csv('C:/Users/cogge/Desktop/Investigación/Genes/genes-metabolismo.csv', header = None)
Genes3 = pd.read_excel('C:/Users/cogge/Desktop/Investigación/Genes/genes2Nuevo.xlsx',header = None)

#BASES DE DATOS GSE149312

Data2_corona_intestine_exp1_ndata0 = pd.read_csv('C:/Users/cogge/Desktop/Investigación/GSE149312/GSE149312_corona_intestine_exp1_ndata.csv')
Filanueva0 = list(Data2_corona_intestine_exp1_ndata0['Unnamed: 0'].apply(lambda string : string[:string.find('_')]))
Data2_corona_intestine_exp1_ndata0.index = list(Data2_corona_intestine_exp1_ndata0['Unnamed: 0'])
Data2_corona_intestine_exp1_ndata0.drop(columns = ['Unnamed: 0'],inplace = True)

Data2_corona_intestina_exp2_ndata0 = pd.read_csv('C:/Users/cogge/Desktop/Investigación/GSE149312/GSE149312_corona_intestine_exp2_ndata.csv')
Filanueva1 = list(Data2_corona_intestina_exp2_ndata0['Unnamed: 0'].apply(lambda string : string[:string.find('_')]))
Data2_corona_intestina_exp2_ndata0.index = list(Data2_corona_intestina_exp2_ndata0['Unnamed: 0'])
Data2_corona_intestina_exp2_ndata0.drop(columns = ['Unnamed: 0'],inplace = True)

#CREAMOS UNA COPIA CON LA CUAL VAMOS A TRABAJAR

Data2_corona_intestine_exp1_ndata = Data2_corona_intestine_exp1_ndata0.copy()
Data2_corona_intestina_exp2_ndata = Data2_corona_intestina_exp2_ndata0.copy()

#GENES DE GSE149312

Genes_GSE149312_0 = list(Data2_corona_intestine_exp1_ndata.index)
Genes_GSE149312_1 = list(Data2_corona_intestina_exp2_ndata.index)

#NORMALIZACION DE exp1 y exp2

for GAPDH,Columna in zip(list(Data2_corona_intestine_exp1_ndata.loc['GAPDH__chr12']),list(Data2_corona_intestine_exp1_ndata.columns)) :
    Data2_corona_intestine_exp1_ndata[Columna] = Data2_corona_intestine_exp1_ndata[Columna]/GAPDH
    
for GAPDH0,Columna0 in zip(list(Data2_corona_intestina_exp2_ndata.loc['GAPDH__chr12']),list(Data2_corona_intestina_exp2_ndata.columns)) :
    Data2_corona_intestina_exp2_ndata[Columna0] = Data2_corona_intestina_exp2_ndata[Columna0]/GAPDH0
    
#DIVIDIMOS LOS CONTROLES PARA exp1

Control1_Data2_corona_intestine_exp1_ndata = Data2_corona_intestine_exp1_ndata.loc[:,'1_1':'1b60_2'].copy()
Control2_Data2_corona_intestine_exp1_ndata = Data2_corona_intestine_exp1_ndata.loc[:,'6_1':'6b60_2'].copy()

#PROMEDIOS

Control1_Data2_corona_intestine_exp1_ndata['Promedio Control 1'] = list(sts.mean(Control1_Data2_corona_intestine_exp1_ndata.loc[Control11,['1_1','1_2']]) for Control11 in Genes_GSE149312_0)
Control1_Data2_corona_intestine_exp1_ndata['Promedio 1b24'] = list(sts.mean(Control1_Data2_corona_intestine_exp1_ndata.loc[Control12,['1b24_1','1b24_2']]) for Control12 in Genes_GSE149312_0)
Control1_Data2_corona_intestine_exp1_ndata['Promedio 1b60'] = list(sts.mean(Control1_Data2_corona_intestine_exp1_ndata.loc[Control13,['1b60_1','1b60_2']]) for Control13 in Genes_GSE149312_0)

Control2_Data2_corona_intestine_exp1_ndata['Promedio Control 2'] = list(sts.mean(Control2_Data2_corona_intestine_exp1_ndata.loc[Control21,['6_1','6_2']]) for Control21 in Genes_GSE149312_0)
Control2_Data2_corona_intestine_exp1_ndata['Promedio 6a24'] = list(sts.mean(Control2_Data2_corona_intestine_exp1_ndata.loc[Control22,['6a24_1','6a24_2']]) for Control22 in Genes_GSE149312_0)
Control2_Data2_corona_intestine_exp1_ndata['Promedio 6a60'] = list(sts.mean(Control2_Data2_corona_intestine_exp1_ndata.loc[Control23,['6a60_1','6a60_2']]) for Control23 in Genes_GSE149312_0)
Control2_Data2_corona_intestine_exp1_ndata['Promedio 6b24'] = list(sts.mean(Control2_Data2_corona_intestine_exp1_ndata.loc[Control24,['6b24_1','6b24_2']]) for Control24 in Genes_GSE149312_0)
Control2_Data2_corona_intestine_exp1_ndata['Promedio 6b60'] = list(sts.mean(Control2_Data2_corona_intestine_exp1_ndata.loc[Control25,['6b60_1','6b60_2']]) for Control25 in Genes_GSE149312_0)

Data2_corona_intestina_exp2_ndata['Promedio Control 7'] = list(sts.mean(Data2_corona_intestina_exp2_ndata.loc[Control71,['7_1','7_2']]) for Control71 in Genes_GSE149312_1)
Data2_corona_intestina_exp2_ndata['Promedio 7a72'] = list(sts.mean(Data2_corona_intestina_exp2_ndata.loc[Control72,['7a72_1','7a72_2']]) for Control72 in Genes_GSE149312_1)
Data2_corona_intestina_exp2_ndata['Promedio 7b72'] = list(sts.mean(Data2_corona_intestina_exp2_ndata.loc[Control73,['7b72_1','7b72_2']]) for Control73 in Genes_GSE149312_1)

#fold-change

Control1_Data2_corona_intestine_exp1_ndata['fold-change Control1 1b24'] = fold_change(Control1_Data2_corona_intestine_exp1_ndata['Promedio Control 1'],Control1_Data2_corona_intestine_exp1_ndata['Promedio 1b24'])
Control1_Data2_corona_intestine_exp1_ndata['fold-change Control1 1b60'] = fold_change(Control1_Data2_corona_intestine_exp1_ndata['Promedio Control 1'],Control1_Data2_corona_intestine_exp1_ndata['Promedio 1b60'])

Control2_Data2_corona_intestine_exp1_ndata['fold-change Control2 6a24'] = fold_change(Control2_Data2_corona_intestine_exp1_ndata['Promedio Control 2'],Control2_Data2_corona_intestine_exp1_ndata['Promedio 6a24'])
Control2_Data2_corona_intestine_exp1_ndata['fold-change Control2 6a60'] = fold_change(Control2_Data2_corona_intestine_exp1_ndata['Promedio Control 2'],Control2_Data2_corona_intestine_exp1_ndata['Promedio 6a60'])
Control2_Data2_corona_intestine_exp1_ndata['fold-change Control2 6b24'] = fold_change(Control2_Data2_corona_intestine_exp1_ndata['Promedio Control 2'],Control2_Data2_corona_intestine_exp1_ndata['Promedio 6b24'])
Control2_Data2_corona_intestine_exp1_ndata['fold-change Control2 6b60'] = fold_change(Control2_Data2_corona_intestine_exp1_ndata['Promedio Control 2'],Control2_Data2_corona_intestine_exp1_ndata['Promedio 6b60'])

Data2_corona_intestina_exp2_ndata['fold-change Control7 7a72'] = fold_change(Data2_corona_intestina_exp2_ndata['Promedio Control 7'],Data2_corona_intestina_exp2_ndata['Promedio 7a72'])
Data2_corona_intestina_exp2_ndata['fold-change Control7 7b72'] = fold_change(Data2_corona_intestina_exp2_ndata['Promedio Control 7'],Data2_corona_intestina_exp2_ndata['Promedio 7b72'])

## ANEXAMOS LOS CONTROLES A exp1

Data2_corona_intestine_exp1_ndata = pd.concat([Data2_corona_intestine_exp1_ndata,Control1_Data2_corona_intestine_exp1_ndata.loc[:,['Promedio Control 1','Promedio 1b24','Promedio 1b60','fold-change Control1 1b24','fold-change Control1 1b60']]], axis = 1)
Data2_corona_intestine_exp1_ndata = pd.concat([Data2_corona_intestine_exp1_ndata,Control2_Data2_corona_intestine_exp1_ndata.loc[:,['Promedio Control 2','Promedio 6a24','Promedio 6a60','Promedio 6b24','Promedio 6b60','fold-change Control2 6a24','fold-change Control2 6a60','fold-change Control2 6b24','fold-change Control2 6b60']]], axis = 1)

#MODIFICAMOS LOS NOMBRES DE LAS FILAS

Data2_corona_intestine_exp1_ndata.index = Filanueva0
Data2_corona_intestina_exp2_ndata.index = Filanueva1

#DF fold-change

fold_change_data_2_1 = Data2_corona_intestine_exp1_ndata.loc[:,['fold-change Control1 1b24','fold-change Control1 1b60','fold-change Control2 6a24','fold-change Control2 6a60','fold-change Control2 6b24','fold-change Control2 6b60']]
fold_change_data_2_2 = Data2_corona_intestina_exp2_ndata.loc[:,['fold-change Control7 7a72','fold-change Control7 7b72']]

#METEMOS LOS GENES QUE NO EXISTEN EN GSE149312

Arreglo1 = pd.DataFrame({_:[nan] for _ in list(fold_change_data_2_1.columns)}, index = ['IFITM5', 'IFNA1', 'IFNA2', 'IRF4', 'TRIM17', 'SERPINA5', 'KNG1', 'CSF1R', 'CSF2RA', 'CSF2RB', 'IL17A', 'IL6', 'TLR8', 'STAT4', 'CXCL12', 'IGF2', 'CRP', 'ESR1', 'P2RY6'])
Arreglo2 = pd.DataFrame({_:[nan] for _ in list(fold_change_data_2_1.columns)}, index = ['DBH', 'HPD', 'PAH', 'TPH1', 'TPH2', 'AMPD1', 'IBSP', 'THY1', 'TYR'])
Arreglo5 = pd.DataFrame({_:[nan] for _ in list(fold_change_data_2_1.columns)}, index = ['CR', 'DHNTP', 'APBA2', 'TDO2', 'TG', 'HTR1A', 'HTR6'])
Arreglo3 = pd.DataFrame({_:[nan] for _ in list(fold_change_data_2_2.columns)}, index = ['IFITM5', 'IFNA1', 'IFNA2', 'IRF4', 'F2', 'GJA1', 'PLG', 'SERPIND1', 'VEGFC', 'COL5A1', 'CSF1R', 'CSF2RA', 'IL17A', 'IL6', 'MSN', 'HLA-DQA2', 'IGF1', 'AGTR1', 'AGTR2', 'MRC1', 'TLR1', 'TLR8', 'TLR9', 'IGF2', 'CRP', 'ESR1', 'P2RY6'])
Arreglo4 = pd.DataFrame({_:[nan] for _ in list(fold_change_data_2_2.columns)}, index = ['DBH', 'HPD', 'PAH', 'TPH1', 'TPH2', 'ALPL', 'AMPD1', 'IBSP', 'THY1', 'TYR']) 
Arreglo6 = pd.DataFrame({_:[nan] for _ in list(fold_change_data_2_2.columns)}, index = ['CR', 'DHNTP', 'APBA2', 'FOXP2', 'TDO2', 'TG', 'HTR2C', 'HTR1A', 'HTR6'])

fold_change_data_2_1 = pd.concat([fold_change_data_2_1,Arreglo1])
fold_change_data_2_1 = pd.concat([fold_change_data_2_1,Arreglo2])
fold_change_data_2_1 = pd.concat([fold_change_data_2_1,Arreglo5])

fold_change_data_2_2 = pd.concat([fold_change_data_2_2,Arreglo3])
fold_change_data_2_2 = pd.concat([fold_change_data_2_2,Arreglo4])
fold_change_data_2_2 = pd.concat([fold_change_data_2_2,Arreglo6])

#FILTRAMOS LOS GENES DESEADOS

fold_change_data_2_1.loc[list(Genes.loc[:,0])].to_csv('fold_change_genes_1_GSE149312_exp1.csv')
fold_change_data_2_1.loc[list(Genes2.loc[:,0])].to_csv('fold_change_genes_2_GSE149312_exp1.csv')
fold_change_data_2_1.loc[list(Genes3.loc[:,0])].to_csv('fold_change_genes_3_GSE149312_exp1.csv')
fold_change_data_2_2.loc[list(Genes.loc[:,0])].to_csv('fold_change_genes_1_GSE149312_exp2.csv')
fold_change_data_2_2.loc[list(Genes2.loc[:,0])].to_csv('fold_change_genes_2_GSE149312_exp2.csv')
fold_change_data_2_2.loc[list(Genes3.loc[:,0])].to_csv('fold_change_genes_3_GSE149312_exp2.csv')

Data2_corona_intestine_exp1_ndata.to_csv('Resultados_Analisis_GSE149312_exp1.csv')
Data2_corona_intestina_exp2_ndata.to_csv('Resultados_Analisis_GSE149312_exp2.csv')