import pandas as pd
import statistics as sts
from numpy import nan
from Funciones import fold_change

#IMPORTAMOS LOS GENES A FILTRAR

Genes = pd.read_csv('C:/Users/cogge/Desktop/Investigación/Genes/genes.csv', header = None)
Genes2 = pd.read_csv('C:/Users/cogge/Desktop/Investigación/Genes/genes-metabolismo.csv', header = None)
Genes3 = pd.read_excel('C:/Users/cogge/Desktop/Investigación/Genes/genes2Nuevo.xlsx',header = None)

#BASE DE DATOS GSE147507

Data0 = pd.read_excel('GSE147507/GSE147507.xlsx')
Data0.index = list(Data0['Unnamed: 0'])
Data0.drop(columns = ['Unnamed: 0']+[_ for _ in list(Data0.columns) if _[6] in ['3','4','8']],
          inplace = True)
Data0.drop(columns = list(Data0.loc[:,'Series9_NHBE_IAV_1':'Series9_NHBE_IAV_4'])+list(Data0.loc[:,'Series9_NHBE_IAVdNS1_1':'Series9_NHBE_IAVdNS1_4']), 
          inplace = True )

#CREAMOS UNA COPIA CON LA CUAL VAMOS A TRABAJAR

Data = Data0.copy()

#TOTAL DE GENES EN GSE147507

Total_Gen = len(Data.index)

#NORMALIZACION

for GAPDH,Columna in zip(list(Data.loc['GAPDH']),list(Data.columns)) :
    Data[Columna] = Data[Columna]/GAPDH
    
series_ids = ["1", "2", "5", "6", "7", "9", "15", "16"]

Series = {}

for id in series_ids:
    Series[id] = Data.loc[:,list(S for S in list(Data.columns) if id in S[6:9])].copy()

    
#DIVIDIMOS POR SERIES

Series1 = Data.loc[:,list(S1 for S1 in list(Data.columns) if S1[6] == '1' and S1[7] == '_')].copy()
Series2 = Data.loc[:,list(S2 for S2 in list(Data.columns) if S2[6] == '2' and S2[7] == '_')].copy()
Series5 = Data.loc[:,list(S5 for S5 in list(Data.columns) if S5[6] == '5' and S5[7] == '_')].copy()
Series6 = Data.loc[:,list(S6 for S6 in list(Data.columns) if S6[6] == '6' and S6[7] == '_')].copy()
Series7 = Data.loc[:,list(S7 for S7 in list(Data.columns) if S7[6] == '7' and S7[7] == '_')].copy()
Series9 = Data.loc[:,list(S9 for S9 in list(Data.columns) if S9[6] == '9' and S9[7] == '_')].copy()
Series15 = Data.loc[:,list(S15 for S15 in list(Data.columns) if S15[6] == '1' and S15[7] == '5')].copy()
Series16 = Data.loc[:,list(S16 for S16 in list(Data.columns) if S16[6] == '1' and S16[7] == '6')].copy()


#PROMEDIOS

Series1['Promedio_Series1_NHBE_Mock'] = list(sts.mean(Series1.iloc[:,0:3].iloc[Fila11]) for Fila11 in range (Total_Gen))
Series1['Promedio_Series1_NHBE_SARS_CoV-2'] = list(sts.mean(Series1.iloc[:,3:6].iloc[Fila12]) for Fila12 in range (Total_Gen))

Series2['Promedio_Series2_A549_Mock'] = list(sts.mean(Series2.iloc[:,0:3].iloc[Fila21]) for Fila21 in range (Total_Gen))
Series2['Promedio_Series2_A549_SARS_CoV-2'] = list(sts.mean(Series2.iloc[:,3:6].iloc[Fila22]) for Fila22 in range (Total_Gen))

Series5['Promedio_Series5_A549_Mock'] = list(sts.mean(Series5.iloc[:,0:3].iloc[Fila51]) for Fila51 in range (Total_Gen))
Series5['Promedio_Series5_A549_SARS_CoV-2'] = list(sts.mean(Series5.iloc[:,3:6].iloc[Fila52]) for Fila52 in range (Total_Gen))

Series6['Promedio_Series6_ACE2_Mock'] = list(sts.mean(Series6.iloc[:,0:3].iloc[Fila61]) for Fila61 in range (Total_Gen))
Series6['Promedio_Series6_ACE2_SARS_CoV-2'] = list(sts.mean(Series6.iloc[:,3:6].iloc[Fila62]) for Fila62 in range (Total_Gen))

Series7['Promedio_Series7_Calu3_Mock'] = list(sts.mean(Series7.iloc[:,0:3].iloc[Fila71]) for Fila71 in range (Total_Gen))
Series7['Promedio_Series7_Calu3_SARS-CoV-2'] = list(sts.mean(Series7.iloc[:,3:6].iloc[Fila72]) for Fila72 in range (Total_Gen))

Series9['Promedio_Series9_NHBE_Mock'] = list(sts.mean(Series9.loc[:,'Series9_NHBE_Mock_1':'Series9_NHBE_Mock_4'].iloc[Fila91]) for Fila91 in range (Total_Gen))
Series9['Promedio_Series9_NHBE_IFNB_4h'] = list(sts.mean(Series9.loc[:,['Series9_NHBE_IFNB_4h_1','Series9_NHBE_IFNB_4h_2']].iloc[Fila94]) for Fila94 in range (Total_Gen))
Series9['Promedio_Series9_NHBE_IFNB_6h'] = list(sts.mean(Series9.loc[:,['Series9_NHBE_IFNB_6h_1','Series9_NHBE_IFNB_6h_2']].iloc[Fila95]) for Fila95 in range (Total_Gen))
Series9['Promedio_Series9_NHBE_IFNB_12h'] = list(sts.mean(Series9.loc[:,['Series9_NHBE_IFNB_12h_1','Series9_NHBE_IFNB_12h_2']].iloc[Fila96])for Fila96 in range (Total_Gen))

Series15['Promedio_Series15_HealthyLungBiopsy'] = list(sts.mean(Series15.iloc[:,[0,1]].iloc[Fila151]) for Fila151 in range (Total_Gen))
Series15['Promedio_Series15_COVID19Lung'] = list(sts.mean(Series15.iloc[:,[2,3]].iloc[Fila152]) for Fila152 in range (Total_Gen))

Series16['Promedio_Series16_A549-ACE2_Mock'] = list(sts.mean(Series16.iloc[:,[0,1,2]].iloc[Fila161]) for Fila161 in range (Total_Gen))
Series16['Promedio_Series16_A549-ACE2_SARS-CoV-2'] = list(sts.mean(Series16.iloc[:,[3,4,5]].iloc[Fila162]) for Fila162 in range (Total_Gen))
Series16['Promedio_Series16_A549-ACE2_SARS-CoV-2_Rux'] = list(sts.mean(Series16.iloc[:,[6,7,8]].iloc[Fila163]) for Fila163 in range (Total_Gen))

#fold-change

Series1['fold-change Series1'] = fold_change(Series1['Promedio_Series1_NHBE_Mock'],Series1['Promedio_Series1_NHBE_SARS_CoV-2'])

Series2['fold-change Series2'] = fold_change(Series2['Promedio_Series2_A549_Mock'],Series2['Promedio_Series2_A549_SARS_CoV-2'])

Series5['fold-change Series5'] = fold_change(Series5['Promedio_Series5_A549_Mock'],Series5['Promedio_Series5_A549_SARS_CoV-2'])

Series6['fold-change Series6'] = fold_change(Series6['Promedio_Series6_ACE2_Mock'],Series6['Promedio_Series6_ACE2_SARS_CoV-2'])

Series7['fold-change Series7'] = fold_change(Series7['Promedio_Series7_Calu3_Mock'],Series7['Promedio_Series7_Calu3_SARS-CoV-2'])

Series9['fold-change Series9 NHBE_IFNB_4h'] = fold_change(Series9['Promedio_Series9_NHBE_Mock'],Series9['Promedio_Series9_NHBE_IFNB_4h'])
Series9['fold-change Series9 NHBE_IFNB_6h'] = fold_change(Series9['Promedio_Series9_NHBE_Mock'],Series9['Promedio_Series9_NHBE_IFNB_6h'])
Series9['fold-change Series9 NHBE_IFNB_12h'] = fold_change(Series9['Promedio_Series9_NHBE_Mock'],Series9['Promedio_Series9_NHBE_IFNB_12h'])

Series15['fold-change Series15'] = fold_change(Series15['Promedio_Series15_HealthyLungBiopsy'],Series15['Promedio_Series15_COVID19Lung'])

Series16['fold-change Series16 A549-ACE2_SARS-CoV-2'] = fold_change(Series16['Promedio_Series16_A549-ACE2_Mock'],Series16['Promedio_Series16_A549-ACE2_SARS-CoV-2'])
Series16['fold-change Series16 A549-ACE2_SARS-CoV-2_Rux'] = fold_change(Series16['Promedio_Series16_A549-ACE2_Mock'],Series16['Promedio_Series16_A549-ACE2_SARS-CoV-2_Rux'])

#CONJUNTO TOTAL

Data = pd.concat([Data,Series1.loc[:,['Promedio_Series1_NHBE_Mock','Promedio_Series1_NHBE_SARS_CoV-2','fold-change Series1']]], axis = 1)
Data = pd.concat([Data,Series2.loc[:,['Promedio_Series2_A549_Mock','Promedio_Series2_A549_SARS_CoV-2','fold-change Series2']]], axis = 1)
Data = pd.concat([Data,Series5.loc[:,['Promedio_Series5_A549_Mock','Promedio_Series5_A549_SARS_CoV-2','fold-change Series5']]], axis = 1)
Data = pd.concat([Data,Series6.loc[:,['Promedio_Series6_ACE2_Mock','Promedio_Series6_ACE2_SARS_CoV-2','fold-change Series6']]], axis = 1)
Data = pd.concat([Data,Series7.loc[:,['Promedio_Series7_Calu3_Mock','Promedio_Series7_Calu3_SARS-CoV-2','fold-change Series7']]], axis = 1)
Data = pd.concat([Data,Series9.loc[:,'Promedio_Series9_NHBE_Mock':'fold-change Series9 NHBE_IFNB_12h']], axis = 1)
Data = pd.concat([Data,Series15.loc[:,'Promedio_Series15_HealthyLungBiopsy':'fold-change Series15']], axis = 1)
Data = pd.concat([Data,Series16.loc[:,'Promedio_Series16_A549-ACE2_Mock':'fold-change Series16 A549-ACE2_SARS-CoV-2_Rux']], axis = 1)

#DF fold-change

fold_change_data = Data.loc[:,list(_ for _ in list(Data.columns) if _.startswith('fold-change'))]

#METEMOS LAS FILAS DE LOS GENES QUE NO EXISTEN EN GSE147507 PERO SI EN EL FILTRADO DE GENES DESEADOS

ALPPL2 = pd.DataFrame({_:[nan] for _ in list(fold_change_data.columns)}, index = ['ALPPL2'])
CR = pd.DataFrame({_:[nan] for _ in list(fold_change_data.columns)}, index = ['CR'])
DHNTP = pd.DataFrame({_:[nan] for _ in list(fold_change_data.columns)}, index = ['DHNTP'])
fold_change_data = pd.concat([fold_change_data,ALPPL2])
fold_change_data = pd.concat([fold_change_data,CR])
fold_change_data = pd.concat([fold_change_data,DHNTP])

#FILTRAMOS LOS GENES DESEADOS DEL DATA FRAME fold_change_data EN FUNCION DE LOS VALORES EN Genes, Genes2 Y Genes3

fold_change_data.loc[list(Genes.loc[:,0])].to_csv('fold_change_genes_1_GSE147507.csv')
fold_change_data.loc[list(Genes2.loc[:,0])].to_csv('fold_change_genes_2_GSE147507.csv')
fold_change_data.loc[list(Genes3.loc[:,0])].to_csv('fold_change_genes_3_GSE147507.csv')
Data.to_csv('Resultados_Analisis_GSE147507.csv')