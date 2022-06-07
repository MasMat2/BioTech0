import pandas as pd
import statistics as sts
from numpy import nan
from Funciones import fold_change

#IMPORTAMOS LOS GENES QUE NOS INTERESAN

Genes1 = pd.read_csv('./Genes/genes-metabolismo.csv', header = None)
Genes2 = pd.read_csv('./Genes/genes-metabolismo.csv', header = None)
Genes3 = pd.read_excel('./Genes/genes2Nuevo.xlsx',header = None)
Genes = set(Genes1[0]).union(set(Genes2[0])).union(Genes3[0])

# Leemos la base de datos
# Eliminamos series 3, 4, 8
# Eliminanos Series9_NHBE_IAV_1 - Series9_NHBE_IAV_4
# Eliminanos Series9_NHBE_IAVdNS1_1 - Series9_NHBE_IAVdNS1_4
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
    Series[id] = Data.loc[:,list(S for S in list(Data.columns) if f"Series{id}_" in S[:9])].copy()

#PROMEDIOS

for key in Series.keys():
    measurement_types = set(map(lambda x: x[:-2], Series[key]))
    for typ in measurement_types:
        Series[key][f"Promedio_{typ}"] = list(sts.mean(Series[key].filter(like=typ, axis=1).iloc[fila]) for fila in range (Total_Gen))
    
#fold-change
for key in Series.keys():
    # Remove else False
    measurement_types = [x for x in Series[key] if "Promedio" in x]
    if key != "15":
        mock = Series[key].filter(regex="^Promedio.*Mock$").columns[0]
    else:
        mock = Series[key].filter(regex="Promedio_Series15_HealthyLungBiopsy").columns[0]
    measurement_types.remove(mock)

    if len(measurement_types) > 1:
        for typ in measurement_types:
            name = typ.split(f"_Series{key}_")[1].replace("_", " ")
            Series[key][f"fold-change Series{key} {name}"] = fold_change(Series[key][mock], Series[key][typ])
    else:
        typ = measurement_types[0]
        Series[key][f"fold-change Series{key}"] = fold_change(Series[key][mock], Series[key][typ])

#CONJUNTO TOTAL
for key in Series.keys():
    Data = pd.concat([Data, Series[key].filter(like="Promedio_Series", axis=1)], axis=1)
    Data = pd.concat([Data, Series[key].filter(like="fold-change", axis=1)], axis=1)

#DF fold-change
fold_change_data = Data.filter(like='fold-change')

#METEMOS LAS FILAS DE LOS GENES QUE NO EXISTEN EN GSE147507 PERO SI EN EL FILTRADO DE GENES DESEADOS
missing_genes = Genes-set(fold_change_data.index)
fold_change_data = pd.DataFrame(fold_change_data, index=list(fold_change_data.index)+list(missing_genes), columns=fold_change_data.columns)

#FILTRAMOS LOS GENES DESEADOS DEL DATA FRAME fold_change_data EN FUNCION DE LOS VALORES EN Genes, Genes2 Y Genes3

fold_change_data.loc[Genes1[0]].to_csv('fold_change_genes_1_GSE147507.csv')
fold_change_data.loc[Genes2[0]].to_csv('fold_change_genes_2_GSE147507.csv')
fold_change_data.loc[Genes3[0]].to_csv('fold_change_genes_3_GSE147507.csv')
Data.to_csv('Resultados_Analisis_GSE147507.csv')