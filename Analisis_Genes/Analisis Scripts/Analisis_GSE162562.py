import pandas as pd
from numpy import nan
from Funciones import fold_change, Promedio

#IMPORTAMOS LOS GENES A FILTRAR

Genes = pd.read_csv('C:/Users/cogge/Desktop/Investigación/Genes/genes.csv', header = None)
Genes2 = pd.read_csv('C:/Users/cogge/Desktop/Investigación/Genes/genes-metabolismo.csv', header = None)
Genes3 = pd.read_excel('C:/Users/cogge/Desktop/Investigación/Genes/genes2Nuevo.xlsx',header = None)

#BASE DE DATOS

Matrix = pd.read_excel('C:/Users/cogge/Desktop/Investigación/GSE162562/GSE162562_series_matrix.txt/GSE162562_series_matrix.xlsx',
                       header = None)
Matrix.index = list(Matrix.loc[:,0])
Matrix.drop(columns = [0],inplace = True)
Matrix = Matrix.transpose()
Matrix.columns = ['Muestra','Edad','Sexo','Infeccion']
Matrix['Edad'] = list(int(str0[-2:]) for str0 in list(Matrix['Edad']))
Matrix['Sexo'] = list('Hombre' if str1.endswith('Male') == True else 'Mujer' for str1 in list(Matrix['Sexo']))
l0 = []
for str2 in list(Matrix['Infeccion']):
    if str2.endswith('asymptom') == True:
        l0.append('Asintomatico')
    elif str2.endswith('highly exposed') == True:
        l0.append('Control')
    elif str2.endswith('Mild') == True:
        l0.append('Leve')
Matrix['Infeccion'] = l0
Muestra = list(Matrix['Muestra'])

Data0 = pd.DataFrame()

for sample in Muestra:
    D = pd.read_excel('C:/Users/cogge/Desktop/Investigación/GSE162562/GSE162562_RAW/'+sample+'.xlsx',header = None)
    D.index = list(D.loc[:,0])
    D.drop(columns = [0],inplace = True)
    Data0 = pd.concat([Data0,D],axis = 1)
    Data0.rename(columns = {1:sample},inplace = True)

#CREAMOS UNA COPIA DE LA BASE DE DATOS PARA PODER MANIPULARLA

Data = Data0.copy()

#NORMALIZAMOS

for GAPDH,Columna in zip(list(Data.loc['GAPDH']),list(Data.columns)):
    Data[Columna] = Data[Columna]/GAPDH

#SEPARAMOS POR GRUPOS

Edad_0_12 = Data.loc[:,list(Matrix[ (Matrix['Edad']>=0) & (Matrix['Edad']<=12) ].loc[:,'Muestra'])].copy()
Edad_13_19 = Data.loc[:,list(Matrix[ (Matrix['Edad']>=13) & (Matrix['Edad']<=19)].loc[:,'Muestra'])].copy()
Edad_20_mas = Data.loc[:,list(Matrix[ Matrix['Edad']>=20 ].loc[:,'Muestra'])].copy()

Mujer = Data.loc[:,list(Matrix[ Matrix['Sexo'] == 'Mujer' ].loc[:,'Muestra'])].copy()
Hombre = Data.loc[:,list(Matrix[ Matrix['Sexo'] == 'Hombre' ].loc[:,'Muestra'])].copy()

Mujer_Edad_0_12 = Data.loc[:,list(Matrix[ (Matrix['Edad']>=0) & (Matrix['Edad']<=12) & (Matrix['Sexo'] == 'Mujer')].loc[:,'Muestra'])].copy()
Mujer_Edad_13_19 = Data.loc[:,list(Matrix[ (Matrix['Edad']>=13) & (Matrix['Edad']<=19) & (Matrix['Sexo'] == 'Mujer')].loc[:,'Muestra'])].copy()
Mujer_Edad_20_mas = Data.loc[:,list(Matrix[ (Matrix['Edad']>=20) & (Matrix['Sexo'] == 'Mujer') ].loc[:,'Muestra'])].copy()
Mujer_Edad_20_mas.drop(columns = ['A_118_Asymptom','A_932_Asymptom'],inplace = True)
A_118_Asymptom = Data.loc[:,'A_118_Asymptom'].to_frame().copy() #MUJER EN EDAD DE 20 O MAS
A_932_Asymptom = Data.loc[:,'A_932_Asymptom'].to_frame().copy() #MUJER EN EDAD DE 20 O MAS

Hombre_Edad_0_12 = Data.loc[:,list(Matrix[ (Matrix['Edad']>=0) & (Matrix['Edad']<=12) & (Matrix['Sexo'] == 'Hombre')].loc[:,'Muestra'])].copy()
Hombre_Edad_13_19 = Data.loc[:,list(Matrix[ (Matrix['Edad']>=13) & (Matrix['Edad']<=19) & (Matrix['Sexo'] == 'Hombre')].loc[:,'Muestra'])].copy()
Hombre_Edad_20_mas = Data.loc[:,list(Matrix[ (Matrix['Edad']>=20) & (Matrix['Sexo'] == 'Hombre')].loc[:,'Muestra'])].copy()

#CONTROLES DE CADA GRUPO

Control_Edad_0_12 = Data.loc[:,list(Matrix[ (Matrix['Infeccion'] == 'Control') & ((Matrix['Edad']>=0) & (Matrix['Edad']<=12)) ].loc[:,'Muestra'])].copy()
Control_Edad_13_19 = Data.loc[:,list(Matrix[ (Matrix['Infeccion'] == 'Control') & ((Matrix['Edad']>=13) & (Matrix['Edad']<=19)) ].loc[:,'Muestra'])].copy()
Control_Edad_20_mas = Data.loc[:,list(Matrix[ (Matrix['Infeccion'] == 'Control') & (Matrix['Edad']>=20) ].loc[:,'Muestra'])].copy()

Control_Mujer = Data.loc[:,list(Matrix[ (Matrix['Infeccion'] == 'Control') & (Matrix['Sexo'] == 'Mujer') ].loc[:,'Muestra'])].copy()
Control_Hombre = Data.loc[:,list(Matrix[ (Matrix['Infeccion'] == 'Control') & (Matrix['Sexo'] == 'Hombre') ].loc[:,'Muestra'])].copy()

Control_Mujer_Edad_0_12 = Data.loc[:,list(Matrix[ (Matrix['Infeccion'] == 'Control') & (Matrix['Sexo'] == 'Mujer') & ((Matrix['Edad']>=0) & (Matrix['Edad']<=12)) ].loc[:,'Muestra'])].copy()
Control_Mujer_Edad_13_19 = Data.loc[:,list(Matrix[ (Matrix['Infeccion'] == 'Control') & (Matrix['Sexo'] == 'Mujer') & ((Matrix['Edad']>=13) & (Matrix['Edad']<=19))].loc[:,'Muestra'])].copy()
Control_Mujer_Edad_20_mas = Data.loc[:,list(Matrix[ (Matrix['Infeccion'] == 'Control') & (Matrix['Sexo'] == 'Mujer') & (Matrix['Edad']>=20) ].loc[:,'Muestra'])].copy()

Control_Hombre_Edad_0_12 = Data.loc[:,list(Matrix[ (Matrix['Infeccion'] == 'Control') & (Matrix['Sexo'] == 'Hombre') & ((Matrix['Edad']>=0) & (Matrix['Edad']<=12)) ].loc[:,'Muestra'])].copy()
Control_Hombre_Edad_13_19 = Data.loc[:,list(Matrix[ (Matrix['Infeccion'] == 'Control') & (Matrix['Sexo'] == 'Hombre') & ((Matrix['Edad']>=13) & (Matrix['Edad']<=19)) ].loc[:,'Muestra'])].copy()
Control_Hombre_Edad_20_mas = Data.loc[:,list(Matrix[ (Matrix['Infeccion'] == 'Control') & (Matrix['Sexo'] == 'Hombre') & (Matrix['Edad']>=20) ].loc[:,'Muestra'])].copy()

#PROMEDIOS

Promedio(Edad_0_12,'edad [0,12]')
Promedio(Edad_13_19,'edad [13,19]')
Promedio(Edad_20_mas,'edad [20,99]')
Promedio(Mujer,'mujer')
Promedio(Hombre,'hombre')
Promedio(Mujer_Edad_0_12,'mujeres en edad [0,12]')
Promedio(Mujer_Edad_13_19,'mujeres en edad [13,19]')
Promedio(Mujer_Edad_20_mas,'mujeres en edad [20,99]')
Promedio(Hombre_Edad_0_12,'hombres en edad [0,12]')
Promedio(Hombre_Edad_13_19,'hombres en edad [13,19]')
Promedio(Hombre_Edad_20_mas,'hombres en edad [20,99]')

Promedio(Control_Edad_0_12,'control edad [0,12]')
Promedio(Control_Edad_13_19,'control edad [13,19]')
Promedio(Control_Edad_20_mas,'control edad [20,99]')
Promedio(Control_Mujer,'control mujer')
Promedio(Control_Hombre,'control hombre')
Promedio(Control_Mujer_Edad_0_12,'control mujeres en edad [0,12]')
Promedio(Control_Mujer_Edad_13_19,'control mujeres en edad [13,19]')
Promedio(Control_Mujer_Edad_20_mas,'control mujeres en edad [20,99]')
Promedio(Control_Hombre_Edad_0_12,'control hombres en edad [0,12]')
Promedio(Control_Hombre_Edad_13_19,'control hombres en edad [13,19]')
Promedio(Control_Hombre_Edad_20_mas,'control hombres en edad [20,99]')

#fold-change

Edad_0_12['fold-change edad [0,12]'] = fold_change(Control_Edad_0_12['Promedio control edad [0,12]'],Edad_0_12['Promedio edad [0,12]'])
Edad_13_19['fold-change edad [13,19]'] = fold_change(Control_Edad_13_19['Promedio control edad [13,19]'],Edad_13_19['Promedio edad [13,19]'])
Edad_20_mas['fold-change edad [20,99]'] = fold_change(Control_Edad_20_mas['Promedio control edad [20,99]'],Edad_20_mas['Promedio edad [20,99]'])
Mujer['fold-change mujer'] = fold_change(Control_Mujer['Promedio control mujer'],Mujer['Promedio mujer'])
Hombre['fold-change hombre'] = fold_change(Control_Hombre['Promedio control hombre'],Hombre['Promedio hombre'])
Mujer_Edad_0_12['fold-change mujer edad [0,12]'] = fold_change(Control_Mujer_Edad_0_12['Promedio control mujeres en edad [0,12]'],Mujer_Edad_0_12['Promedio mujeres en edad [0,12]'])
Mujer_Edad_13_19['fold-change mujer edad [13,19]'] = fold_change(Control_Mujer_Edad_13_19['Promedio control mujeres en edad [13,19]'],Mujer_Edad_13_19['Promedio mujeres en edad [13,19]'])
Mujer_Edad_20_mas['fold-change mujer edad [20,99]'] = fold_change(Control_Mujer_Edad_20_mas['Promedio control mujeres en edad [20,99]'],Mujer_Edad_20_mas['Promedio mujeres en edad [20,99]'])
Hombre_Edad_0_12['fold-change hombre edad [0,12]'] = fold_change(Control_Hombre_Edad_0_12['Promedio control hombres en edad [0,12]'],Hombre_Edad_0_12['Promedio hombres en edad [0,12]'])
Hombre_Edad_13_19['fold-change hombre edad [13,19]'] = fold_change(Control_Hombre_Edad_13_19['Promedio control hombres en edad [13,19]'],Hombre_Edad_13_19['Promedio hombres en edad [13,19]'])
Hombre_Edad_20_mas['fold-change hombre edad [20,99]'] = fold_change(Control_Hombre_Edad_20_mas['Promedio control hombres en edad [20,99]'],Hombre_Edad_20_mas['Promedio hombres en edad [20,99]'])

A_118_Asymptom['fold-change A_118_Asymptom'] = fold_change(Control_Mujer_Edad_20_mas['Promedio control mujeres en edad [20,99]'],A_118_Asymptom['A_118_Asymptom'])
A_932_Asymptom['fold-change A_932_Asymptom'] = fold_change(Control_Mujer_Edad_20_mas['Promedio control mujeres en edad [20,99]'],A_932_Asymptom['A_932_Asymptom'])

#JUNTAMOS LA INFORMACION

Data = pd.concat([Data,Edad_0_12.loc[:,['Promedio edad [0,12]','fold-change edad [0,12]']]],axis = 1)
Data = pd.concat([Data,Edad_13_19.loc[:,['Promedio edad [13,19]','fold-change edad [13,19]']]],axis = 1)
Data = pd.concat([Data,Edad_20_mas.loc[:,['Promedio edad [20,99]','fold-change edad [20,99]']]],axis = 1)

Data = pd.concat([Data,Mujer.loc[:,['Promedio mujer','fold-change mujer']]],axis = 1)
Data = pd.concat([Data,Hombre.loc[:,['Promedio hombre','fold-change hombre']]],axis = 1)

Data = pd.concat([Data,Mujer_Edad_0_12.loc[:,['Promedio mujeres en edad [0,12]','fold-change mujer edad [0,12]']]],axis = 1)
Data = pd.concat([Data,Mujer_Edad_13_19.loc[:,['Promedio mujeres en edad [13,19]','fold-change mujer edad [13,19]']]],axis = 1)
Data = pd.concat([Data,Mujer_Edad_20_mas.loc[:,['Promedio mujeres en edad [20,99]','fold-change mujer edad [20,99]']]],axis = 1)

Data = pd.concat([Data,Hombre_Edad_0_12.loc[:,['Promedio hombres en edad [0,12]','fold-change hombre edad [0,12]']]], axis = 1)
Data = pd.concat([Data,Hombre_Edad_13_19.loc[:,['Promedio hombres en edad [13,19]','fold-change hombre edad [13,19]']]],axis = 1)
Data = pd.concat([Data,Hombre_Edad_20_mas.loc[:,['Promedio hombres en edad [20,99]','fold-change hombre edad [20,99]']]],axis = 1)

Data = pd.concat([Data,A_118_Asymptom.loc[:,'fold-change A_118_Asymptom'].to_frame()],axis = 1)
Data = pd.concat([Data,A_932_Asymptom.loc[:,'fold-change A_932_Asymptom'].to_frame()],axis = 1)

#EXPORTAMOS LA INFORMACION

f_c = ['fold-change edad [0,12]','fold-change edad [13,19]',
       'fold-change edad [20,99]','fold-change mujer',
       'fold-change hombre','fold-change mujer edad [0,12]',
       'fold-change mujer edad [13,19]','fold-change mujer edad [20,99]',
       'fold-change hombre edad [0,12]','fold-change hombre edad [13,19]','fold-change hombre edad [20,99]']

Data.loc[list(Genes.loc[:,0]),f_c].to_csv('fold_change_genes_1_GSE162562.csv')
Data.loc[list(Genes2.loc[:,0]),f_c].to_csv('fold_change_genes_2_GSE162562.csv')

#Arreglo para genes3

Arreglo = pd.DataFrame({_:[nan] for _ in list(Data.columns)},index = ['CR', 'DHNTP'])
Data = pd.concat([Data,Arreglo])

Data.loc[list(Genes3.loc[:,0]),f_c].to_csv('fold_change_genes_3_GSE162562.csv')
Data.to_csv('Resultados_Analisis_GSE162562.csv')