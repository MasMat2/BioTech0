import pandas as pd
from numpy import nan
from Funciones import fold_change, Promedio

#IMPORTAMOS LOS GENES A FILTRAR

Genes = pd.read_csv('C:/Users/cogge/Desktop/Investigación/Genes/genes.csv', header = None)
Genes2 = pd.read_csv('C:/Users/cogge/Desktop/Investigación/Genes/genes-metabolismo.csv', header = None)
Genes3 = pd.read_excel('C:/Users/cogge/Desktop/Investigación/Genes/genes2Nuevo.xlsx',header = None)

#BASES DE DATOS GSE179277

Data0 = pd.read_excel('C:/Users/cogge/Desktop/Investigación/GSE179277/GSE179277_all_adult_and_ped_samples_combined_counts_unfiltered.xlsx')
Data0.rename(columns = {'Unnamed: 0':'gene_name_codigo','gene_name':'gene_name_descriptivo'},inplace = True)

#LA VARIABLE gene_name_codigo Y LA VARIABLE gene_name_descriptivo LES ASIGNAREMOS LAS COLUMNAS CON LOS CORRESPONDIENTES NOMBRES

gene_name_codigo = list(Data0['gene_name_codigo'])
gene_name_descriptivo = list(Data0['gene_name_descriptivo'])

#MATRIZ DE LA BASE DE DATOS

Matrix = pd.read_excel('C:/Users/cogge/Desktop/Investigación/GSE179277/GSE179277_series_matrix.xlsx',header = None)
Matrix.drop(index = list(range(4,82)),inplace = True)
Matrix.drop(columns = [0],inplace = True)
Matrix.index = ['Sexo','Edad','Infeccion','Muestra']
Matrix = Matrix.transpose()
Matrix['Sexo'] = list(string[-1] for string in Matrix['Sexo'])
Matrix['Edad'] = list(string2[-2:] for string2 in Matrix['Edad'])
Matrix['Edad'] = Matrix['Edad'].astype(int)
Infeccion = list()
for _ in list(Matrix['Infeccion']):
    if _.endswith('SC2') == True:
        Infeccion.append('SC2')
    elif _.endswith('no_virus') == True:
        Infeccion.append('no_virus')
    else:
        Infeccion.append('other_virus')
Matrix['Infeccion'] = Infeccion

#CREAMOS UNA COPIA DE LA BASE DE DATOS PARA TRABAJAR CON ELLA

Data = Data0.copy()
Data.index = gene_name_descriptivo
Data.drop(columns = ['gene_name_codigo','gene_name_descriptivo'],inplace = True)

#NORMALIZAMOS

for GAPDH,Columna in zip(list(Data.loc['GAPDH']),list(Data.columns)):
    Data[Columna] = Data[Columna]/GAPDH

#GUARDAREMOS LA FILA PRKY YA QUE CON ELLA SABREMOS SI LA COLUMNA CORRESPONDIENTE ES PROVENIENTE DE UN HOMBRE O UNA MUJER

PRKY = list(Data.loc['PRKY'])

#NORMALIZAMOS RESPECTO AL GEN 'GAPDH'

for GAPDH, Columna in zip (list(Data.loc['GAPDH']),list(Data.columns)):
    Data[Columna] = Data[Columna]/GAPDH

# SC2

#VAMOS A SEPARAR LOS GRUPOS A ESTUDIAR

Edad_0_12_SC2 = Data.loc[:,list(Matrix[(Matrix['Edad']>=0) & (Matrix['Edad']<=12) & (Matrix['Infeccion'] == 'SC2')]['Muestra'])].copy()
Edad_13_19_SC2 = Data.loc[:,list(Matrix[(Matrix['Edad']>=13) & (Matrix['Edad']<=19) & (Matrix['Infeccion'] == 'SC2')]['Muestra'])].copy()
Edad_20_mas_SC2 = Data.loc[:,list(Matrix[(Matrix['Edad']>=20) & (Matrix['Infeccion'] == 'SC2')]['Muestra'])].copy()

Mujeres_SC2 = Data.loc[:,list(Matrix[(Matrix['Sexo'] == 'F') & (Matrix['Infeccion'] == 'SC2')]['Muestra'])].copy()
Hombres_SC2 = Data.loc[:,list(Matrix[(Matrix['Sexo'] == 'M') & (Matrix['Infeccion'] == 'SC2')]['Muestra'])].copy()

Mujeres_edad_0_12_SC2 = Data.loc[:,list(Matrix[(Matrix['Sexo'] == 'F') & ((Matrix['Edad']>=0) & (Matrix['Edad']<=12)) & (Matrix['Infeccion'] == 'SC2')]['Muestra'])].copy()
Mujeres_edad_13_19_SC2 = Data.loc[:,list(Matrix[(Matrix['Sexo'] == 'F') & ((Matrix['Edad']>=13) & (Matrix['Edad']<=19)) & (Matrix['Infeccion'] == 'SC2')]['Muestra'])].copy()
Mujeres_edad_20_mas_SC2 = Data.loc[:,list(Matrix[(Matrix['Sexo'] == 'F') & (Matrix['Edad']>=20) & (Matrix['Infeccion'] == 'SC2')]['Muestra'])].copy()
Hombres_edad_0_12_SC2 = Data.loc[:,list(Matrix[(Matrix['Sexo'] == 'M') & ((Matrix['Edad']>=0) & (Matrix['Edad']<=12)) & (Matrix['Infeccion'] == 'SC2')]['Muestra'])].copy()
Hombres_edad_13_19_SC2 = Data.loc[:,list(Matrix[(Matrix['Sexo'] == 'M') & ((Matrix['Edad']>=13) & (Matrix['Edad']<=19)) & (Matrix['Infeccion'] == 'SC2')]['Muestra'])].copy()
Hombres_edad_20_mas_SC2 = Data.loc[:,list(Matrix[(Matrix['Sexo'] == 'M') & (Matrix['Edad']>=20) & (Matrix['Infeccion'] == 'SC2')]['Muestra'])].copy()

# Other_virus

#VAMOS A SEPARAR LOS GRUPOS A ESTUDIAR

Edad_0_12_other_virus = Data.loc[:,list(Matrix[(Matrix['Edad']>=0) & (Matrix['Edad']<=12) & (Matrix['Infeccion'] == 'other_virus')]['Muestra'])].copy()
Edad_13_19_other_virus = Data.loc[:,list(Matrix[(Matrix['Edad']>=13) & (Matrix['Edad']<=19) & (Matrix['Infeccion'] == 'other_virus')]['Muestra'])].copy()
Edad_20_mas_other_virus = Data.loc[:,list(Matrix[(Matrix['Edad']>=20) & (Matrix['Infeccion'] == 'other_virus')]['Muestra'])].copy()

Mujeres_other_virus = Data.loc[:,list(Matrix[(Matrix['Sexo'] == 'F') & (Matrix['Infeccion'] == 'other_virus')]['Muestra'])].copy()
Hombres_other_virus = Data.loc[:,list(Matrix[(Matrix['Sexo'] == 'M') & (Matrix['Infeccion'] == 'other_virus')]['Muestra'])].copy()

Mujeres_edad_0_12_other_virus = Data.loc[:,list(Matrix[(Matrix['Sexo'] == 'F') & ((Matrix['Edad']>=0) & (Matrix['Edad']<=12)) & (Matrix['Infeccion'] == 'other_virus')]['Muestra'])].copy()
Mujeres_edad_13_19_other_virus = Data.loc[:,list(Matrix[(Matrix['Sexo'] == 'F') & ((Matrix['Edad']>=13) & (Matrix['Edad']<=19)) & (Matrix['Infeccion'] == 'other_virus')]['Muestra'])].copy()
Mujeres_edad_20_mas_other_virus = Data.loc[:,list(Matrix[(Matrix['Sexo'] == 'F') & (Matrix['Edad']>=20) & (Matrix['Infeccion'] == 'other_virus')]['Muestra'])].copy()
Hombres_edad_0_12_other_virus = Data.loc[:,list(Matrix[(Matrix['Sexo'] == 'M') & ((Matrix['Edad']>=0) & (Matrix['Edad']<=12)) & (Matrix['Infeccion'] == 'other_virus')]['Muestra'])].copy()
Hombres_edad_13_19_other_virus = Data.loc[:,list(Matrix[(Matrix['Sexo'] == 'M') & ((Matrix['Edad']>=13) & (Matrix['Edad']<=19)) & (Matrix['Infeccion'] == 'other_virus')]['Muestra'])].copy()
Hombres_edad_20_mas_other_virus = Data.loc[:,list(Matrix[(Matrix['Sexo'] == 'M') & (Matrix['Edad']>=20) & (Matrix['Infeccion'] == 'other_virus')]['Muestra'])].copy()


#VAMOS A SEPARAR LOS CONTROLES DE CADA GRUPO

Control_Edad_0_12 = Data.loc[:,list(Matrix[(Matrix['Infeccion'] == 'no_virus') & ((Matrix['Edad']>=0) & (Matrix['Edad']<=12))]['Muestra'])].copy()
Control_Edad_13_19 = Data.loc[:,list(Matrix[(Matrix['Infeccion'] == 'no_virus') & ((Matrix['Edad']>=13) & (Matrix['Edad']<=19))]['Muestra'])].copy()
Control_Edad_20_mas = Data.loc[:,list(Matrix[(Matrix['Infeccion'] == 'no_virus') & (Matrix['Edad']>=20)]['Muestra'])].copy()

Control_Mujeres = Data.loc[:,list(Matrix[(Matrix['Infeccion'] == 'no_virus') & (Matrix['Sexo'] == 'F')]['Muestra'])].copy()
Control_Hombres = Data.loc[:,list(Matrix[(Matrix['Infeccion'] == 'no_virus') & (Matrix['Sexo'] == 'M')]['Muestra'])].copy()

Control_Mujeres_edad_0_12 = Data.loc[:,list(Matrix[(Matrix['Infeccion'] == 'no_virus') & ((Matrix['Sexo'] == 'F') & ((Matrix['Edad']>=0) & (Matrix['Edad']<=12)))]['Muestra'])].copy()
Control_Mujeres_edad_13_19 = Data.loc[:,list(Matrix[(Matrix['Infeccion'] == 'no_virus') & ((Matrix['Sexo'] == 'F') & ((Matrix['Edad']>=13) & (Matrix['Edad']<=19)))]['Muestra'])].copy()
Control_Mujeres_edad_20_mas = Data.loc[:,list(Matrix[(Matrix['Infeccion'] == 'no_virus') & ((Matrix['Sexo'] == 'F') & (Matrix['Edad']>=20))]['Muestra'])].copy()
Control_Hombres_edad_0_12 = Data.loc[:,list(Matrix[(Matrix['Infeccion'] == 'no_virus') & ((Matrix['Sexo'] == 'M') & ((Matrix['Edad']>=0) & (Matrix['Edad']<=12)))]['Muestra'])].copy()
Control_Hombres_edad_13_19 = Data.loc[:,list(Matrix[(Matrix['Infeccion'] == 'no_virus') & ((Matrix['Sexo'] == 'M') & ((Matrix['Edad']>=13) & (Matrix['Edad']<=19)))]['Muestra'])].copy()
Control_Hombres_edad_20_mas = Data.loc[:,list(Matrix[(Matrix['Infeccion'] == 'no_virus')& ((Matrix['Sexo'] == 'M') & (Matrix['Edad']>=20))]['Muestra'])].copy()

#OBTENEMOS LOS PROMEDIOS

#CREAMOS UNA FUNCION QUE LO HAGA PUES EN ESTA BASE DE DATOS YA CLASIFICAMOS POR SEPARADO LOS CONTROLES Y LOS GRUPOS A ESTUDIAR

# SC2

Promedio(Edad_0_12_SC2,'edad [0,12] - SC2') 
Promedio(Edad_13_19_SC2,'edad [13,19] - SC2') 
Promedio(Edad_20_mas_SC2,'edad [20,99] - SC2') 
Promedio(Mujeres_SC2,'mujeres - SC2') 
Promedio(Hombres_SC2,'hombres - SC2') 
Promedio(Mujeres_edad_0_12_SC2,'mujeres en edad [0,12] - SC2') 
Promedio(Mujeres_edad_13_19_SC2,'mujeres en edad [13,19] - SC2')
Promedio(Mujeres_edad_20_mas_SC2,'mujeres en edad [20,99] - SC2')
Promedio(Hombres_edad_0_12_SC2,'hombres en edad [0,12] - SC2') 
Promedio(Hombres_edad_13_19_SC2,'hombres en edad [13,19] - SC2')
Promedio(Hombres_edad_20_mas_SC2,'hombres en edad [20,99] - SC2')

# other_virus

Promedio(Edad_0_12_other_virus,'edad [0,12] - other virus') 
Promedio(Edad_13_19_other_virus,'edad [13,19] - other virus') 
Promedio(Edad_20_mas_other_virus,'edad [20,99] - other virus') 
Promedio(Mujeres_other_virus,'mujeres - other virus') 
Promedio(Hombres_other_virus,'hombres - other virus') 
Promedio(Mujeres_edad_0_12_other_virus,'mujeres en edad [0,12] - other virus') 
#Promedio(Mujeres_edad_13_19_other_virus,'mujeres en edad [13,19] - other virus') NO HAY MUJERES CON OTRO VIRUS EN ESE RANGO DE EDAD
Promedio(Mujeres_edad_20_mas_other_virus,'mujeres en edad [20,99] - other virus')
Promedio(Hombres_edad_0_12_other_virus,'hombres en edad [0,12] - other virus') 
Promedio(Hombres_edad_13_19_other_virus,'hombres en edad [13,19] - other virus')
Promedio(Hombres_edad_20_mas_other_virus,'hombres en edad [20,99] - other virus')

# Controles

Promedio(Control_Edad_0_12,'control edad [0,12]') 
Promedio(Control_Edad_13_19,'control edad [13,19]') 
Promedio(Control_Edad_20_mas,'control edad [20,99]')
Promedio(Control_Mujeres,'control mujeres') 
Promedio(Control_Hombres,'control hombres') 
Promedio(Control_Mujeres_edad_0_12,'control mujeres en edad [0,12]')
Promedio(Control_Mujeres_edad_13_19,'control mujeres en edad [13,19]')
Promedio(Control_Mujeres_edad_20_mas,'control mujeres en edad [20,99]')
Promedio(Control_Hombres_edad_0_12,'control hombres en edad [0,12]') 
Promedio(Control_Hombres_edad_13_19,'control hombres en edad [13,19]') 
Promedio(Control_Hombres_edad_20_mas,'control hombres en edad [20,99]')

#fold-change

# SC2

Edad_0_12_SC2['fold-change edad [0,12] - SC2'] = fold_change(Control_Edad_0_12['Promedio control edad [0,12]'],Edad_0_12_SC2['Promedio edad [0,12] - SC2'])
Edad_13_19_SC2['fold-change edad [13,19] - SC2'] = fold_change(Control_Edad_13_19['Promedio control edad [13,19]'],Edad_13_19_SC2['Promedio edad [13,19] - SC2'])
Edad_20_mas_SC2['fold-change edad [20,99] - SC2'] = fold_change(Control_Edad_20_mas['Promedio control edad [20,99]'],Edad_20_mas_SC2['Promedio edad [20,99] - SC2'])
Mujeres_SC2['fold-change mujeres - SC2'] = fold_change(Control_Mujeres['Promedio control mujeres'],Mujeres_SC2['Promedio mujeres - SC2'])
Hombres_SC2['fold-change hombres - SC2'] = fold_change(Control_Hombres['Promedio control hombres'],Hombres_SC2['Promedio hombres - SC2'])
Mujeres_edad_0_12_SC2['fold-change mujeres en edad [0,12] - SC2'] = fold_change(Control_Mujeres_edad_0_12['Promedio control mujeres en edad [0,12]'],Mujeres_edad_0_12_SC2['Promedio mujeres en edad [0,12] - SC2'])
Mujeres_edad_13_19_SC2['fold-change mujeres en edad [13,19] - SC2'] = fold_change(Control_Mujeres_edad_13_19['Promedio control mujeres en edad [13,19]'],Mujeres_edad_13_19_SC2['Promedio mujeres en edad [13,19] - SC2'])
Mujeres_edad_20_mas_SC2['fold-change mujeres en edad [20,99] - SC2'] = fold_change(Control_Mujeres_edad_20_mas['Promedio control mujeres en edad [20,99]'],Mujeres_edad_20_mas_SC2['Promedio mujeres en edad [20,99] - SC2'])
Hombres_edad_0_12_SC2['fold-change hombres en edad [0,12] - SC2'] = fold_change(Control_Hombres_edad_0_12['Promedio control hombres en edad [0,12]'],Hombres_edad_0_12_SC2['Promedio hombres en edad [0,12] - SC2'])
Hombres_edad_13_19_SC2['fold-change hombres en edad [13,19] - SC2'] = fold_change(Control_Hombres_edad_13_19['Promedio control hombres en edad [13,19]'],Hombres_edad_13_19_SC2['Promedio hombres en edad [13,19] - SC2'])
Hombres_edad_20_mas_SC2['fold-change hombres en edad [20,99] - SC2'] = fold_change(Control_Hombres_edad_20_mas['Promedio control hombres en edad [20,99]'],Hombres_edad_20_mas_SC2['Promedio hombres en edad [20,99] - SC2'])

# Other-virus

Edad_0_12_other_virus['fold-change edad [0,12] - other virus'] = fold_change(Control_Edad_0_12['Promedio control edad [0,12]'],Edad_0_12_other_virus['Promedio edad [0,12] - other virus'])
Edad_13_19_other_virus['fold-change edad [13,19] - other virus'] = fold_change(Control_Edad_13_19['Promedio control edad [13,19]'],Edad_13_19_other_virus['Promedio edad [13,19] - other virus'])
Edad_20_mas_other_virus['fold-change edad [20,99] - other virus'] = fold_change(Control_Edad_20_mas['Promedio control edad [20,99]'],Edad_20_mas_other_virus['Promedio edad [20,99] - other virus'])
Mujeres_other_virus['fold-change mujeres - other virus'] = fold_change(Control_Mujeres['Promedio control mujeres'],Mujeres_other_virus['Promedio mujeres - other virus'])
Hombres_other_virus['fold-change hombres - other virus'] = fold_change(Control_Hombres['Promedio control hombres'],Hombres_other_virus['Promedio hombres - other virus'])
Mujeres_edad_0_12_other_virus['fold-change mujeres en edad [0,12] - other virus'] = fold_change(Control_Mujeres_edad_0_12['Promedio control mujeres en edad [0,12]'],Mujeres_edad_0_12_other_virus['Promedio mujeres en edad [0,12] - other virus'])
#Mujeres_edad_13_19_other_virus['fold-change mujeres en edad [13,19] - other virus'] = fold_change(Control_Mujeres_edad_13_19['Promedio control mujeres en edad [13,19]'],Mujeres_edad_13_19_other_virus['Promedio mujeres en edad [13,19] - other virus'])
Mujeres_edad_20_mas_other_virus['fold-change mujeres en edad [20,99] - other virus'] = fold_change(Control_Mujeres_edad_20_mas['Promedio control mujeres en edad [20,99]'],Mujeres_edad_20_mas_other_virus['Promedio mujeres en edad [20,99] - other virus'])
Hombres_edad_0_12_other_virus['fold-change hombres en edad [0,12] - other virus'] = fold_change(Control_Hombres_edad_0_12['Promedio control hombres en edad [0,12]'],Hombres_edad_0_12_other_virus['Promedio hombres en edad [0,12] - other virus'])
Hombres_edad_13_19_other_virus['fold-change hombres en edad [13,19] - other virus'] = fold_change(Control_Hombres_edad_13_19['Promedio control hombres en edad [13,19]'],Hombres_edad_13_19_other_virus['Promedio hombres en edad [13,19] - other virus'])
Hombres_edad_20_mas_other_virus['fold-change hombres en edad [20,99] - other virus'] = fold_change(Control_Hombres_edad_20_mas['Promedio control hombres en edad [20,99]'],Hombres_edad_20_mas_other_virus['Promedio hombres en edad [20,99] - other virus'])

#UNIMOS LA INFORMACION OBTENIDA EN LA BASE DE DATOS

# SC2

Data = pd.concat([Data,Edad_0_12_SC2.loc[:,['Promedio edad [0,12] - SC2','fold-change edad [0,12] - SC2']]],axis = 1)
Data = pd.concat([Data,Edad_13_19_SC2.loc[:,['Promedio edad [13,19] - SC2','fold-change edad [13,19] - SC2']]],axis = 1)
Data = pd.concat([Data,Edad_20_mas_SC2.loc[:,['Promedio edad [20,99] - SC2','fold-change edad [20,99] - SC2']]],axis = 1)
Data = pd.concat([Data,Mujeres_SC2.loc[:,['Promedio mujeres - SC2','fold-change mujeres - SC2']]],axis = 1)
Data = pd.concat([Data,Hombres_SC2.loc[:,['Promedio hombres - SC2','fold-change hombres - SC2']]],axis = 1)
Data = pd.concat([Data,Mujeres_edad_0_12_SC2.loc[:,['Promedio mujeres en edad [0,12] - SC2','fold-change mujeres en edad [0,12] - SC2']]],axis = 1)
Data = pd.concat([Data,Mujeres_edad_13_19_SC2.loc[:,['Promedio mujeres en edad [13,19] - SC2','fold-change mujeres en edad [13,19] - SC2']]],axis = 1)
Data = pd.concat([Data,Mujeres_edad_20_mas_SC2.loc[:,['Promedio mujeres en edad [20,99] - SC2','fold-change mujeres en edad [20,99] - SC2']]],axis = 1)
Data = pd.concat([Data,Hombres_edad_0_12_SC2.loc[:,['Promedio hombres en edad [0,12] - SC2','fold-change hombres en edad [0,12] - SC2']]],axis = 1)
Data = pd.concat([Data,Hombres_edad_13_19_SC2.loc[:,['Promedio hombres en edad [13,19] - SC2','fold-change hombres en edad [13,19] - SC2']]],axis = 1)
Data = pd.concat([Data,Hombres_edad_20_mas_SC2.loc[:,['Promedio hombres en edad [20,99] - SC2','fold-change hombres en edad [20,99] - SC2']]],axis = 1)

# Other virus

Data = pd.concat([Data,Edad_0_12_other_virus.loc[:,['Promedio edad [0,12] - other virus','fold-change edad [0,12] - other virus']]],axis = 1)
Data = pd.concat([Data,Edad_13_19_other_virus.loc[:,['Promedio edad [13,19] - other virus','fold-change edad [13,19] - other virus']]],axis = 1)
Data = pd.concat([Data,Edad_20_mas_other_virus.loc[:,['Promedio edad [20,99] - other virus','fold-change edad [20,99] - other virus']]],axis = 1)
Data = pd.concat([Data,Mujeres_other_virus.loc[:,['Promedio mujeres - other virus','fold-change mujeres - other virus']]],axis = 1)
Data = pd.concat([Data,Hombres_other_virus.loc[:,['Promedio hombres - other virus','fold-change hombres - other virus']]],axis = 1)
Data = pd.concat([Data,Mujeres_edad_0_12_other_virus.loc[:,['Promedio mujeres en edad [0,12] - other virus','fold-change mujeres en edad [0,12] - other virus']]],axis = 1)
#Data = pd.concat([Data,Mujeres_edad_13_19_other_virus.loc[:,['Promedio mujeres en edad [13,19] - other virus','fold-change mujeres en edad [13,19] - other virus']]],axis = 1)
Data = pd.concat([Data,Mujeres_edad_20_mas_other_virus.loc[:,['Promedio mujeres en edad [20,99] - other virus','fold-change mujeres en edad [20,99] - other virus']]],axis = 1)
Data = pd.concat([Data,Hombres_edad_0_12_other_virus.loc[:,['Promedio hombres en edad [0,12] - other virus','fold-change hombres en edad [0,12] - other virus']]],axis = 1)
Data = pd.concat([Data,Hombres_edad_13_19_other_virus.loc[:,['Promedio hombres en edad [13,19] - other virus','fold-change hombres en edad [13,19] - other virus']]],axis = 1)
Data = pd.concat([Data,Hombres_edad_20_mas_other_virus.loc[:,['Promedio hombres en edad [20,99] - other virus','fold-change hombres en edad [20,99] - other virus']]],axis = 1)

# Controles

Data = pd.concat([Data,Control_Edad_0_12['Promedio control edad [0,12]'].to_frame()],axis = 1)
Data = pd.concat([Data,Control_Edad_13_19['Promedio control edad [13,19]'].to_frame()],axis = 1)
Data = pd.concat([Data,Control_Edad_20_mas['Promedio control edad [20,99]'].to_frame()],axis = 1)
Data = pd.concat([Data,Control_Mujeres['Promedio control mujeres'].to_frame()],axis = 1)
Data = pd.concat([Data,Control_Hombres['Promedio control hombres'].to_frame()],axis = 1)
Data = pd.concat([Data,Control_Mujeres_edad_0_12['Promedio control mujeres en edad [0,12]'].to_frame()],axis = 1)
Data = pd.concat([Data,Control_Mujeres_edad_13_19['Promedio control mujeres en edad [13,19]'].to_frame()], axis = 1)
Data = pd.concat([Data,Control_Mujeres_edad_20_mas['Promedio control mujeres en edad [20,99]'].to_frame()],axis = 1)
Data = pd.concat([Data,Control_Hombres_edad_0_12['Promedio control hombres en edad [0,12]'].to_frame()],axis = 1)
Data = pd.concat([Data,Control_Hombres_edad_13_19['Promedio control hombres en edad [13,19]'].to_frame()], axis =1)
Data = pd.concat([Data,Control_Hombres_edad_20_mas['Promedio control hombres en edad [20,99]'].to_frame()],axis = 1)

#EXPORTAMOS LA INFORMACION

FC_SC2 = ['fold-change edad [0,12] - SC2','fold-change edad [13,19] - SC2',
      'fold-change edad [20,99] - SC2','fold-change mujeres - SC2',
      'fold-change hombres - SC2','fold-change mujeres en edad [0,12] - SC2',
      'fold-change mujeres en edad [13,19] - SC2','fold-change mujeres en edad [20,99] - SC2',
      'fold-change hombres en edad [0,12] - SC2','fold-change hombres en edad [13,19] - SC2',
      'fold-change hombres en edad [20,99] - SC2']

FC_other_virus = ['fold-change edad [0,12] - other virus','fold-change edad [13,19] - other virus',
      'fold-change edad [20,99] - other virus','fold-change mujeres - other virus',
      'fold-change hombres - other virus','fold-change mujeres en edad [0,12] - other virus',
      'fold-change mujeres en edad [20,99] - other virus','fold-change hombres en edad [0,12] - other virus',
      'fold-change hombres en edad [13,19] - other virus','fold-change hombres en edad [20,99] - other virus']

#Primer arreglo para Genes

Arreglo0 = pd.DataFrame({_:[nan] for _ in list(Data.columns)}, 
                        index = ['VEGFA', 'HLA-F-AS1', 'HLA-H', 'HLA-J', 'HLA-L', 'VDR', 'TLR8-AS1', 'SOCS2-AS1'])

Data = pd.concat([Data,Arreglo0])

#SC2

Data.loc[list(Genes.loc[:,0]),FC_SC2].to_csv('fold_change_genes_1_SC2_GSE179277.csv')

#Other-virus

Data.loc[list(Genes.loc[:,0]),FC_other_virus].to_csv('fold_change_genes_1_other_virus_GSE179277.csv')

#Segundo arreglo para Genes2

Arreglo1 = pd.DataFrame({_:[nan] for _ in list(Data.columns)}, 
                        index = ['MDH1', 'ADSS', 'ADSSL1', 'ALPPL2'])

Data = pd.concat([Data,Arreglo1])

#SC2

Data.loc[list(Genes2.loc[:,0]),FC_SC2].to_csv('fold_change_genes_2_SC2_GSE179277.csv')

#Other-virus

Data.loc[list(Genes2.loc[:,0]),FC_other_virus].to_csv('fold_change_genes_2_other_virus_GSE179277.csv')

#Tercer arreglo para Genes3

Arreglo2 = pd.DataFrame({_:[nan] for _ in list(Data.columns)},
                        index = ['MDH1', 'ADSS', 'ADSSL1', 'ALPPL2', 'CR', 'DHNTP', 'GPX1'])

Data = pd.concat([Data,Arreglo2])

#SC2

Data.loc[list(Genes3.loc[:,0]),FC_SC2].to_csv('fold_change_genes_3_SC2_GSE179277.csv')

#Other-virus

Data.loc[list(Genes3.loc[:,0]),FC_other_virus].to_csv('fold_change_genes_3_other_virus_GSE179277.csv')

Data.to_csv('Resultados_Analisis_GSE179277.csv')