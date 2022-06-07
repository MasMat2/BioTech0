import pandas as pd
import statistics as sts
from numpy import nan
from Funciones import fold_change

#IMPORTAMOS LOS GENES A FILTRAR

Genes = pd.read_csv('C:/Users/cogge/Desktop/Investigación/Genes/genes.csv', header = None)
Genes2 = pd.read_csv('C:/Users/cogge/Desktop/Investigación/Genes/genes-metabolismo.csv', header = None)
Genes3 = pd.read_excel('C:/Users/cogge/Desktop/Investigación/Genes/genes2Nuevo.xlsx',header = None)

#BASES DE DATOS GSE150316

Data3_0 = pd.read_excel('C:/Users/cogge/Desktop/Investigación/GSE150316/GSE150316_RPMNormCounts_final.xlsx')
Data3_0.index = list(Data3_0['Unnamed: 0'])
Data3_0.drop(columns = ['Unnamed: 0'], inplace = True)

#CREAMOS UNA COPIA PARA TRABAJAR EL ANALISIS

Data3 = Data3_0.copy()

#NORMALIZAMOS LA BD CON Y SIN GENERO

for GAPDH,Columna in zip(list(Data3.loc['GAPDH']),list(Data3.columns)) :
    Data3[Columna] = Data3[Columna]/GAPDH

#FILTRAMOS LOS CONTROLES y GENES EN LA BD GSE150316

controles = Data3.loc[:,'NegControl1':'NegControl5'].copy()
GenData3 = list(Data3_0.index)

##############################################################################
#ANALISIS SIN GENERO

#SEPARAMOS POR CASOS

case1 = Data3.loc[:,'case1-lung1':'case1-heart1'].copy()
case2 = Data3.loc[:,'case2-lung1':'case2-heart1'].copy()
case3 = Data3.loc[:,'case3-lung1':'case3-lung2'].copy()
case3 = pd.concat([case3,Data3.loc[:,'case3-liver1']], axis = 1)
case4 = Data3.loc[:,'case4-lung1':'case4-bowel1'].copy()
case5 = Data3.loc[:,'case5-lung1':'case5-marrow1'].copy()

#PROMEDIOS

case1['Promedio case1 lung'] = list(sts.mean(case1.loc[case1lung,['case1-lung1','case1-lung2','case1-lung3','case1-lung4']]) for case1lung in GenData3)
case1['Promedio case1 heart'] = list(case1['case1-heart1'])

case2['Promedio case2 lung'] = list(sts.mean(case2.loc[case2lung,['case2-lung1','case2-lung2','case2-lung3']]) for case2lung in GenData3)
case2['Promedio case2 jejunum'] = list(case2['case2-jejunum1'])
case2['Promedio case2 heart'] = list(case2['case2-heart1'])

case3['Promedio case3 lung'] = list(sts.mean(case3.loc[case3lung,['case3-lung1','case3-lung2']]) for case3lung in GenData3)
case3['Promedio case3 heart'] = list(case3['case3-heart1'])
case3['Promedio case3 liver'] = list(case3['case3-liver1'])

case4['Promedio case4 lung'] = list(sts.mean(case4.loc[case4lung,['case4-lung1','case4-lung2']]) for case4lung in GenData3)
case4['Promedio case4 heart'] = list(sts.mean(case4.loc[case4heart,['case4-heart1','case4-heart2']]) for case4heart in GenData3)
case4['Promedio case4 liver'] = list(case4['case4-liver1'])
case4['Promedio case4 kidney'] = list(case4['case4-kidney1'])
case4['Promedio case4 bowel'] = list(case4['case4-bowel1'])

case5['Promedio case5 lung'] = list(sts.mean(case5.loc[case5lung,['case5-lung1','case5-lung2','case5-lung3','case5-lung4','case5-lung5']]) for case5lung in GenData3)
case5['Promedio case5 fat'] = list(case5['case5-fat1'])
case5['Promedio case5 skin'] = list(case5['case5-skin1'])
case5['Promedio case5 bowel'] = list(case5['case5-bowel1'])
case5['Promedio case5 liver'] = list(case5['case5-liver1'])
case5['Promedio case5 kidney'] = list(case5['case5-kidney1'])
case5['Promedio case5 heart'] = list(case5['case5-heart1'])
case5['Promedio case5 marrow'] = list(case5['case5-marrow1'])

#fold-change

case1['fold-change control1 case1 lung'] = fold_change(controles['NegControl1'],case1['Promedio case1 lung'])
case1['fold-change control1 case1 heart'] = fold_change(controles['NegControl1'],case1['Promedio case1 heart'])

case2['fold-change control2 case2 lung'] = fold_change(controles['NegControl2'],case2['Promedio case2 lung'])
case2['fold-change control2 case2 jejunum'] = fold_change(controles['NegControl2'],case2['Promedio case2 jejunum'])
case2['fold-change control2 case2 heart'] = fold_change(controles['NegControl2'],case2['Promedio case2 heart'])

case3['fold-change control3 case3 lung'] = fold_change(controles['NegControl3'],case3['Promedio case3 lung'])
case3['fold-change control3 case3 heart'] = fold_change(controles['NegControl3'],case3['Promedio case3 heart'])
case3['fold-change control3 case3 liver'] = fold_change(controles['NegControl3'],case3['Promedio case3 liver'])

case4['fold-change control4 case4 lung'] = fold_change(controles['NegControl4'],case4['Promedio case4 lung'])
case4['fold-change control4 case4 heart'] = fold_change(controles['NegControl4'],case4['Promedio case4 heart'])
case4['fold-change control4 case4 liver'] = fold_change(controles['NegControl4'],case4['Promedio case4 liver'])
case4['fold-change control4 case4 kidney'] = list(case4['case4-kidney1'])
case4['fold-change control4 case4 bowel'] = fold_change(controles['NegControl4'],case4['Promedio case4 bowel'])

case5['fold-change control5 case5 lung'] = fold_change(controles['NegControl5'],case5['Promedio case5 lung'])
case5['fold-change control5 case5 fat'] = fold_change(controles['NegControl5'],case5['Promedio case5 fat'])
case5['fold-change control5 case5 skin'] = fold_change(controles['NegControl5'],case5['Promedio case5 skin'])
case5['fold-change control5 case5 bowel'] = fold_change(controles['NegControl5'],case5['Promedio case5 bowel'])
case5['fold-change control5 case5 liver'] = fold_change(controles['NegControl5'],case5['Promedio case5 liver'])
case5['fold-change control5 case5 kidney'] = fold_change(controles['NegControl5'],case5['Promedio case5 kidney'])
case5['fold-change control5 case5 heart'] = fold_change(controles['NegControl5'],case5['Promedio case5 heart'])
case5['fold-change control5 case5 marrow'] = fold_change(controles['NegControl5'],case5['Promedio case5 marrow'])

#ANEXAMOS LOS fold-change A Data3

Data3 = pd.concat([Data3,case1.loc[:,['fold-change control1 case1 lung','fold-change control1 case1 heart']]], axis = 1)
Data3 = pd.concat([Data3,case2.loc[:,['fold-change control2 case2 lung','fold-change control2 case2 jejunum','fold-change control2 case2 heart']]], axis = 1)
Data3 = pd.concat([Data3,case3.loc[:,['fold-change control3 case3 lung','fold-change control3 case3 liver']]], axis = 1)
Data3 = pd.concat([Data3,case4.loc[:,'fold-change control4 case4 lung':'fold-change control4 case4 bowel']],axis = 1)
Data3 = pd.concat([Data3,case5.loc[:,'fold-change control5 case5 lung':'fold-change control5 case5 marrow']], axis = 1)

#############################################################################
#ANALISIS CON GENERO

#SEPARAMOS LOS CONTROLES POR GENERO

controles['Control Mujer']  = list(sts.mean(controles.loc[ControlM,['NegControl3','NegControl4']]) for ControlM in GenData3)
controles['Control Hombre'] = list(sts.mean(controles.loc[ControlH,['NegControl1','NegControl2','NegControl5']]) for ControlH in GenData3)

#SEPARAMOS POR CASOS PARA EL ESTUDIO CON GENERO

case6 = Data3.loc[:,'case6-lung1':'case6-lung5'].copy()
case7 = Data3.loc[:,'case7-lung1':'case7-lung5'].copy()
case8 = Data3.loc[:,'case8-lung1':'case8-heart1'].copy()
case9 = Data3.loc[:,'case9-lung1':'case9-lung5'].copy()
case10 = Data3.loc[:,'case10-lung1':'case10-lung3'].copy()
case10 = pd.concat([case10,Data3.loc[:,'case10-liver1']], axis = 1)
case11 = Data3.loc[:,'case11-lung1':'case11-bowel1'].copy()
case12 = Data3.loc[:,'case12-liver1'].to_frame().copy()
caseA = Data3.loc[:,'caseA-lung-NYC'].to_frame().copy()
caseB = Data3.loc[:,'caseB-lung-NYC'].to_frame().copy()
caseC = Data3.loc[:,'caseC-lung-NYC'].to_frame().copy()
caseD = Data3.loc[:,'caseD-lung-NYC'].to_frame().copy()
caseF = Data3.loc[:,'caseF-lung-NYC'].to_frame().copy()
caseE = Data3.loc[:,'caseE-lung-NYC'].to_frame().copy()
caseG = Data3.loc[:,'caseG-lung-NYC'].to_frame().copy()
caseH = Data3.loc[:,'caseH-lung-NYC'].to_frame().copy()
caseI = Data3.loc[:,'caseI-lung-NYC'].to_frame().copy()
caseJ = Data3.loc[:,'caseJ-lung-NYC'].to_frame().copy()
caseP1 = Data3.loc[:,'caseP1-placenta1':'caseP1-placenta2'].copy()
caseP2 = Data3.loc[:,'caseP2-placenta1':'caseP2-placenta2'].copy()
caseP3 = Data3.loc[:,'caseP3-placenta1':'caseP3-placenta2'].copy()
caseP4 = Data3.loc[:,'caseP4-placenta1'].to_frame().copy()

#PROMEDIOS DE LOS CASOS PARA ESTUDIAR CON GENERO

case6['Promedio case6 lung'] = list(sts.mean(case6.loc[case6lung]) for case6lung in GenData3)

case7['Promedio case7 lung'] = list(sts.mean(case7.loc[case7lung]) for case7lung in GenData3)

case8['Promedio case8 lung'] = list(sts.mean(case8.loc[case8lung,'case8-lung1':'case8-lung5']) for case8lung in GenData3)
case8['Promedio case8 liver'] = list(case8['case8-liver1'])
case8['Promedio case8 bowel'] = list(case8['case8-bowel1'])
case8['Promedio case8 heart'] = list(case8['case8-heart1'])

case9['Promedio case9 lung'] = list(sts.mean(case9.loc[case9lung]) for case9lung in GenData3)

case10['Promedio case10 lung'] = list(sts.mean(case10.loc[case10lung,'case10-lung1':'case10-lung3']) for case10lung in GenData3)
case10['Promedio case10 liver'] = list(case10['case10-liver1'])

case11['Promedio case11 lung'] = list(sts.mean(case11.loc[case11lung,'case11-lung1':'case11-lung3']) for case11lung in GenData3)
case11['Promedio case11 kidney'] = list(case11['case11-kidney1'])
case11['Promedio case11 bowel'] = list(case11['case11-bowel1'])

case12['Promedio case12 liver'] = list(case12['case12-liver1'])

casosletras = [caseA,caseB,caseC,caseD,
               caseF,caseE,caseG,caseH,
               caseI,caseJ]

for c,l in zip(casosletras,list('ABCDFEGHIJ')):
    c['Promedio case'+l+' lung'] = list(c['case'+l+'-lung-NYC'])
    
caseP1['Promedio caseP1 placenta'] = list(sts.mean(caseP1.loc[caseP1placenta]) for caseP1placenta in GenData3)

caseP2['Promedio caseP2 placenta'] = list(sts.mean(caseP2.loc[caseP2placenta]) for caseP2placenta in GenData3)

caseP3['Promedio caseP3 placenta'] = list(sts.mean(caseP3.loc[caseP3placenta]) for caseP3placenta in GenData3)

caseP4['Promedio caseP4 placenta'] = list(caseP4['caseP4-placenta1'])

#fold-change

case1['fold-change control mujer case1 lung'] = fold_change(controles['Control Mujer'],case1['Promedio case1 lung'])
case1['fold-change control mujer case1 heart'] = fold_change(controles['Control Mujer'],case1['Promedio case1 heart'])

case2['fold-change control hombre case2 lung'] = fold_change(controles['Control Hombre'],case2['Promedio case2 lung'])
case2['fold-change control hombre case2 jejunum'] = fold_change(controles['Control Hombre'],case2['Promedio case2 jejunum'])
case2['fold-change control hombre case2 heart'] = fold_change(controles['Control Hombre'],case2['Promedio case2 heart'])

case3['fold-change control hombre case3 lung'] = fold_change(controles['Control Hombre'],case3['Promedio case3 lung'])
case3['fold-change control hombre case3 heart'] = fold_change(controles['Control Hombre'],case3['Promedio case3 heart'])

case4['fold-change control mujer case4 lung'] = fold_change(controles['Control Mujer'],case4['Promedio case4 lung'])
case4['fold-change control mujer case4 heart'] = fold_change(controles['Control Mujer'],case4['Promedio case4 heart'])
case4['fold-change control mujer case4 liver'] = fold_change(controles['Control Mujer'],case4['Promedio case4 liver'])
case4['fold-change control mujer case4 kidney'] = list(case4['case4-kidney1'])
case4['fold-change control mujer case4 bowel'] = fold_change(controles['Control Mujer'],case4['Promedio case4 bowel'])

case5['fold-change control hombre case5 lung'] = fold_change(controles['Control Hombre'],case5['Promedio case5 lung'])
case5['fold-change control hombre case5 fat'] = fold_change(controles['Control Hombre'],case5['Promedio case5 fat'])
case5['fold-change control hombre case5 skin'] = fold_change(controles['Control Hombre'],case5['Promedio case5 skin'])
case5['fold-change control hombre case5 bowel'] = fold_change(controles['Control Hombre'],case5['Promedio case5 bowel'])
case5['fold-change control hombre case5 liver'] = fold_change(controles['Control Hombre'],case5['Promedio case5 liver'])
case5['fold-change control hombre case5 kidney'] = fold_change(controles['Control Hombre'],case5['Promedio case5 kidney'])
case5['fold-change control hombre case5 heart'] = fold_change(controles['Control Hombre'],case5['Promedio case5 heart'])
case5['fold-change control hombre case5 marrow'] = fold_change(controles['Control Hombre'],case5['Promedio case5 marrow'])

case6['fold-change control hombre case6 lung'] = fold_change(controles['Control Hombre'],case6['Promedio case6 lung'])

case7['fold-change control mujer case7 lung'] = fold_change(controles['Control Mujer'],case7['Promedio case7 lung'])

case8['fold-change control hombre case8 lung'] = fold_change(controles['Control Hombre'],case8['Promedio case8 lung'])
case8['fold-change control hombre case8 liver'] = fold_change(controles['Control Hombre'],case8['Promedio case8 liver'])
case8['fold-change control hombre case8 bowel'] = fold_change(controles['Control Hombre'],case8['Promedio case8 bowel'])
case8['fold-change control hombre case8 heart'] = fold_change(controles['Control Hombre'],case8['Promedio case8 heart'])

case9['fold-change control mujer case9 lung'] = fold_change(controles['Control Mujer'],case9['Promedio case9 lung'])

case10['fold-change control mujer case10 lung'] = fold_change(controles['Control Mujer'],case10['Promedio case10 lung'])
case10['fold-change control mujer case10 liver'] = fold_change(controles['Control Mujer'],case10['Promedio case10 liver'])

case11['fold-change control hombre case11 lung'] = fold_change(controles['Control Hombre'],case11['Promedio case11 lung'])
case11['fold-change control hombre case11 kidney'] = fold_change(controles['Control Hombre'],case11['Promedio case11 kidney'])
case11['fold-change control hombre case11 bowel'] = fold_change(controles['Control Hombre'],case11['Promedio case11 bowel'])

case12['fold-change control mujer case12 liver'] = fold_change(controles['Control Mujer'],case12['Promedio case12 liver'])

for c1,l1 in zip([caseA,caseB,caseC,caseF,caseE,caseG,caseH,caseI],list('ABCFEGHI')):
    c1['fold-change control hombre case'+l1+' lung'] = fold_change(controles['Control Hombre'],c1['Promedio case'+l1+' lung'])

caseD['fold-change control mujer caseD lung'] = fold_change(controles['Control Mujer'],caseD['Promedio caseD lung'])

caseJ['fold-change control mujer caseJ lung'] = fold_change(controles['Control Mujer'],caseJ['Promedio caseJ lung'])

caseP1['fold-change control mujer caseP1 placenta'] = fold_change(controles['Control Mujer'],caseP1['Promedio caseP1 placenta'])

caseP2['fold-change control mujer caseP2 placenta'] = fold_change(controles['Control Mujer'],caseP2['Promedio caseP2 placenta'])

caseP3['fold-change control hombre caseP3 placenta'] = fold_change(controles['Control Hombre'],caseP3['Promedio caseP3 placenta'])

caseP4['fold-change control hombre caseP4 placenta'] = fold_change(controles['Control Hombre'],caseP4['Promedio caseP4 placenta'])

#ANEXAMOS LOS fold-change A Data3 CON GENERO

Data3 = pd.concat([Data3,case1.loc[:,['fold-change control mujer case1 lung','fold-change control mujer case1 heart']]], axis = 1)
Data3 = pd.concat([Data3,case2.loc[:,['fold-change control hombre case2 lung','fold-change control hombre case2 jejunum','fold-change control hombre case2 heart']]], axis = 1)
Data3 = pd.concat([Data3,case3.loc[:,['fold-change control hombre case3 lung','fold-change control hombre case3 heart']]], axis = 1)
Data3 = pd.concat([Data3,case4.loc[:,'fold-change control mujer case4 lung':'fold-change control mujer case4 bowel']],axis = 1)
Data3 = pd.concat([Data3,case5.loc[:,'fold-change control hombre case5 lung':'fold-change control hombre case5 marrow']], axis = 1)
Data3 = pd.concat([Data3,case6.loc[:,'fold-change control hombre case6 lung']], axis = 1)
Data3 = pd.concat([Data3,case7.loc[:,'fold-change control mujer case7 lung']], axis = 1)
Data3 = pd.concat([Data3,case8.loc[:,'fold-change control hombre case8 lung':'fold-change control hombre case8 heart']],axis = 1)
Data3 = pd.concat([Data3,case9.loc[:,'fold-change control mujer case9 lung']],axis = 1)
Data3 = pd.concat([Data3,case10.loc[:,['fold-change control mujer case10 lung','fold-change control mujer case10 liver']]], axis = 1)
Data3 = pd.concat([Data3,case11.loc[:,'fold-change control hombre case11 lung':'fold-change control hombre case11 bowel']],axis = 1)
Data3 = pd.concat([Data3,case12.loc[:,'fold-change control mujer case12 liver']], axis  = 1)
for c2,l2 in zip([caseA,caseB,caseC,caseF,caseE,caseG,caseH,caseI],list('ABCFEGHI')):
    Data3 = pd.concat([Data3,c2.loc[:,'fold-change control hombre case'+l2+' lung']], axis = 1)
Data3 = pd.concat([Data3,caseD.loc[:,'fold-change control mujer caseD lung']],axis = 1)
Data3 = pd.concat([Data3,caseJ.loc[:,'fold-change control mujer caseJ lung']], axis = 1)
Data3 = pd.concat([Data3,caseP1.loc[:,'fold-change control mujer caseP1 placenta']], axis = 1)
Data3 = pd.concat([Data3,caseP2.loc[:,'fold-change control mujer caseP2 placenta']], axis = 1)
Data3 = pd.concat([Data3,caseP3.loc[:,'fold-change control hombre caseP3 placenta']], axis = 1)
Data3 = pd.concat([Data3,caseP4.loc[:,'fold-change control hombre caseP4 placenta']], axis = 1)

#FILTRAMOS LOS GENES DESEADOS CON GENERO Y SIN GENERO

fold_change_data_3 = Data3.loc[:,'fold-change control1 case1 lung':'fold-change control5 case5 marrow'].copy()
fold_change_data_3.loc[list(Genes.loc[:,0])].to_csv('fold_change_genes_1_General_GSE150316.csv')
fold_change_data_3.loc[list(Genes2.loc[:,0])].to_csv('fold_change_genes_2_General_GSE150316.csv')

Arreglo0 = pd.DataFrame({_:[nan] for _ in list(fold_change_data_3.columns)}, index = ['CR', 'DHNTP'])
fold_change_data_3 = pd.concat([fold_change_data_3,Arreglo0])

fold_change_data_3.loc[list(Genes3.loc[:,0])].to_csv('fold_change_genes_3_General_GSE150316.csv')

fold_change_data_3_genero = Data3.loc[:,'fold-change control mujer case1 lung':].copy()
fold_change_data_3_genero.loc[list(Genes.loc[:,0])].to_csv('fold_change_genes_1_Genero_GSE150316.csv')
fold_change_data_3_genero.loc[list(Genes2.loc[:,0])].to_csv('fold_change_genes_2_Genero_GSE150316.csv')

Arreglo1 = pd.DataFrame({_:[nan] for _ in list(fold_change_data_3_genero.columns)}, index = ['CR', 'DHNTP'])
fold_change_data_3_genero = pd.concat([fold_change_data_3_genero,Arreglo1])

fold_change_data_3_genero.loc[list(Genes3.loc[:,0])].to_csv('fold_change_genes_3_Genero_GSE150316.csv')

Data3.to_csv('Resultados_Analisis_GSE150316.csv')