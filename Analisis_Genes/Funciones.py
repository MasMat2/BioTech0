def fold_change (Columna_Mock,Columna_Experimento2):
    from numpy import nan
    from math import log
    L = []
    for Mock, Experimento in zip (list(Columna_Mock),list(Columna_Experimento2)):
        if (Experimento == 0 and Mock < 0) or (Experimento < 0 and Mock == 0):
            L.append('-inf')
        elif (Experimento == 0 and Mock > 0) or (Experimento > 0 and Mock == 0):
            L.append('inf')
        elif Experimento == 0 and Mock == 0:
            L.append(0)
        elif Experimento == nan or Mock == nan:
            L.append(nan)
        else:
            L.append(log(Experimento/Mock,2))
    return L

def Promedio (DF,Nombre):
    from statistics import mean
    DF['Promedio '+Nombre] = list(mean(DF.iloc[_]) for _ in range(len(DF.index)))