1. Los archivos que inician con nombre 'Analisis_GSE...' y son archivos de Python se utilizan
para obtener los fold-change de los respectivos GSE.

	1.1 La obtención de los datos GSE correspondientes al archivo se encuentra en la carpeta
	    del mismo nombre, es decir, para el análisis de GSE179277 usamos el Python script 
            'Analisis_GSE179277.py' y la data se almacena en la carpeta 'GSE179277'.

2. El archivo 'Funciones.py' es un archivo de Python el cual contiene funciones que se utilizan en 
los archivos 'Analisis_GSE...', por ejemplo : 

	*Cálculo de fold-change
	*Obtención de la media

3. Los resultados fold-change se almacenan en la carpeta 'fold-change'

	3.1 El análisis general se almacena en la carpeta 'Resultados'

4. Los archivos 'Genes_descripcion_codigo' almacena dos manera distintas de la descripción de un mismo gen,
esto se tiene ya que las bases de datos no manejan un estandar en la descripción de un gen.

5. La carpeta 'Genes' almacena los genes deseados a analizar.