Este repositorio contiene los archivos generados a partir del análisis exploratorio del estudio ST000291 del Metabolomics Workbench, que consta de dos análisis experimentales: AN00064 y AN00065.

## Contenidos

- **datos_ST000291-AN000464.txt**: Matriz de datos del análisis AN00064.
- **datos_ST000291-AN000465.txt**: Matriz de datos del análisis AN00065.
- **metadatos_ST000291.md**: Metadatos del estudio, donde está la información sobre los experimentos, los metabolitos de ambos análisis y las muestras.
- **informe.pdf**: Informe que describe todo el proceso realizado, desde la selección del dataset, pasando por el análisis exploratorio del mismo, hasta la interpretación de los resultados.
- **informe.rmd**: Mismo contenido que el informe.pdf, pero aquí se puede ver como se han construido las tablas y los gráficos para el documento PDF.
- **analisis_exploratorio.R**: Código R utilizado tanto para la descarga del dataset como para el análisis exploratorio de este.
- **ST000291_datos_y_metadatos.Rda**: Objeto que contiene los datos y metadatos de los dos análisis, cada uno constituyendo un objeto de clase `SummarizedExperiment`.
