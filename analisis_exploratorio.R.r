# Cargamos el paquete de Metabolomic Workbech
library(metabolomicsWorkbenchR)
library(SummarizedExperiment)
library(POMA)
library(ggtext)
library(ggplot2)

# Descargamos el dataset seleccionado como output un objeto SummarizedExperiment
SE <- do_query(
  context = "study",
  input_item = "study_id",
  input_value = "ST000291",
  output_item = "SummarizedExperiment"
)

#### Análisis exploratorio de la estructura del dataset
# 1. Mostramos qué contiene el objeto
print(SE)

# Vemos que contiene 2 análisis (ambos son objetos SummarizedExperiment):
# AN000464 y AN000465. Haremos la misma exploración para ambos

# 2. Mostramos los metadatos de cada análisis
metadata_AN000464 <- metadata(SE$AN000464)
print(metadata_AN000464)

metadata_AN000465 <- metadata(SE$AN000465)
print(metadata_AN000465)

# Creamos una tabla con los metadatos de cada análisis
metadata_exp <- data.frame(matrix(rep(NA, 3 * length(metadata_AN000464)), ncol = 3))
for (i in 1:length(metadata_AN000464)) {
  metadata_exp[i, 1] <- names(metadata_AN000464)[i]
  metadata_exp[i, 2] <- metadata_AN000464[i]
  metadata_exp[i, 3] <- metadata_AN000465[i]
}
colnames(metadata_exp) <- c("Campo", "AN000464", "AN000465")

# 3. Mostramos las matrices de datos
data_AN000464 <- assay(SE$AN000464)
print(data_AN000464)

data_AN000465 <- assay(SE$AN000465)
print(data_AN000465)

# 4. Mostramos la información de los metabolitos (filas de la
# matriz de datos)
metabolite_AN000464 <- rowData(SE$AN000464)
print(metabolite_AN000464)

metabolite_AN000465 <- rowData(SE$AN000465)
print(metabolite_AN000465)

# 5. Mostramos la información de las muestras (columnas de la
# matriz de datos) que es la misma para los dos análisis, por
# lo que solo es necesario inspeccionar una
muestras <- colData(SE$AN000465)
print(muestras)

### Análisis descriptivo del dataset
# Número de muestras por condición experimental
muestras_condicion <- table(muestras$Treatment_)
print(muestras_condicion)

# Número de valores faltantes, por muestra (columna), en la matriz
# de datos
NAs_table <- data.frame(
  AN000464 = apply(as.matrix(data_AN000464), 2, function(x) sum(is.na(x))),
  AN000465 = apply(as.matrix(data_AN000465), 2, function(x) sum(is.na(x)))
)
print(NAs_table)

# Ante la presencia de valores faltantes realizaremos su inputación
# utilizando POMA
set.seed(123)
AN000464_no_NAs <- PomaImpute(SE$AN000464)
AN000465_no_NAs <- PomaImpute(SE$AN000465)

# Resumen numérico univariante de cada muestra en la matriz
# de datos (media, mediana, desviación estándar y coeficiente
# de variación)
data_AN000464_no_NAs <- as.matrix(assay(AN000464_no_NAs))
data_AN000465_no_NAs <- as.matrix(assay(AN000465_no_NAs))

table_AN000464 <- data.frame(
  Media = apply(data_AN000464_no_NAs, 2, function(x) mean(x)),
  Mediana = apply(data_AN000464_no_NAs, 2, function(x) median(x)),
  SD = apply(data_AN000464_no_NAs, 2, function(x) sd(x)),
  VC = apply(data_AN000464_no_NAs, 2, function(x) sd(x) / mean(x))
)
print(table_AN000464)

table_AN000465 <- data.frame(
  Media = apply(data_AN000465_no_NAs, 2, function(x) mean(x)),
  Mediana = apply(data_AN000465_no_NAs, 2, function(x) median(x)),
  SD = apply(data_AN000465_no_NAs, 2, function(x) sd(x)),
  VC = apply(data_AN000465_no_NAs, 2, function(x) sd(x) / mean(x))
)
print(table_AN000465)

# A este tipo de datos es conveniente aplicarles una normalización
AN000464_no_NAs_norm <- PomaNorm(AN000464_no_NAs, method = "log_scaling")
AN000465_no_NAs_norm <- PomaNorm(AN000465_no_NAs, method = "log_scaling")

# Box plot y density plot de las muestras, usando POMA, para ver el efecto de
# la normalización
PomaBoxplots(AN000464_no_NAs) +
  ggplot2::ggtitle("Sin normalizar")
PomaBoxplots(AN000465_no_NAs) +
  ggplot2::ggtitle("Sin normalizar")
PomaBoxplots(AN000464_no_NAs_norm) +
  ggplot2::ggtitle("Con normalización")
PomaBoxplots(AN000465_no_NAs_norm) +
  ggplot2::ggtitle("Con normalización")

PomaDensity(AN000464_no_NAs) +
  ggplot2::ggtitle("Sin normalizar")
PomaDensity(AN000465_no_NAs) +
  ggplot2::ggtitle("Sin normalizar")
PomaDensity(AN000464_no_NAs_norm) +
  ggplot2::ggtitle("Con normalización")
PomaDensity(AN000465_no_NAs_norm) +
  ggplot2::ggtitle("Con normalización")

# Detección de outliers en nuestros datos utilizando la función
# PomaOutliers
outliers_AN000464 <- PomaOutliers(AN000464_no_NAs, outcome = "Treatment_")
outliers_AN000464$polygon_plot
outliers_AN000464$outliers

outliers_AN000465 <- PomaOutliers(AN000465_no_NAs, outcome = "Treatment_")
outliers_AN000465$polygon_plot
outliers_AN00046$outliers

# Obervar que tras normalizar los datos, las muestras antes consideradas
# outliers ya no lo son
outliers_AN000464_norm <- PomaOutliers(AN000464_no_NAs_norm, outcome = "Treatment_")
outliers_AN000464_norm$polygon_plot
outliers_AN000464_norm$outliers

outliers_AN000465_norm <- PomaOutliers(AN000465_no_NAs_norm, outcome = "Treatment_")
outliers_AN000465_norm$polygon_plot
outliers_AN00046_norm$outliers

# Finalmente, realizamos un PCA para visualizar la estructura de los datos
# normalizados y comprobar la existencia de posibles patrones
PCA_AN000464 <- prcomp(t(assay(AN000464_no_NAs_norm)), center = FALSE)

# Vamos a calcular el porcentaje de varianza explicada por cada componente
# para añadirlo al gráfico
var_exp_1 <- round(100 * PCA_AN000464$sdev^2 / sum(PCA_AN000464$sdev^2), 2)

# Hacemos el gráfico
PCA_AN000464_df <- as.data.frame(PCA_AN000464$x)
PCA_AN000464_df$group <- colData(AN000464_no_NAs_norm)$Treatment_
ggplot(PCA_AN000464_df, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    x = paste0("PC1", "(", var_exp_1[1], "%", ")"),
    y = paste0("PC2", "(", var_exp_2[2], "%", ")")
  ) +
  theme_classic()

# Lo mismo para AN000465
PCA_AN000465 <- prcomp(t(assay(AN000465_no_NAs_norm)), center = FALSE)

var_exp_2 <- round(100 * PCA_AN000465$sdev^2 / sum(PCA_AN000465$sdev^2), 2)

PCA_AN000465_df <- as.data.frame(PCA_AN000465$x)
PCA_AN000465_df$group <- colData(AN000465_no_NAs_norm)$Treatment_
ggplot(PCA_AN000465_df, aes(x = PC1, y = PC2, color = group, shape = group)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    x = paste0("PC1", "(", var_exp_2[1], "%", ")"),
    y = paste0("PC2", "(", var_exp_2[2], "%", ")")
  ) +
  theme_classic()
