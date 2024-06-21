library(R6)

# Definición de la clase NeuralNetwork
NeuralNetwork <- R6Class("NeuralNetwork",
                         public = list(
                           inputs = NULL,
                           hidden = NULL,
                           outputs = NULL,
                           weights1 = NULL,
                           weights2 = NULL,
                           learning_rate = NULL,
                           layer1 = NULL,
                           output = NULL,
                           
                           initialize = function(inputs, hidden, outputs, learning_rate = 0.1) {
                             self$inputs <- inputs
                             self$hidden <- hidden
                             self$outputs <- outputs
                             self$learning_rate <- learning_rate
                             
                             # Inicializar pesos aleatoriamente
                             self$weights1 <- matrix(runif(inputs * hidden, min = -1, max = 1), nrow = inputs, ncol = hidden)
                             self$weights2 <- matrix(runif(hidden * outputs, min = -1, max = 1), nrow = hidden, ncol = outputs)
                           },
                           
                           sigmoid = function(x) {
                             result <- 1 / (1 + exp(-x))
                             return(result)
                           },
                           
                           sigmoid_derivative = function(x) {
                             result <- x * (1 - x)
                             return(result)
                           },
                           
                           feedforward = function(X) {
                             self$layer1 <- self$sigmoid(X %*% self$weights1)
                             # Aplicar sigmoide en la capa de salida
                             self$output <- self$sigmoid(self$layer1 %*% self$weights2)  
                             return(self$output)
                           },
                           
                           backpropagate = function(X, y) {
                             output_error <- y - self$output
                             d_weights2 <- t(self$layer1) %*% (output_error * self$sigmoid_derivative(self$output))
                             
                             hidden_error <- (output_error %*% t(self$weights2))
                             d_weights1 <- t(X) %*% (hidden_error * self$sigmoid_derivative(self$layer1))
                             
                             # Actualizar pesos
                             self$weights1 <- self$weights1 + self$learning_rate * d_weights1
                             self$weights2 <- self$weights2 + self$learning_rate * d_weights2
                           },
                           
                           train = function(X, y, epochs = 10000) {
                             for (i in 1:epochs) {
                               self$feedforward(X)
                               self$backpropagate(X, y)
                             }
                           }
                         )
)

# Cargar librerías necesarias
library(ncdf4)
library(xts)
library(R6)

# Función para normalizar los datos
normalize <- function(x) {
  min_val <- min(x)
  max_val <- max(x)
  return((x - min_val) / (max_val - min_val))
}

# Rutas de los archivos NetCDF
archivo_tmax <- "G:/scripts/CR2MET_tmax_v2.0_mon_1979_2019_005deg.nc"
archivo_pr <- "G:/scripts/CR2MET_pr_v2.0_mon_1979_2019_005deg.nc"

# Abrir los archivos NetCDF
nc_data_tmax <- nc_open(archivo_tmax)
nc_data_pr <- nc_open(archivo_pr)

# Obtener las dimensiones y variables de temperatura máxima y precipitación
lon <- ncvar_get(nc_data_tmax, "lon")
lat <- ncvar_get(nc_data_tmax, "lat")
tmax <- ncvar_get(nc_data_tmax, "tmax")
pr <- ncvar_get(nc_data_pr, "pr")

# Calcular el índice del valor central para lon y lat
indice_central_lon <- length(lon) %/% 2
indice_central_lat <- length(lat) %/% 2

# Obtener las temperaturas máximas y precipitaciones para la coordenada central
tmax_central <- tmax[indice_central_lon, indice_central_lat, ]
pr_central <- pr[indice_central_lon, indice_central_lat, ]

# Crear una serie temporal para temperatura máxima y precipitación
fecha_inicio <- as.Date("1978-12-15")
fechas <- seq(fecha_inicio, by = "months", length.out = length(tmax_central))
serie_temporal_tmax <- xts(tmax_central, order.by = fechas)
serie_temporal_pr <- xts(pr_central, order.by = fechas)

# Cerrar los archivos NetCDF
nc_close(nc_data_tmax)
nc_close(nc_data_pr)

# Normalizar los datos de temperatura máxima y precipitación
tmax_norm <- normalize(tmax_central)
pr_norm <- normalize(pr_central)

# Crear el dataframe con índices y temperatura para las entradas, y precipitación para las salidas
df <- data.frame(
  index = 1:length(tmax_norm),
  temperature = tmax_norm,
  precipitation = pr_norm
)

# Dividir los datos en conjuntos de entrenamiento y prueba (70% - 30%)
set.seed(123)  # Para reproducibilidad
n <- nrow(df)
n_train <- round(0.7 * n)

# Conjunto de entrenamiento
datos_entrenamiento <- df[1:n_train, ]

# Conjunto de prueba
datos_prueba <- df[(n_train + 1):n, ]

# Preparar los datos de entrada y salida para la red neuronal
X_train <- as.matrix(datos_entrenamiento[, c("index", "temperature")])
y_train <- as.matrix(datos_entrenamiento[, "precipitation"])

X_test <- as.matrix(datos_prueba[, c("index", "temperature")])
y_test <- as.matrix(datos_prueba[, "precipitation"])

# Crear una instancia de la red neuronal
nn <- NeuralNetwork$new(inputs = 2, hidden = 30, outputs = 1, learning_rate = 0.00001)

# Entrenar la red neuronal
nn$train(X_train, y_train, epochs = 1000)

# Hacer predicciones sobre el conjunto de prueba
predicciones <- nn$feedforward(X_test)

# Comparar las predicciones con las precipitaciones reales
resultados <- data.frame(
  Index = datos_prueba$index,
  Temperature = datos_prueba$temperature,
  Real = datos_prueba$precipitation,
  Prediction = predicciones
)

# Graficar las precipitaciones reales y las predicciones en el conjunto de prueba
plot(resultados$Real, type = "l", col = "blue", xlab = "Índice", ylab = "Precipitación Normalizada", main = "Precipitaciones Reales vs. Predicciones (Conjunto de Prueba)")
lines(resultados$Prediction, col = "red")
legend("topright", legend = c("Precipitación Real", "Predicción"), col = c("blue", "red"), lty = 1)

# Calcular errores
errores <- resultados$Real - resultados$Prediction

# Calcular el Error Cuadrático Medio (MSE)
mse <- mean(errores^2)

# Calcular el Error Absoluto Medio (MAE)
mae <- mean(abs(errores))

# Calcular el Coeficiente de Determinación (R^2)
ss_res <- sum(errores^2)
ss_tot <- sum((resultados$Real - mean(resultados$Real))^2)
r_cuadrado <- 1 - ss_res/ss_tot


# Calcular PRESS
press <- sum((y_test - nn$feedforward(X_test))^2)

# Calcular SSTotal
sstotal <- sum((y_test - mean(y_train))^2)

# Calcular R^2(pred)
r2_pred <- 1 - (press / sstotal)

# Reportar las métricas
cat("Error Cuadrático Medio (MSE):", mse, "\n")
cat("Error Absoluto Medio (MAE):", mae, "\n")
cat("Coeficiente de Determinación (R^2):", r_cuadrado, "\n")
cat("Coeficiente de Determinación de Predicción (R^2 pred):", r2_pred, "\n")

