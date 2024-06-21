# Cargar librerías necesarias
library(ncdf4)
library(xts)
library(R6)

# Definición de la clase NeuralNetwork con dos capas ocultas
NeuralNetwork <- R6Class("NeuralNetwork",
                         public = list(
                           inputs = NULL,
                           hidden1 = NULL,
                           hidden2 = NULL,
                           outputs = NULL,
                           weights1 = NULL,
                           weights2 = NULL,
                           weights3 = NULL,
                           learning_rate = NULL,
                           layer1 = NULL,
                           layer2 = NULL,
                           output = NULL,
                           
                           initialize = function(inputs, hidden1, hidden2, outputs, learning_rate = 0.01) {
                             self$inputs <- inputs
                             self$hidden1 <- hidden1
                             self$hidden2 <- hidden2
                             self$outputs <- outputs
                             self$learning_rate <- learning_rate
                             
                             # Inicializar pesos aleatoriamente
                             self$weights1 <- matrix(runif(inputs * hidden1, min = -1, max = 1), nrow = inputs, ncol = hidden1)
                             self$weights2 <- matrix(runif(hidden1 * hidden2, min = -1, max = 1), nrow = hidden1, ncol = hidden2)
                             self$weights3 <- matrix(runif(hidden2 * outputs, min = -1, max = 1), nrow = hidden2, ncol = outputs)
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
                             self$layer2 <- self$sigmoid(self$layer1 %*% self$weights2)
                             self$output <- self$sigmoid(self$layer2 %*% self$weights3)  
                             return(self$output)
                           },
                           
                           backpropagate = function(X, y) {
                             output_error <- y - self$output
                             d_weights3 <- t(self$layer2) %*% (output_error * self$sigmoid_derivative(self$output))
                             
                             hidden2_error <- (output_error %*% t(self$weights3))
                             d_weights2 <- t(self$layer1) %*% (hidden2_error * self$sigmoid_derivative(self$layer2))
                             
                             hidden1_error <- (hidden2_error %*% t(self$weights2))
                             d_weights1 <- t(X) %*% (hidden1_error * self$sigmoid_derivative(self$layer1))
                             
                             # Actualizar pesos
                             self$weights1 <- self$weights1 + self$learning_rate * d_weights1
                             self$weights2 <- self$weights2 + self$learning_rate * d_weights2
                             self$weights3 <- self$weights3 + self$learning_rate * d_weights3
                           },
                           
                           train = function(X, y, epochs = 10000) {
                             for (i in 1:epochs) {
                               self$feedforward(X)
                               self$backpropagate(X, y)
                             }
                           }
                         )
)

# Función para normalizar los datos
normalize <- function(x, min_val = NULL, max_val = NULL) {
  if (is.null(min_val)) min_val <- min(x)
  if (is.null(max_val)) max_val <- max(x)
  return((x - min_val) / (max_val - min_val))
}

# Rutas de los archivos NetCDF
archivo_tmax <- "D:/Documentos/tesis/CR2MET_tmax_v2.0_mon_1979_2019_005deg.nc"
archivo_pr <- "D:/Documentos/tesis/CR2MET_pr_v2.0_mon_1979_2019_005deg.nc"

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

# Obtener los valores mínimos y máximos de los datos de entrenamiento
n_train <- round(0.7 * length(tmax_central))
tmax_min <- min(tmax_central[1:n_train])
tmax_max <- max(tmax_central[1:n_train])
pr_min <- min(pr_central[1:n_train])
pr_max <- max(pr_central[1:n_train])

# Normalizar los datos de entrenamiento y prueba utilizando los mismos valores mínimos y máximos
tmax_norm <- normalize(tmax_central, tmax_min, tmax_max)
pr_norm <- normalize(pr_central, pr_min, pr_max)

# Crear ventanas deslizantes para las series temporales
create_sliding_window <- function(data, window_size) {
  X <- NULL
  y <- NULL
  for (i in seq_len(length(data) - window_size)) {
    X <- rbind(X, data[i:(i + window_size - 1)])
    y <- c(y, data[i + window_size])
  }
  return(list(X = X, y = y))
}

# Definir el tamaño de la ventana
window_size <- 12

# Crear ventanas deslizantes para la temperatura máxima y la precipitación
windows_tmax <- create_sliding_window(tmax_norm, window_size)
windows_pr <- create_sliding_window(pr_norm, window_size)

# Crear el dataframe con las ventanas deslizantes
df <- data.frame(
  windows_tmax$X,
  precipitation = windows_pr$y
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
X_train <- as.matrix(datos_entrenamiento[, -ncol(datos_entrenamiento)])
y_train <- as.matrix(datos_entrenamiento[, ncol(datos_entrenamiento)])

X_test <- as.matrix(datos_prueba[, -ncol(datos_prueba)])
y_test <- as.matrix(datos_prueba[, ncol(datos_prueba)])

# Crear una instancia de la red neuronal con dos capas ocultas
nn <- NeuralNetwork$new(inputs = window_size, hidden1 = 50, hidden2 = 30, outputs = 1, learning_rate = 0.01)

# Entrenar la red neuronal
nn$train(X_train, y_train, epochs = 10000)

# Hacer predicciones sobre el conjunto de prueba
predicciones <- nn$feedforward(X_test)

# Comparar las predicciones con las precipitaciones reales
resultados <- data.frame(
  Index = seq_len(nrow(datos_prueba)),
  Real = y_test,
  Prediction = predicciones
)

# Graficar las precipitaciones reales y las predicciones en el conjunto de prueba
plot(resultados$Real, type = "l", col = "blue", xlab = "Índice", ylab = "Precipitación Normalizada", main = "Precipitaciones Reales vs. Predicciones (Conjunto de Prueba)")
lines(resultados$Prediction, col = "red")
legend("topright", legend = c("Precipitación Real", "Predicción"), col = c("blue", "red"), lty = 1)
