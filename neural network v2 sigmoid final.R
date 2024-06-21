library(ncdf4)
library(xts)
library(R6)

NeuralNetwork <- R6Class(
  "NeuralNetwork",
  public = list(
    inputs = NULL,
    hidden1 = NULL,
    hidden2 = NULL,
    outputs = NULL,
    weights1 = NULL,
    weights2 = NULL,
    weights3 = NULL,
    bias1 = NULL,
    bias2 = NULL,
    bias3 = NULL,
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
      
      # Inicializar bias aleatoriamente
      self$bias1 <- runif(hidden1, min = -1, max = 1)
      self$bias2 <- runif(hidden2, min = -1, max = 1)
      self$bias3 <- runif(outputs, min = -1, max = 1)
    },
    
    # Función de activación
    sigmoid = function(x) {
      result <- 1 / (1 + exp(-x))
      return(result)
    },
    
    sigmoid_derivative = function(x) {
      result <- x * (1 - x)
      return(result)
    },
    
    feedforward = function(X) {
      self$layer1 <- self$sigmoid(X %*% self$weights1 + matrix(self$bias1, nrow = nrow(X), ncol = length(self$bias1), byrow = TRUE))
      self$layer2 <- self$sigmoid(self$layer1 %*% self$weights2 + matrix(self$bias2, nrow = nrow(self$layer1), ncol = length(self$bias2), byrow = TRUE))
      self$output <- self$sigmoid(self$layer2 %*% self$weights3 + matrix(self$bias3, nrow = nrow(self$layer2), ncol = length(self$bias3), byrow = TRUE))
      return(self$output)
    },
    
    # Ajustar pesos y bias en base al error
    backpropagate = function(X, y) {
      output_error <- y - self$output
      d_output <- output_error * self$sigmoid_derivative(self$output)
      d_weights3 <- t(self$layer2) %*% d_output
      d_bias3 <- colSums(d_output)
      
      hidden2_error <- d_output %*% t(self$weights3)
      d_hidden2 <- hidden2_error * self$sigmoid_derivative(self$layer2)
      d_weights2 <- t(self$layer1) %*% d_hidden2
      d_bias2 <- colSums(d_hidden2)
      
      hidden1_error <- d_hidden2 %*% t(self$weights2)
      d_hidden1 <- hidden1_error * self$sigmoid_derivative(self$layer1)
      d_weights1 <- t(X) %*% d_hidden1
      d_bias1 <- colSums(d_hidden1)
      
      # Actualizar pesos y bias
      self$weights1 <- self$weights1 + self$learning_rate * d_weights1
      self$weights2 <- self$weights2 + self$learning_rate * d_weights2
      self$weights3 <- self$weights3 + self$learning_rate * d_weights3
      
      self$bias1 <- self$bias1 + self$learning_rate * d_bias1
      self$bias2 <- self$bias2 + self$learning_rate * d_bias2
      self$bias3 <- self$bias3 + self$learning_rate * d_bias3
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
normalize <- function(x,
                      min_val = NULL,
                      max_val = NULL) {
  if (is.null(min_val))
    min_val <- min(x)
  if (is.null(max_val))
    max_val <- max(x)
  return((x - min_val) / (max_val - min_val))
}

#Cargar datos
archivo_tmax <- "CR2MET_tmax_v2.0_mon_1979_2019_005deg.nc"
archivo_pr <- "CR2MET_pr_v2.0_mon_1979_2019_005deg.nc"

#Abrir archivo series
nc_data_tmax <- nc_open(archivo_tmax)
nc_data_pr <- nc_open(archivo_pr)

# Obtener las dimensiones y variables de temperatura máxima y precipitación
lon <- ncvar_get(nc_data_tmax, "lon")
lat <- ncvar_get(nc_data_tmax, "lat")
tmax <- ncvar_get(nc_data_tmax, "tmax")
pr <- ncvar_get(nc_data_pr, "pr")

# Obtener coordenadas centrales del mejor país de Chile
indice_central_lon <- length(lon) %/% 2
indice_central_lat <- length(lat) %/% 2

# Obtener las temperaturas máximas y precipitaciones para la coordenada central
tmax_central <- tmax[indice_central_lon, indice_central_lat, ]
pr_central <- pr[indice_central_lon, indice_central_lat, ]

# Crear una serie temporal para temperatura máxima y precipitación
fecha_inicio <- as.Date("1978-12-15")
fechas <- seq(fecha_inicio,
              by = "months",
              length.out = length(tmax_central))
serie_temporal_tmax <- xts(tmax_central, order.by = fechas)
serie_temporal_pr <- xts(pr_central, order.by = fechas)

# Cerrar los archivos
nc_close(nc_data_tmax)
nc_close(nc_data_pr)

# Obtener los valores mínimos y máximos de los datos de entrenamiento para normalizar
n_train <- round(0.7 * length(tmax_central))
tmax_min <- min(tmax_central[1:n_train])
tmax_max <- max(tmax_central[1:n_train])
pr_min <- min(pr_central[1:n_train])
pr_max <- max(pr_central[1:n_train])

# Normalizar los datos de entrenamiento y prueba utilizando los mismos valores mínimos y máximos
tmax_norm <- normalize(tmax_central, tmax_min, tmax_max)
pr_norm <- normalize(pr_central, pr_min, pr_max)

# Función para crear ventanas deslizantes
create_sliding_windows <- function(data, windows_size) {
  X <- NULL
  y <- NULL
  for (i in seq_len(length(data) - windows_size)) {
    X <- rbind(X, data[i:(i + windows_size - 1)])
    y <- c(y, data[i + windows_size])
  }
  return(list(X = X, y = y))
}

# Definir el tamaño de la ventana
windows_size <- round(0.25 * length(tmax_central))

# Crear vectores de temperaturas pasadas y precipitacion objetivo
windows_tmax <- create_sliding_windows(tmax_norm, windows_size)
windows_pr <- create_sliding_windows(pr_norm, windows_size)

# Extraer las matrices X e y directamente de los resultados de la función create_sliding_windows
X <- windows_tmax$X
y <- windows_tmax$y

#print(y)

# Dividir los datos en conjuntos de entrenamiento y prueba (70% - 30%)
set.seed(123)
n <- nrow(X)
n_train <- round(0.7 * n)

# Conjunto de entrenamiento
X_train <- X[1:n_train, ]
y_train <- y[1:n_train]

# Conjunto de prueba
X_test <- X[(n_train + 1):n, ]
y_test <- y[(n_train + 1):n]

#Añadir columnas para lags
X_train <- cbind(X_train,
                 lag1 = rep(NA, nrow(X_train)),
                 lag2 = rep(NA, nrow(X_train)))
X_test <- cbind(X_test, lag1 = rep(NA, nrow(X_test)), lag2 = rep(NA, nrow(X_test)))

#Primeros 12 lags son los mismos valores de y (temp reales)
X_train[, "lag1"] <- y_train
X_train[, "lag2"] <- y_train

X_test[, "lag1"] <- y_test
X_test[, "lag2"] <- y_test

# Agregar lag1 y lag2
for (i in 13:nrow(X_train)) {
  X_train[i, "lag1"] <- y_train[i] - y_train[i - 12]
}

for (i in 25:nrow(X_train)) {
  X_train[i, "lag2"] <- y_train[i] - y_train[i - 24]
}


for (i in 13:nrow(X_test)) {
  X_test[i, "lag1"] <- y_test[i] - y_test[i - 12]
}

for (i in 25:nrow(X_test)) {
  X_test[i, "lag2"] <- y_test[i] - y_test[i - 24]
}



#print(X_train[, "lag1", drop = FALSE])

# Combinar X_train[, "lag1"] y y_train
combined_data <- cbind(X_train[, "lag1"], X_train[, "lag2"], y_train)

# Imprimir lado a lado en la consola
#print(combined_data)

# Crear una lista de diferentes valores de sesgo para probar
#bias_values <- seq(from = -1, to = 1, by = 0.1)
bias_value <- 1

# Inicializar variables para guardar el mejor sesgo, mejor R² y el mejor modelo
mejor_sesgo <- NULL
mejor_r2 <- -Inf
mejor_nn <- NULL
mejor_X_train_b <- NULL
mejor_X_test_b <- NULL


# Bucle para probar diferentes valores de sesgo
#for (bias_value in bias_values) {
# Añadir una columna de sesgo a X_train y X_test
#X_train_b <- cbind(bias_value, X_train)
#X_test_b <- cbind(bias_value, X_test)
global_hidden1 = 30
global_hidden2 = 30

# Crear y entrenar la red neuronal
nn <- NeuralNetwork$new(
  inputs = ncol(X_train),
  global_hidden1,
  global_hidden2,
  outputs = 1,
  learning_rate = 0.01
)
nn$train(X_train, y_train, epochs = 1000)

# Hacer predicciones sobre el conjunto de entrenamiento
predicciones_entrenamiento <- nn$feedforward(X_train)

# Calcular el R² para las predicciones
errores <- y_train - predicciones_entrenamiento
ss_res <- sum(errores ^ 2)
ss_tot <- sum((y_train - mean(y_train)) ^ 2)
r_cuadrado <- 1 - (ss_res / ss_tot)

# Guardar el mejor sesgo y R²
if (r_cuadrado > mejor_r2) {
  mejor_r2 <- r_cuadrado
  mejor_sesgo <- bias_value
  mejor_X_train <- X_train
  mejor_X_test_b <- X_test_b
  mejor_nn <- nn
}
#}

nn <- mejor_nn
X_train <- mejor_X_train
X_test_b <- mejor_X_test_b

# Calcular errores
errores <- y_test - nn$feedforward(X_test)

# Calcular el Error Cuadrático Medio (MSE)
mse <- mean(errores ^ 2)
# Calcular el Error Absoluto Medio (MAE)
mae <- mean(abs(errores))
# Calcular el R² para las predicciones de entrenamiento
errores_entrenamiento <- y_train - nn$feedforward(X_train)
ss_res_entrenamiento <- sum(errores_entrenamiento ^ 2)
ss_tot_entrenamiento <- sum((y_train - mean(y_train)) ^ 2)
r_cuadrado_entrenamiento <- 1 - (ss_res_entrenamiento / ss_tot_entrenamiento)

# Calcular el R² para las predicciones de prueba
errores_prueba<- y_test - nn$feedforward(X_test)
ss_res_prueba <- sum(errores_prueba ^ 2)
ss_tot_prueba <- sum((y_test - mean(y_test)) ^ 2)
r_cuadrado_prueba <- 1 - (ss_res_prueba / ss_tot_prueba)


# Imprimir el mejor sesgo y R²
cat("Mejor Sesgo:", mejor_sesgo, "\n")
#mse
cat("Error Cuadrático Medio:", mse, "\n")
#mae
cat("Error Absoluto Medio:", mae, "\n")
# Reportar el coeficiente de determinación (R²) para los datos de entrenamiento
cat(
  "Coeficiente de Determinación (R^2) para los datos de entrenamiento:",
  r_cuadrado_entrenamiento,
  "\n"
)
cat("Coeficiente de Predicción R²:", r_cuadrado_prueba, "\n")

#print(X_test_b)

# Crear vector de errores
errores <- y_train - nn$feedforward(X_train)
errores_test <- runif(nrow(X_test), min(errores), max(errores))

#Añadir errores al conjunto de entrenamiento
X_train_con_errores <- cbind(X_train, errores)
#ruido blanco para conjunto de prueba
X_test_con_errores <- cbind(X_test, errores_test)

#print(X_train_con_errores)
#print(X_test_con_errores)


#Entrenar nuevamente la red con los errores calculados
# Crear una lista de diferentes valores de sesgo para probar
#bias_values <- seq(from = -1, to = 1, by = 0.1)
bias_value <- 1

# Inicializar variables para guardar el mejor sesgo, mejor R² y el mejor modelo
mejor_sesgo <- NULL
mejor_r2 <- -Inf
mejor_nn <- NULL
mejor_X_train_con_errores <- NULL
mejor_X_test_con_errores <- NULL

# Bucle para probar diferentes valores de sesgo
#for (bias_value in bias_values) {
#X_train_con_errores_b <- cbind(bias_value, X_train_con_errores)
#X_test_con_errores_b <- cbind(bias_value, X_test_con_errores)

# Crear y entrenar la red neuronal
nn <- NeuralNetwork$new(
  inputs = ncol(X_train_con_errores),
  global_hidden1,
  global_hidden2,
  outputs = 1,
  learning_rate = 0.01
)
nn$train(X_train_con_errores, y_train, epochs = 1000)

# Hacer predicciones sobre el conjunto de prueba
predicciones_entrenamiento <- nn$feedforward(X_train_con_errores)

# Calcular el R² para las predicciones
errores <- y_train - predicciones_entrenamiento
ss_res <- sum(errores ^ 2)
ss_tot <- sum((y_train - mean(y_train)) ^ 2)
r_cuadrado <- 1 - (ss_res / ss_tot)

# Guardar el mejor sesgo y R²
if (r_cuadrado > mejor_r2) {
  mejor_r2 <- r_cuadrado
  mejor_sesgo <- bias_value
  mejor_X_train_con_errores <- X_train_con_errores
  mejor_X_test_con_errores <- X_test_con_errores
  mejor_nn <- nn
}
#}

nn <- mejor_nn
#X_train_con_errores <- mejor_X_train_con_errores
#X_test_con_errores <- mejor_X_test_con_errores


# Calcular errores
errores <- y_test - nn$feedforward(X_test_con_errores)
# Calcular el Error Cuadrático Medio (MSE)
mse <- mean(errores ^ 2)
# Calcular el Error Absoluto Medio (MAE)
mae <- mean(abs(errores))
# Calcular el R² para las predicciones de entrenamiento
errores_entrenamiento <- y_train - nn$feedforward(X_train_con_errores)
ss_res_entrenamiento <- sum(errores_entrenamiento ^ 2)
ss_tot_entrenamiento <- sum((y_train - mean(y_train)) ^ 2)
r_cuadrado_entrenamiento <- 1 - (ss_res_entrenamiento / ss_tot_entrenamiento)

# Calcular el R² para las predicciones de prueba
errores_prueba<- y_test - nn$feedforward(X_test_con_errores)
ss_res_prueba <- sum(errores_prueba ^ 2)
ss_tot_prueba <- sum((y_test - mean(y_test)) ^ 2)
r_cuadrado_prueba <- 1 - (ss_res_prueba / ss_tot_prueba)


# Imprimir el mejor sesgo y R²
cat("Mejor Sesgo:", mejor_sesgo, "\n")
#mse
cat("Error Cuadrático Medio:", mse, "\n")
#mae
cat("Error Absoluto Medio:", mae, "\n")
# Reportar el coeficiente de determinación (R²) para los datos de entrenamiento
cat(
  "Coeficiente de Determinación (R^2) para los datos de entrenamiento:",
  r_cuadrado_entrenamiento,
  "\n"
)
cat("Coeficiente de Predicción R²:", r_cuadrado_prueba, "\n")

# Hacer predicciones sobre el conjunto de prueba
predicciones_prueba <- nn$feedforward(X_test_con_errores)


# Crear un vector con NA para el tamaño total de los datos
predicciones_total <- rep(NA, n)

# Colocar las predicciones del conjunto de prueba en el vector total de predicciones
predicciones_total[(n_train + 1):n] <- predicciones_prueba


#print(predicciones_total)

# plot(df$temperature, type = "l", col = "blue", lwd=2, xlab = "Index", ylab = "Precipitación Normalizada", main = "Precipitaciones Reales vs. Predicciones")
# lines(nn$feedforward(X_train), col = "red",lwd=2)
# lines(predicciones_total, col = "green",lwd=2)
#
# legend("topright", legend = c("Precipitación Real", "Ajuste de Entrenamiento (70%)", "Predicciones (30% de Prueba)"), col = c("blue", "red", "green"), lty = 1)

# Plotear con las fechas ajustadas
par(mar = c(5, 4, 4, 9) + 0.1)
fecha_inicio <- fecha_inicio + 123
print(fecha_inicio)
fechas_ajustadas <- seq(fecha_inicio, by = "months", length.out = length(y))
print(length(y))

plot(
  fechas_ajustadas,
  y,
  type = "l",
  col = "blue",
  lwd = 1,
  xlab = "Fecha",
  ylab = "Temperaturas Máximas Normalizada",
  main = "Temperaturas máximas Reales vs. Predicciones"
)

# Obtener las predicciones para el conjunto de entrenamiento y prueba
predicciones_entrenamiento <- nn$feedforward(X_train_con_errores)
predicciones_prueba <- predicciones_total[(n_train + 1):n]


# Crear un vector de fechas para el conjunto de entrenamiento
fechas_entrenamiento <- fechas_ajustadas[1:n_train]

# Crear un vector de fechas para el conjunto de prueba
fechas_prueba <- fechas_ajustadas[(n_train + 1):n]

# Agregar líneas para el ajuste de entrenamiento y las predicciones
lines(
  fechas_entrenamiento,
  predicciones_entrenamiento,
  col = rgb(0, 128, 0, maxColorValue = 255),
  lwd = 2
)
lines(
  fechas_prueba,
  predicciones_prueba,
  col = rgb(255, 0, 0, maxColorValue = 255),
  lwd = 2
)

# Agregar la leyenda
legend(
  "topright",
  legend = c("Real", "Ajuste", "Predicción"),
  col = c("blue", "green", "red"),
  lty = 1,
  xpd = TRUE,
  inset = c(-0.25, 0)
)

# Calcular los residuos
residuos <- y_test - predicciones_prueba

# Calcular los residuos estandarizados
residuos_estandarizados <- residuos / sd(residuos)


par(mfrow=c(3,1), mar=c(4,4,4,1)+.1)
# Graficar los residuos estandarizados
plot(residuos_estandarizados, type = "h", col = "black", lwd = 1, xlab = "Índice", ylab = "Residuos Estandarizados", main = "Residuos Estandarizados")
abline(h = 0, col = "black", lty = 1)
# Calcular la función de autocorrelación de los residuos
acf_residuos <- acf(residuos, main = "ACF de Residuos", lag.max = 20)

# Calcular la estadística de prueba Q_k para cada lag k
Q_statistic <- acf_residuos$acf[-1]^2 * length(residuos)

# Calcular los valores críticos para el test de Ljung-Box
valores_criticos <- qchisq(0.95, df = acf_residuos$lag[-1])

# Calcular los p-values
p_values <- 1 - pchisq(Q_statistic, df = acf_residuos$lag[-1])

plot(p_values, xlab = "Lag", ylab = "p-value", 
     main = "p-values del Test de Ljung-Box",
     type = "p")  # Cambiado de "o" a "p"
abline(h = 0.05, col = "blue", lty = 2)  # Agrega una línea en 0.05 para referencia

# Añadir etiquetas al eje x
axis(1, at = seq_along(p_values), labels = seq_along(p_values))


#print(residuos_estandarizados)
#print(acf_residuos)
#print(p_values)

