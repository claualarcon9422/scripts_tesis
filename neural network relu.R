
library(keras)

# Generar datos
set.seed(42)
x <- seq(-2 * pi, 2 * pi, length.out = 1000)
y <- sin(x)

# Dividir los datos en conjuntos de entrenamiento y prueba
set.seed(42)
train_indices <- sample(1:length(x), size = 0.8 * length(x))
x_train <- x[train_indices]
y_train <- y[train_indices]
x_test <- x[-train_indices]
y_test <- y[-train_indices]

# Normalizar los datos
x_train <- array_reshape(x_train, c(length(x_train), 1))
x_test <- array_reshape(x_test, c(length(x_test), 1))

# Definir el modelo de red neuronal
model <- keras_model_sequential()

model %>%
  layer_dense(units = 64, activation = 'relu', input_shape = c(1)) %>%
  layer_dense(units = 64, activation = 'relu') %>%
  layer_dense(units = 1)

# Compilar el modelo
model %>% compile(
  loss = 'mse',
  optimizer = optimizer_adam(),
  metrics = c('mae')
)

# Entrenar el modelo
history <- model %>% fit(
  x_train, y_train,
  epochs = 100,
  batch_size = 32,
  validation_split = 0.2
)

# Predecir con los datos de prueba
y_pred <- model %>% predict(x_test)

# Trazar la función real y la predicción de la red
plot(x_test, y_test, type = "l", col = "blue", lwd = 2, 
     main = "Función Seno Real vs Predicción de la Red Neuronal",
     xlab = "x", ylab = "y")
lines(x_test, y_pred, col = "red", lwd = 2)
legend("topright", legend = c("Función Real", "Predicción"),
       col = c("blue", "red"), lwd = 2)

# Evaluar el modelo
score <- model %>% evaluate(x_test, y_test)
cat('Mean Absolute Error en datos de prueba:', score$mae, '\n')

