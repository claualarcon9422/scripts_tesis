
n <- 120
ultima_fila <- X_test_con_errores_b[nrow(X_test_con_errores_b), ]
nuevos_datos <- matrix(nrow = 0, ncol = ncol(X_test_con_errores_b))  # Inicializar matriz vacía

for(i in 1:n){
  #print(ultima_fila)
  nueva_fila <- ultima_fila[3:(length(ultima_fila)-3)]
  #print(nueva_fila)
  prediccion <- nn$feedforward(ultima_fila)
  #print(prediccion)
  lag1 <- prediccion - nueva_fila[length(ultima_fila) - 11]
  lag2 <- prediccion - nueva_fila[length(ultima_fila) - 23]
  error_prediccion <- runif(1, min(errores), max(errores))
  #print(lag1)
  #print(lag2)
  nueva_fila <- c(nueva_fila,prediccion)
  #print(nueva_fila)
  nueva_fila <- c("bias"=mejor_sesgo, nueva_fila, "lag1"=lag1, "lag2"=lag2,"error"=error_prediccion)
  #print(length(nueva_fila))
  nuevos_datos <- rbind(nuevos_datos, nueva_fila)  # Agregar nueva fila a nuevos_datos
  #print(nueva_fila)
  #print(nuevos_datos)
  ultima_fila <- nueva_fila
  #print(ultima_fila)
}

#print(nn$feedforward(nuevos_datos))

# Hacer predicciones sobre los nuevos datos
predicciones_nuevas <- nn$feedforward(nuevos_datos)

# Crear un vector de fechas para las nuevas predicciones
fechas_nuevas <- seq(fechas_prueba[length(fechas_prueba)], by = "months", length.out = n)

# Agregar las nuevas predicciones al gráfico existente
lines(
  fechas_nuevas,
  predicciones_nuevas,
  col = "orange",  # Color de las nuevas predicciones
  lwd = 2,        # Grosor de la línea
  type = "l"      # Tipo de línea
)

# Agregar una leyenda actualizada
legend(
  "topleft",
  legend = c("Real", "Ajuste", "Predicción (Prueba)", "Predicción (Nuevos Datos)"),
  col = c("blue", "green", "red", "orange"),
  lty = 1,
  xpd = TRUE,
  inset = c(-0.25, 0)
)