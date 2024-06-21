library(ncdf4)
library(dplyr)  
library(ggplot2)
library(grid)
library(forecast)  # Para suavizado exponencial
library(xts)
library(tseries)

# Ruta del archivo NetCDF
archivo <- "CR2MET_tmax_v2.0_mon_1979_2019_005deg.nc"

# Abrir el archivo NetCDF
nc_data <- nc_open(archivo)

# Obtener las dimensiones y variables
lon <- ncvar_get(nc_data, "lon")
lat <- ncvar_get(nc_data, "lat")
tmax <- ncvar_get(nc_data, "tmax")

# Calcular el índice del valor central para lon y lat
indice_central_lon <- length(lon) %/% 2
indice_central_lat <- length(lat) %/% 2

# Obtener las temperaturas máximas para la coordenada central
tmax_central <- tmax[indice_central_lon, indice_central_lat, ]

# Crear una serie temporal
fecha_inicio <- as.Date("1978-12-15")
fechas <- seq(fecha_inicio, by = "months", length.out = length(tmax_central))
serie_temporal <- xts(tmax_central, order.by = fechas)

# Cerrar el archivo NetCDF
nc_close(nc_data)

# Ajustar modelo ARIMA manualmente
modelo1 <- arima(serie_temporal, order = c(12, 0, 1))

# Predecir usando el modelo ajustado
predicciones_modelo1 <- forecast(modelo1, h = length(serie_temporal))$mean
ajustes_modelo1 <- fitted(modelo1)

# Crear dataframe para ggplot2
df <- data.frame(
  Fecha = index(serie_temporal),
  Tmax = coredata(serie_temporal),
  Ajustes_modelo1 = ajustes_modelo1,
  Predicciones_modelo1 = predicciones_modelo1
)

intervalos_confianza <- forecast(modelo1, h = length(serie_temporal), level = 0.95)
# Calcular intervalos de confianza
intervalos_df <- data.frame(
  Fecha = index(serie_temporal),
  Lower = forecast(modelo1, h = length(serie_temporal))$lower,
  Upper = forecast(modelo1, h = length(serie_temporal))$upper
)

# Verificar la estructura de intervalos_df
str(intervalos_df)

# Renombrar las columnas en intervalos_df
names(intervalos_df) <- c("Fecha", "Lower", "Upper")

# plot
grafico <- ggplot(data = df, aes(x = Fecha)) +
  geom_line(aes(y = Tmax, color = "Original"), linewidth = 0.70) +
  geom_line(aes(y = Ajustes_modelo1, color = "Ajuste Modelo 1"), linewidth = 0.70, linetype = "dashed") +
  geom_line(aes(y = Predicciones_modelo1, color = "Predicción Modelo 1"), linewidth = 0.70, linetype = "dashed") +
  geom_ribbon(data = intervalos_df, aes(ymin = Lower, ymax = Upper), fill = "lightblue", alpha = 0.5) +
  labs(x = "Fecha", y = "Temperatura Máxima (°C)", 
       title = "Serie Temporal de Temperaturas Máximas y Predicciones con Modelo ARIMA 1",
       color = "Variables") +
  scale_x_date(date_breaks = "48 month", date_labels = "%b %Y") +
  theme_bw() +
  theme(legend.position = "topright") +
  guides(
    color = guide_legend(order = 1),  # Orden de la leyenda de las líneas
    fill = guide_legend(order = 2)    # Orden de la leyenda de la cinta
  )




# Realizar el test ADF
adf_result <- adf.test(serie_temporal, alternative = "stationary")

# Extraer el p-value
p_value <- adf_result$p.value

# Imprimir el resultado en consola
if (p_value > 0.05) {
  print("La serie NO es estacionaria (p-value > 0.05)")
} else {
  print("La serie es estacionaria (p-value <= 0.05)")
}

par(mfrow=c(2,1), mar=c(4,4,4,1)+.1)

# Gráfico de autocorrelación
acf(serie_temporal, main="Autocorrelación de la serie temporal, frecuencia = 12", cex.main=1)


# Gráfico de autocorrelación parcial
pacf(serie_temporal, main="Autocorrelación parcial de la serie temporal, frecuencia = 12", cex.main=1)


# Dibujar el gráfico
print(grafico)

# Calcular las funciones de autocorrelación
acf_result <- acf(serie_temporal, plot = FALSE)
pacf_result <- pacf(serie_temporal, plot = FALSE)


# Calcular errores
residuos <- resid(modelo1)

# Calcular métricas
mse <- mean(residuos^2)
mae <- mean(abs(residuos))
mape <- mean(abs(residuos / serie_temporal)) * 100
r_cuadrado <- 1 - sum(residuos^2) / sum((serie_temporal - mean(serie_temporal))^2)

# Imprimir métricas
cat("Error Cuadrático Medio (MSE):", mse, "\n")
cat("Error Absoluto Medio (MAE):", mae, "\n")
cat("Error Porcentual Absoluto Medio (MAPE):", mape, "%\n")
cat("Coeficiente de Determinación (R^2):", r_cuadrado, "\n")


