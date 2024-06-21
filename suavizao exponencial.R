library(ncdf4)
library(dplyr)  
library(ggplot2)
library(grid)
library(xts)
library(tseries)
library(forecast)  # Para suavizado exponencial
library(viridis)
library(lubridate)

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

# Dividir los datos en conjuntos de entrenamiento y prueba (70% - 30%)
n <- length(serie_temporal)
n_train <- round(0.7 * n)

# Conjunto de entrenamiento
datos_entrenamiento <- serie_temporal[1:n_train]

# Conjunto de prueba
datos_prueba <- serie_temporal[(n_train + 1):n]

# Entrenar el modelo de suavizado exponencial
modelo <- ets(datos_entrenamiento)

# Generar predicciones para el conjunto de prueba
predicciones <- forecast(modelo, h = nrow(datos_prueba))

# Calcular métricas de precisión
metricas <- accuracy(predicciones, datos_prueba)


# Calcular intervalos de confianza
intervalos_confianza <- forecast:::forecast.ets(modelo, h = nrow(datos_prueba), level = 0.95)

# Crear un data frame con los intervalos de confianza
intervalos_df <- data.frame(
  Fecha = index(datos_prueba),
  Lower = intervalos_confianza$lower,
  Upper = intervalos_confianza$upper
)
# Renombrar las columnas
colnames(intervalos_df) <- c("Fecha", "Lower", "Upper")

parametros <- coef(modelo)
alfa <- parametros["alpha"]
beta <- parametros["beta"]
gamma <- parametros["gamma"]

texto_metricas <- paste("RMSE:", round(metricas[1, "RMSE"], 4),
                        "\nMAE:", round(metricas[1, "MAE"], 4),
                        "\nMAPE:", round(metricas[1, "MAPE"], 4),
                        "\nMASE:", round(metricas[1, "MASE"], 4))

texto_alfa <- paste("\n\u03B1:", round(alfa, 4))


print(modelo)


# plot
grafico <- ggplot(data = as.data.frame(serie_temporal), aes(x = index(serie_temporal))) +
  geom_line(aes(y = serie_temporal[,1], color = "Original"), linewidth = 0.70) +
  geom_line(data = data.frame(Fecha = index(datos_entrenamiento), Prediccion = fitted(modelo)),
            aes(x = Fecha, y = Prediccion, color = "Ajuste Modelo 1"), linewidth = 0.70, linetype = "dashed") +
  geom_line(data = data.frame(Fecha = index(datos_prueba), Prediccion = predicciones$mean),
            aes(x = Fecha, y = Prediccion, color = "Predicción Modelo 1"), linewidth = 0.70, linetype = "dashed") +
  geom_ribbon(data = intervalos_df, aes(x = Fecha, ymin = Lower, ymax = Upper, fill = "Intervalo de confianza"), alpha = 0.5)+
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

# Añadir un cuadradito manual para el intervalo de confianza en la leyenda
grafico <- grafico + 
  scale_fill_manual(name = "", values = "lightblue",
                    labels = c("Intervalo de confianza" = "Intervalo de confianza"))

# Dibujar el gráfico
grid.newpage()
grid.draw(grafico)
xoff <-0.825
yoff <- -0.2

grid.text(label = "Métricas", x = unit(xoff, "npc"), y = unit(0.74 + yoff, "npc"),
          just = "left", gp = gpar(fontsize = 12, col = "black"))
grid.text(label = texto_metricas, x = unit(xoff+0.01, "npc"), y = unit(0.65 + yoff, "npc"),
          just = "left", gp = gpar(col = "black", fontsize = 10))

grid.text(label = "Factor de suavizado", x = unit(xoff, "npc"), y = unit(0.5 + yoff, "npc"),
          just = "left", gp = gpar(fontsize = 12, col = "black"))
grid.text(label = texto_alfa, x = unit(xoff+0.01, "npc"), y = unit(0.48 + yoff, "npc"),
          just = "left", gp = gpar(col = "black", fontsize = 10))

