library(ncdf4)
library(dplyr)  
library(ggplot2)
library(grid)
library(forecast)  # Para suavizado exponencial

# Función para extraer datos de precipitación para una coordenada específica de un archivo NetCDF
extraer_precipitacion <- function(archivo, lat, lon) {
  nc_data <- nc_open(archivo)
  precipitacion <- ncvar_get(nc_data, "pr")
  lon_index <- which.min(abs(ncvar_get(nc_data, "lon") - lon))
  lat_index <- which.min(abs(ncvar_get(nc_data, "lat") - lat))
  nc_close(nc_data)
  return(precipitacion[lon_index, lat_index, , drop = TRUE])
}

# Carpeta que contiene los archivos NetCDF
carpeta <- "D:/Descargas/CR2MET_pr_v2.5_day_1960-2021_005deg"

# Coordenadas específicas
lat_especifica <- -55.125
lon_especifica <- -73.025

# Lista para almacenar los datos de precipitación de cada archivo
lista_precipitacion <- list()

# Leer cada archivo en la carpeta y extraer los datos de precipitación
archivos <- list.files(path = carpeta, pattern = "\\.nc$", full.names = TRUE)
for (archivo in archivos) {
  precipitacion_mes <- extraer_precipitacion(archivo, lat_especifica, lon_especifica)
  lista_precipitacion[[length(lista_precipitacion) + 1]] <- precipitacion_mes
}

# Combinar los datos de precipitación de todos los meses
precipitacion_total <- do.call(c, lista_precipitacion)

# Crear un data frame para almacenar las fechas y la precipitación
fecha_inicio <- as.Date("2021-01-01")
fechas <- seq(fecha_inicio, by = "day", length.out = length(precipitacion_total))
datos_precipitacion <- data.frame(Fecha = fechas, Precipitacion = precipitacion_total)

# Dividir los datos en conjuntos de entrenamiento y prueba (70% - 30%)
n <- nrow(datos_precipitacion)
n_train <- round(0.7 * n)

# Conjunto de entrenamiento
datos_entrenamiento <- datos_precipitacion[1:n_train, ]

# Conjunto de prueba
datos_prueba <- datos_precipitacion[(n_train + 1):n, ]

# Entrenar el modelo ARIMA
modelo_arima <- auto.arima(datos_entrenamiento$Precipitacion)

# Generar predicciones para el conjunto de prueba
predicciones_arima <- forecast(modelo_arima, h = nrow(datos_prueba))

# Calcular métricas de precisión
metricas_arima <- accuracy(predicciones_arima, datos_prueba$Precipitacion)

# Calcular intervalos de confianza
intervalos_confianza_arima <- forecast:::forecast.Arima(modelo_arima, h = nrow(datos_prueba), level = 0.95)

# Crear un data frame con los intervalos de confianza
intervalos_df_arima <- data.frame(
  Fecha = datos_prueba$Fecha,
  Lower = intervalos_confianza_arima$lower,
  Upper = intervalos_confianza_arima$upper
)
# Renombrar las columnas
colnames(intervalos_df_arima) <- c("Fecha", "Lower", "Upper")

texto_metricas_arima <- paste("RMSE:", round(metricas_arima[1, "RMSE"], 4),
                              "\nMAE:", round(metricas_arima[1, "MAE"], 4),
                              "\nMAPE:", round(metricas_arima[1, "MAPE"], 4),
                              "\nMASE:", round(metricas_arima[1, "MASE"], 4))

# plot
grafico <- ggplot() +
  geom_line(data = datos_precipitacion, aes(x = Fecha, y = Precipitacion, color = "Original"), size = 0.70) +
  geom_line(data = data.frame(Fecha = datos_entrenamiento$Fecha, Prediccion = fitted(modelo_arima)),
            aes(x = Fecha, y = Prediccion, color = "Ajuste ARIMA"), size = 0.70) +
  geom_line(data = data.frame(Fecha = datos_prueba$Fecha, Prediccion = predicciones_arima$mean),
            aes(x = Fecha, y = Prediccion, color = "Predicción ARIMA"), size = 0.70) +
  geom_ribbon(data = intervalos_df_arima, aes(x = Fecha, ymin = Lower, ymax = Upper, fill = "Intervalo de confianza ARIMA"), alpha = 0.5) +
  labs(x = "Fecha", y = "Precipitación (mm/día)", 
       title = "Serie Temporal de Precipitación y Predicciones con Modelo ARIMA",
       color = "Variables") +
  scale_x_date(date_breaks = "2 month", date_labels = "%b %Y") +
  theme_bw() +
  coord_cartesian(clip = "off") +
  theme(legend.position = c(1.01, 1),
        legend.justification = c(0,1))+ 
  theme(plot.margin = margin(5.5, 130, 5.5, 5.5)) +
  guides(
    color = guide_legend(order = 1),  # Orden de la leyenda de las líneas
    fill = guide_legend(order = 2)    # Orden de la leyenda de la cinta
  )

# Añadir un cuadradito manual para el intervalo de confianza en la leyenda
grafico <- grafico + 
  scale_fill_manual(name = "", values = "lightblue",
                    labels = c("Intervalo de confianza ARIMA" = "Intervalo de confianza"))

# Dibujar el gráfico
grid.newpage()
grid.draw(grafico)
xoff <-0.825
yoff <- -0.2

grid.text(label = "Métricas ARIMA", x = unit(xoff, "npc"), y = unit(0.74 + yoff, "npc"),
          just = "left", gp = gpar(fontsize = 12, col = "black"))
grid.text(label = texto_metricas_arima, x = unit(xoff+0.01, "npc"), y = unit(0.65 + yoff, "npc"),
          just = "left", gp = gpar(col = "black", fontsize = 10))
