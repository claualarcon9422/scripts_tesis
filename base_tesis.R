library(ncdf4)
library(dplyr)  
library(ggplot2)
library(grid)
library(forecast)
library(xts)
library(stats)
library(gridExtra)
library(tseries)

# Ruta de archivos NetCDF
archivo_tmax <- "CR2MET_tmax_v2.0_mon_1979_2019_005deg.nc"
archivo_pr <- "CR2MET_pr_v2.0_mon_1979_2019_005deg.nc"

# Función para cargar y procesar los datos
cargar_datos <- function(archivo, variable) {
  # Abrir el archivo NetCDF
  nc_data <- nc_open(archivo)
  
  # Obtener las dimensiones y variables
  lon <- ncvar_get(nc_data, "lon")
  lat <- ncvar_get(nc_data, "lat")
  data <- ncvar_get(nc_data, variable)
  
  # Calcular el índice del valor central para lon y lat
  indice_central_lon <- length(lon) %/% 2
  indice_central_lat <- length(lat) %/% 2
  
  # Obtener los datos para la coordenada central
  datos_centrales <- data[indice_central_lon, indice_central_lat, ]
  
  # Crear una serie temporal
  fecha_inicio <- as.Date("1978-12-15")
  fechas <- seq(fecha_inicio, by = "months", length.out = length(datos_centrales))
  serie_temporal <- xts(datos_centrales, order.by = fechas)
  
  # Cerrar el archivo NetCDF
  nc_close(nc_data)
  
  return(serie_temporal)
}

# Cargar datos de temperatura máxima
serie_temporal_tmax <- cargar_datos(archivo_tmax, "tmax")

# Cargar datos de precipitación
serie_temporal_pr <- cargar_datos(archivo_pr, "pr")

# Función para graficar la serie temporal
graficar_serie_temporal <- function(serie, variable, titulo) {
  # Definir el data frame para ggplot
  df <- data.frame(fecha = index(serie), valor = coredata(serie))
  
  # Graficar la serie temporal
  ggplot(df, aes(x = fecha, y = valor)) +
    geom_line(color = "blue") +
    labs(x = "Fecha", y = variable, title = titulo) +
    theme_minimal()
}


# Función para aplicar el test de Ljung-Box y graficar
aplicar_test_ljung_box <- function(serie, max_lag, titulo) {
  p_values <- sapply(1:max_lag, function(lag) Box.test(serie, lag = lag, type = "Ljung-Box")$p.value)
  
  # Crear un dataframe con los resultados
  resultados <- data.frame(lag = 1:max_lag, p_value = p_values)
  
  # Graficar los p-valores
  gg <- ggplot(resultados, aes(x = lag, y = p_value)) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    labs(title = paste("P-valores del test de Ljung-Box para", titulo),
         x = "Rezago",
         y = "P-valor") +
    theme_minimal()
  
  return(gg)
}

# Graficar serie temporal de temperatura máxima
graficar_serie_temporal(serie_temporal_tmax, "Temperatura Máxima", "Serie Temporal de Temperaturas Máximas")

# Graficar serie temporal de precipitación
graficar_serie_temporal(serie_temporal_pr, "Precipitación", "Serie Temporal de Precipitaciones Máximas")


# Crear gráficos de Ljung-Box y ACF/PACF
par(mfrow=c(4, 1), mar=c(4,4,2,1))  # Ajustar mar y mfrow

# Gráficos de Ljung-Box
ljung_box_tmax <- aplicar_test_ljung_box(serie_temporal_tmax, max_lag = 50, titulo = "Temperatura Máxima")
ljung_box_pr <- aplicar_test_ljung_box(serie_temporal_pr, max_lag = 50, titulo = "Precipitación")

# Mostrar gráficos de Ljung-Box con grid.arrange
grid.arrange(ljung_box_tmax, ljung_box_pr, ncol = 1)


max_lag <- 50
ljung_box_test_tmax <- Box.test(serie_temporal_tmax, lag = max_lag, type = "Ljung-Box")
print(ljung_box_test_tmax)
ljung_box_test_pr <- Box.test(serie_temporal_pr, lag = max_lag, type = "Ljung-Box")
print(ljung_box_test_pr)


par(mfrow=c(2,1), mar=c(4,4,4,1)+.1)
# Gráfico de autocorrelación
acf(serie_temporal_tmax, main="Autocorrelación de la Serie de Temperaturas Máximas", cex.main=1)
# Gráfico de autocorrelación parcial
pacf(serie_temporal_tmax, main="Autocorrelación Parcial de la Serie de Temperaturas Máximas", cex.main=1)

# Gráfico de autocorrelación
acf(serie_temporal_pr, main="Autocorrelación de la Serie de Precipitaciones Máximas", cex.main=1)
# Gráfico de autocorrelación parcial
pacf(serie_temporal_pr, main="Autocorrelación Parcial de la Serie de Precipitaciones Máximas", cex.main=1)

# Función para normalizar una serie temporal entre 0 y 1
normalizar_serie <- function(serie) {
  min_val <- min(serie)
  max_val <- max(serie)
  serie_normalizada <- (serie - min_val) / (max_val - min_val)
  return(serie_normalizada)
}

#Test Dickey-Fuller Aumentada
adf.test(serie_temporal_tmax)
adf.test(serie_temporal_pr)

# Normalizar serie de temperatura máxima
serie_temporal_tmax <- normalizar_serie(serie_temporal_tmax)
# Normalizar serie de precipitación
serie_temporal_pr <- normalizar_serie(serie_temporal_pr)




# Función para ajustar y graficar modelo ARIMA
ajustar_y_graficar_arima <- function(serie, titulo) {
  # Crear un data frame con la serie temporal
  df <- data.frame(Fecha = index(serie), Valor = coredata(serie))
  
  # Dividir los datos en conjuntos de entrenamiento y prueba (70% - 30%)
  n <- nrow(df)
  n_train <- round(0.7 * n)
  
  # Conjunto de entrenamiento
  datos_entrenamiento <- df[1:n_train, ]
  
  # Conjunto de prueba
  datos_prueba <- df[(n_train + 1):n, ]
  
  # Entrenar el modelo ARIMA
  modelo <- auto.arima(datos_entrenamiento$Valor)
  str(modelo)
  
  # Generar predicciones para el conjunto de prueba
  predicciones <- forecast(modelo, h = nrow(datos_prueba))
  
  # Calcular métricas de precisión
  metricas <- accuracy(predicciones, datos_prueba$Valor)
  
  # Calcular el R^2
  r_cuadrado <- 1 - (sum((datos_prueba$Valor - predicciones$mean)^2) / sum((datos_prueba$Valor - mean(datos_prueba$Valor))^2))
  
  # Calcular intervalos de confianza
  intervalos_confianza <- forecast:::forecast.Arima(modelo, h = nrow(datos_prueba), level = 0.95)
  
  # Crear un data frame con los intervalos de confianza
  intervalos_df <- data.frame(
    Fecha = datos_prueba$Fecha,
    Lower = intervalos_confianza$lower,
    Upper = intervalos_confianza$upper
  )
  # Renombrar las columnas
  colnames(intervalos_df) <- c("Fecha", "Lower", "Upper")
  
  texto_metricas <- paste("RMSE:", round(metricas[1, "RMSE"], 4),
                          "\nMAE:", round(metricas[1, "MAE"], 4),
                          "\nMAPE:", round(metricas[1, "MAPE"], 4),
                          "\nMASE:", round(metricas[1, "MASE"], 4),
                          "\nR²:", round(r_cuadrado, 4))
  
  # Gráfico de diagnóstico del modelo ARIMA
  tsdiag(modelo)
  
  # Plot de la serie temporal y predicciones
  grafico <- ggplot() +
    geom_line(data = df, aes(x = Fecha, y = Valor, color = "Original"), linewidth = 0.70) +
    geom_line(data = data.frame(Fecha = datos_entrenamiento$Fecha, Prediccion = fitted(modelo)),
              aes(x = Fecha, y = Prediccion, color = "Ajuste"), linewidth = 0.70) +
    geom_line(data = data.frame(Fecha = datos_prueba$Fecha, Prediccion = predicciones$mean),
              aes(x = Fecha, y = Prediccion, color = "Predicción"), linewidth = 0.70) +
    geom_ribbon(data = intervalos_df, aes(x = Fecha, ymin = Lower, ymax = Upper, fill = "Intervalo de confianza"), alpha = 0.5) +
    labs(x = "Fecha", y = "Valor", 
         title = titulo,
         color = "Variables") +
    scale_x_date(date_breaks = "72 month", date_labels = "%b %Y") +
    theme_bw() +
    coord_cartesian(clip = "off") +
    theme(legend.position = c(1.01, 1),
          legend.justification = c(0,1))+ 
    theme(plot.margin = margin(5.5, 130, 5.5, 5.5)) +
    guides(
      color = guide_legend(order = 1),  # Orden de la leyenda de las líneas
      fill = guide_legend(order = 2)    # Orden de la leyenda de la cinta
    )+
    scale_color_manual(values = c("Original" = "blue", "Ajuste" = rgb(0,128,0,maxColorValue = 255), "Predicción"  = rgb(255,0,0,maxColorValue = 255)))
  
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
  grid.text(label = texto_metricas, x = unit(xoff+0.01, "npc"), y = unit(0.6 + yoff, "npc"),
            just = "left", gp = gpar(col = "black", fontsize = 10))
}


# Ajustar y graficar modelo ARIMA para temperatura máxima normalizada
ajustar_y_graficar_arima(serie_temporal_tmax, "Serie Temporal Normalizada de Temperatura Máxima y Predicciones con ARIMA")

# Ajustar y graficar modelo ARIMA para precipitación normalizada
ajustar_y_graficar_arima(serie_temporal_pr, "Serie Temporal Normalizada de Precipitación y Predicciones con ARIMA")

