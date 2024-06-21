library(ncdf4)
library(dplyr)  
library(ggplot2)
library(grid)
library(forecast)
library(xts)
library(stats)
# Ruta de archivos NetCDF
#archivo <- "CR2MET_tmax_v2.0_mon_1979_2019_005deg.nc"
archivo <- "CR2MET_pr_v2.0_mon_1979_2019_005deg.nc"

# Abrir el archivo NetCDF
nc_data <- nc_open(archivo)

# Obtener las dimensiones y variables
lon <- ncvar_get(nc_data, "lon")
lat <- ncvar_get(nc_data, "lat")
#tmax <- ncvar_get(nc_data, "tmax")
pr <- ncvar_get(nc_data, "pr")

# Calcular el índice del valor central para lon y lat
indice_central_lon <- length(lon) %/% 2
indice_central_lat <- length(lat) %/% 2

# Obtener las coordenadas de la coordenada central
lon_central <- lon[indice_central_lon]
lat_central <- lat[indice_central_lat]

# Obtener las temperaturas máximas para la coordenada central
tmax_central <- tmax[indice_central_lon, indice_central_lat, ]

# Crear una serie temporal
fecha_inicio <- as.Date("1978-12-15")
fechas <- seq(fecha_inicio, by = "months", length.out = length(tmax_central))
serie_temporal <- xts(tmax_central, order.by = fechas)

# Cerrar el archivo NetCDF
nc_close(nc_data)

# Definir el data frame para ggplot
df <- data.frame(fecha = index(serie_temporal), valor = coredata(serie_temporal))

# Graficar la serie temporal de temperaturas
ggplot(df, aes(x = fecha, y = valor)) +
  geom_line(color = "blue") +
  labs(x = "Fecha", y = "Temperatura Máxima", title = "Serie Temporal de Temperatura Máxima") +
  theme_minimal()

# Graficar la serie temporal de precipitaciones
ggplot(df, aes(x = fecha, y = valor)) +
  geom_line(color = "blue") +
  labs(x = "Fecha", y = "Precipitación Máxima", title = "Serie Temporal de Precipitación Máxima") +
  theme_minimal()

par(mfrow=c(2,1), mar=c(4,4,4,1)+.1)

# Gráfico de autocorrelación
acf(serie_temporal, main="Autocorrelación de la serie temporal", cex.main=1)


# Gráfico de autocorrelación parcial
pacf(serie_temporal, main="Autocorrelación parcial de la serie temporal", cex.main=1)

# Aplicar el test de Ljung-Box para varios rezagos
max_lag <- 50

ljung_box_test <- Box.test(serie_temporal, lag = max_lag, type = "Ljung-Box")

# Mostrar el resultado del test
print(ljung_box_test)

p_values <- sapply(1:max_lag, function(lag) Box.test(serie_temporal, lag = lag, type = "Ljung-Box")$p.value)

# Crear un dataframe con los resultados
resultados <- data.frame(lag = 1:max_lag, p_value = p_values)

# Graficar los p-valores
ggplot(resultados, aes(x = lag, y = p_value)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = "P-valores del test de Ljung-Box para diferentes rezagos",
       x = "Rezago",
       y = "P-valor") +
  theme_minimal()

# Crear un data frame con la serie temporal
df <- data.frame(Fecha = index(serie_temporal), Tmax = coredata(serie_temporal))

# Dividir los datos en conjuntos de entrenamiento y prueba (70% - 30%)
n <- nrow(df)
n_train <- round(0.7 * n)

# Conjunto de entrenamiento
datos_entrenamiento <- df[1:n_train, ]

# Conjunto de prueba
datos_prueba <- df[(n_train + 1):n, ]

# # Entrenar el modelo ARMA
 #modelo  = arima(datos_entrenamiento$Tmax, order = c(1,1,0))
 #titulo_grafico <- "Serie Temporal de temperatura maxima y Predicciones con ARMA"

#modelo  = auto.arima(datos_entrenamiento$Tmax)
#titulo_grafico <- "Serie Temporal de temperatura maxima y Predicciones con ARIMA"

 modelo <- Arima(datos_entrenamiento$Tmax, order=c(12,0,10), seasonal=c(12,0,1))
 titulo_grafico <- "Serie Temporal de temperatura maxima y Predicciones con SARIMA"


#Estructura interna del modelo
str(modelo)

#Gráfico de residuos estandarizados, ACF de residuos y valores p de test Ljung-Box
tsdiag(modelo)

# Generar predicciones para el conjunto de prueba
predicciones <- forecast(modelo, h = nrow(datos_prueba))

# Calcular métricas de precisión
metricas <- accuracy(predicciones, datos_prueba$Tmax)
# Calcular el R^2
r_cuadrado <- 1 - (sum((datos_prueba$Tmax - predicciones$mean)^2) / sum((datos_prueba$Tmax - mean(datos_prueba$Tmax))^2))

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

#print(modelo)

# plot
grafico <- ggplot() +
  geom_line(data = df, aes(x = Fecha, y = Tmax, color = "Original"), linewidth = 0.70) +
  geom_line(data = data.frame(Fecha = datos_entrenamiento$Fecha, Prediccion = fitted(modelo)),
            aes(x = Fecha, y = Prediccion, color = "Ajuste"), linewidth = 0.70) +
  geom_line(data = data.frame(Fecha = datos_prueba$Fecha, Prediccion = predicciones$mean),
            aes(x = Fecha, y = Prediccion, color = "Predicción"), linewidth = 0.70) +
  geom_ribbon(data = intervalos_df, aes(x = Fecha, ymin = Lower, ymax = Upper, fill = "Intervalo de confianza"), alpha = 0.5) +
  labs(x = "Fecha", y = "T° Máxima (°C/mes)", 
       title = titulo_grafico,
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




