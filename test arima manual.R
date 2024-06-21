library(ncdf4)
library(ggplot2)
library(grid)
library(xts)
library(tseries)
library(forecast)  # Para suavizado exponencial

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
acf(serie_temporal, main="Autocorrelación de la serie temporal", cex.main=1)


# Gráfico de autocorrelación parcial
pacf(serie_temporal, main="Autocorrelación parcial de la serie temporal", cex.main=1)

# Título para el gráfico de autocorrelación parcial
title(main = "Autocorrelación parcial", line = -2, cex.main = 1)

acf(ts(serie_temporal, frequency=12))
pacf(ts(serie_temporal, frequency=12))

# Calcular las funciones de autocorrelación
acf_result <- acf(serie_temporal, plot = FALSE)
pacf_result <- pacf(serie_temporal, plot = FALSE)

# Encontrar los órdenes AR y MA
ar_orders <- which(acf_result$acf[-1] > 1.96 / sqrt(length(serie_temporal)))
ma_orders <- which(pacf_result$acf[-1] > 1.96 / sqrt(length(serie_temporal)))

# Mostrar los órdenes encontrados
print(paste("Orden AR:", ar_orders))
print(paste("Orden MA:", ma_orders))

modelo1  = arima(serie_temporal, order=c(12,0,1))
modelo1
tsdiag(modelo1)
Box.test(residuals(modelo1), type = "Ljung-Box")

modelo = modelo <- auto.arima(serie_temporal)
Box.test(residuals(modelo), type = "Ljung-Box")

# plot de la serie temporal con ggplot2
grafico <- ggplot(data = fortify(serie_temporal), aes(x = Index, y = serie_temporal)) +
  geom_line() +
  labs(x = "Fecha", y = "Temperatura máxima (°C)", 
       title = "Serie temporal temperaturas máximas") +
  theme_bw()

# Dibujar el gráfico
print(grafico)
