library(ncdf4)
library(forecast)
library(xts)

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

# Ajustar un modelo SARIMA
modelo_sarima <- Arima(serie_temporal, order=c(0,0,1), seasonal=c(1,0,1))
print(modelo_sarima)

# Diagnosticar el modelo
tsdiag(modelo_sarima)

# Prueba de Ljung-Box para los residuos
Box.test(residuals(modelo_sarima), type = "Ljung-Box")

# Plot de la serie temporal y pronósticos
plot(forecast(modelo_sarima))
