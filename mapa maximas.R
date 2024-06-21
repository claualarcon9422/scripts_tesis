library(ncdf4)
library(ggplot2)
library(viridis)
library(lubridate)

# Abrir el archivo NetCDF
archivo <- "G:/scripts/CR2MET_tmax_v2.0_mon_1979_2019_005deg.nc"
nc_data <- nc_open(archivo)

# Obtener las dimensiones y la primera fecha disponible
lat <- ncvar_get(nc_data, "lat")
lon <- ncvar_get(nc_data, "lon")
time <- ncvar_get(nc_data, "time")[1]

# Convertir la dimensión de tiempo a una fecha legible
fecha_inicio <- as.Date("1978-12-15")
fecha_primera <- fecha_inicio %m+% months(time)

# Obtener la variable "tmax" para la primera fecha disponible
tmax <- ncvar_get(nc_data, "tmax", start = c(1, 1, 1), count = c(-1, -1, 1))

# Cerrar el archivo NetCDF
nc_close(nc_data)

# Crear una matriz con las coordenadas y los valores de temperatura para la primera fecha
lon_lat <- expand.grid(lon = lon, lat = lat)
data <- data.frame(lon_lat, tmax = as.vector(tmax))

# Generar el mapa con ggplot2
ggplot(data, aes(x = lon, y = lat, fill = tmax)) +
  geom_tile() +
  scale_fill_viridis(name = "Tmax (°C)", na.value = "transparent") +
  coord_equal() +
  labs(title = paste("Mapa de Temperaturas Máximas en", fecha_primera),
       x = "Longitud",
       y = "Latitud") +
  theme_minimal()
