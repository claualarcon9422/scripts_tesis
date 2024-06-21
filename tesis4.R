# Dividir los datos en conjuntos de entrenamiento y prueba (70% - 30%)
n <- nrow(df)
n_train <- round(0.7 * n)

# Conjunto de entrenamiento
datos_entrenamiento <- df[1:n_train, ]

# Conjunto de prueba
datos_prueba <- df[(n_train + 1):n, ]

# Entrenar el modelo de suavizado exponencial
modelo <- ets(datos_entrenamiento$Tmax)

# Generar predicciones para el conjunto de prueba
predicciones <- forecast(modelo, h = nrow(datos_prueba))

# Calcular métricas de precisión
metricas <- accuracy(predicciones, datos_prueba$Tmax)

# Calcular intervalos de confianza
intervalos_confianza <- forecast:::forecast.ets(modelo, h = nrow(datos_prueba), level = 0.95)

# Crear un data frame con los intervalos de confianza
intervalos_df <- data.frame(
  Fecha = datos_prueba$Fecha,
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
grafico <- ggplot() +
  geom_line(data = df, aes(x = Fecha, y = Tmax, color = "Original"), size = 0.70) +
  geom_line(data = data.frame(Fecha = datos_entrenamiento$Fecha, Prediccion = fitted(modelo)),
            aes(x = Fecha, y = Prediccion, color = "Ajuste"), size = 0.70) +
  geom_line(data = data.frame(Fecha = datos_prueba$Fecha, Prediccion = predicciones$mean),
            aes(x = Fecha, y = Prediccion, color = "Predicción"), size = 0.70) +
  geom_ribbon(data = intervalos_df, aes(x = Fecha, ymin = Lower, ymax = Upper, fill = "Intervalo de confianza"), alpha = 0.5)+
labs(x = "Fecha", y = "Precipitación (mm/día)", 
     title = "Serie Temporal de temperatura maxima y Predicciones con Suavizado Exponencial",
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
