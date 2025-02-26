diagnostic_plots <- function(m) {
  df <- m$data
  if (!is.null(m$formula.mods)) {
    df<- df[complete.cases(df),]
  }
  
  df$residuals      <- residuals(m, type="pearson")  # Resíduos padronizados
  df$ajusted_values <- fitted(m)  # Valores ajustados pelo modelo
  df$color <- ifelse(abs(df$residuals) > 2, "blue", "black") # outliers (for influential cases)
  
  # Plot 1 - Histograma dos resíduos
  p1 <- ggplot(df, aes(x = residuals)) +
    geom_histogram(bins = 30, fill = "lightblue", color = "black") +
    labs(title = "Residuals histogram", x = "Standardized residuals", y = "Frequency") +
    theme_minimal()+
    theme(
      plot.title = element_text(size = 10),       # Title size
      axis.title.x = element_text(size = 10),    # x-axis title size
      axis.title.y = element_text(size = 10)) # y-axis title size
  
  # Plot 2 - QQ Plot dos resíduos
  p2<-ggplot(df, aes(sample = residuals)) +
    stat_qq(aes(color = after_stat(ifelse(abs(sample) > 2, "blue", "black")))) + 
    stat_qq_line(color = "red") +
    scale_color_identity() + # Mantém as cores definidas manualmente
    labs(title = "Residuals QQ-Plot", x = "Theoretical Quantiles", y = "Observed Residuals") +
    theme_minimal()+
    theme(
      plot.title = element_text(size = 10),       # Title size
      axis.title.x = element_text(size = 10),    # x-axis title size
      axis.title.y = element_text(size = 10)) # y-axis title size
  
  # Plot 3 - Observado vs. Resíduos
  p3 <- ggplot(df, aes(x = yi, y = residuals, color = color)) +
    geom_point() + scale_color_identity() +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    labs(title = "Observed vs. Residuals", x = "Observed", y = "Standardized residuals") +
    theme_minimal()+
    theme(
      plot.title = element_text(size = 10),       # Title size
      axis.title.x = element_text(size = 10),    # x-axis title size
      axis.title.y = element_text(size = 10)) # y-axis title size
  
  # Plot 4 - Ajustado vs. Resíduos
  p4 <- ggplot(df, aes(x = ajusted_values, y = residuals, color = color)) +
    geom_point() + scale_color_identity() +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    labs(title = "Fitted vs. Residuals", x = "Fitted values", y = "Standardized residuals") +
    theme_minimal()+
    theme(
      plot.title = element_text(size = 10),       # Title size
      axis.title.x = element_text(size = 10),    # x-axis title size
      axis.title.y = element_text(size = 10)) # y-axis title size
  
  # Plot 5 - Resíduos por estudo (ocupa 2 linhas)
  p5 <- ggplot(data = df, aes(y = Study, x = residuals, color = color)) +
    geom_jitter(height = 0.2, size = 2) + scale_color_identity() + 
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    labs(y = "Study", x = "Standardized residuals", title = "Residuals for each Study") +
    theme_minimal()+
    theme(
      plot.title = element_text(size = 10),       # Title size
      axis.title.x = element_text(size = 10),    # x-axis title size
      axis.title.y = element_text(size = 10)) # y-axis title size
  
  # Organizando os gráficos
  layout <- "
             ABE
             CDE
             "
  
  figure <- p1 + p2 + p3 + p4 + p5 + plot_layout(design = layout)
  return(figure)
}
