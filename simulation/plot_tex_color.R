# Install required packages if not already installed
# install.packages(c("ggplot2", "dplyr", "tibble", "patchwork", "tikzDevice"))

library(ggplot2)
library(dplyr)
library(tibble)
library(patchwork)
library(tikzDevice)

options(tikzLatexPackages = c(getOption("tikzLatexPackages"), 
                              "\\usepackage{amsmath}", 
                              "\\usepackage{amssymb}"))

# ==========================================
# 1. DATA TRANSCRIPTION
# ==========================================

# =====================================================================
# Scenario I: True Model Recovery & Expected Size
# =====================================================================
df_scen1 <- tribble(
  ~Scenario,    ~K,  ~Method,      ~m,  ~Prob,  ~Size,
  
  # K = 6
  "Scenario I",  6,  "LASSO",      10, 0.9990, 3.0845, #
  "Scenario I",  6,  "LASSO+FTBS", 10, 0.9440, 1.1320, #
  "Scenario I",  6,  "FTBS",       10, 0.8340, 0.8990, #
  "Scenario I",  6,  "LASSO",      20, 1.0000, 3.0310, #
  "Scenario I",  6,  "LASSO+FTBS", 20, 1.0000, 1.1680, #
  "Scenario I",  6,  "FTBS",       20, 1.0000, 1.0080, #
  "Scenario I",  6,  "LASSO",      30, 1.0000, 3.0205, #
  "Scenario I",  6,  "LASSO+FTBS", 30, 1.0000, 1.1720, #
  "Scenario I",  6,  "FTBS",       30, 1.0000, 1.0070, #
  
  # K = 10
  "Scenario I", 10,  "LASSO",      10, 0.9950, 4.0625, #
  "Scenario I", 10,  "LASSO+FTBS", 10, 0.9660, 1.9260, #
  "Scenario I", 10,  "FTBS",       10, 0.5510, 0.8215, #
  "Scenario I", 10,  "LASSO",      20, 1.0000, 4.0715, #
  "Scenario I", 10,  "LASSO+FTBS", 20, 1.0000, 1.9255, #
  "Scenario I", 10,  "FTBS",       20, 0.9910, 1.0060, #
  "Scenario I", 10,  "LASSO",      30, 1.0000, 4.0660, #
  "Scenario I", 10,  "LASSO+FTBS", 30, 1.0000, 1.8975, #
  "Scenario I", 10,  "FTBS",       30, 1.0000, 1.0065, #
  
  # K = 20
  "Scenario I", 20,  "LASSO",      10, 0.9790, 4.1525, #
  "Scenario I", 20,  "LASSO+FTBS", 10, 0.9490, 2.9035, #
  "Scenario I", 20,  "FTBS",       10, 0.2480, 0.6880, #
  "Scenario I", 20,  "LASSO",      20, 1.0000, 4.2145, #
  "Scenario I", 20,  "LASSO+FTBS", 20, 0.9990, 2.9740, #
  "Scenario I", 20,  "FTBS",       20, 0.8520, 0.9365, #
  "Scenario I", 20,  "LASSO",      30, 1.0000, 4.1970, #
  "Scenario I", 20,  "LASSO+FTBS", 30, 1.0000, 2.9710, #
  "Scenario I", 20,  "FTBS",       30, 0.9920, 1.0015, #
  
  # K = 50
  "Scenario I", 50,  "LASSO",      10, 0.8250, 2.6425, #
  "Scenario I", 50,  "LASSO+FTBS", 10, 0.8050, 2.4320, #
  "Scenario I", 50,  "FTBS",       10, 0.1250, 0.7240, #
  "Scenario I", 50,  "LASSO",      20, 1.0000, 2.9715, #
  "Scenario I", 50,  "LASSO+FTBS", 20, 1.0000, 2.7245, #
  "Scenario I", 50,  "FTBS",       20, 0.5270, 0.8010, #
  "Scenario I", 50,  "LASSO",      30, 1.0000, 3.2165, #
  "Scenario I", 50,  "LASSO+FTBS", 30, 1.0000, 2.9045, #
  "Scenario I", 50,  "FTBS",       30, 0.9070, 0.9575  #
)

# =====================================================================
# Scenario II: True Model Recovery & Expected Size
# =====================================================================
df_scen2 <- tribble(
  ~Scenario,     ~K,  ~Method,      ~m,  ~Prob,  ~Size,
  
  # K = 6
  "Scenario II",  6,  "LASSO",      10, 0.9980, 2.5967, #
  "Scenario II",  6,  "LASSO+FTBS", 10, 0.7810, 0.9013, #
  "Scenario II",  6,  "FTBS",       10, 0.7260, 0.8000, #
  "Scenario II",  6,  "LASSO",      20, 1.0000, 2.5043, #
  "Scenario II",  6,  "LASSO+FTBS", 20, 0.9540, 1.0730, #
  "Scenario II",  6,  "FTBS",       20, 0.9890, 1.0117, #
  "Scenario II",  6,  "LASSO",      30, 1.0000, 2.4673, #
  "Scenario II",  6,  "LASSO+FTBS", 30, 0.9770, 1.0963, #
  "Scenario II",  6,  "FTBS",       30, 1.0000, 1.0103, #
  
  # K = 10
  "Scenario II", 10,  "LASSO",      10, 0.9950, 3.8753, #
  "Scenario II", 10,  "LASSO+FTBS", 10, 0.9310, 1.6723, #
  "Scenario II", 10,  "FTBS",       10, 0.4000, 0.8593, #
  "Scenario II", 10,  "LASSO",      20, 1.0000, 3.8600, #
  "Scenario II", 10,  "LASSO+FTBS", 20, 0.9960, 1.7393, #
  "Scenario II", 10,  "FTBS",       20, 0.9760, 1.0233, #
  "Scenario II", 10,  "LASSO",      30, 1.0000, 3.8510, #
  "Scenario II", 10,  "LASSO+FTBS", 30, 1.0000, 1.6540, #
  "Scenario II", 10,  "FTBS",       30, 0.9980, 1.0123, #
  
  # K = 20
  "Scenario II", 20,  "LASSO",      10, 0.9730, 4.1723, #
  "Scenario II", 20,  "LASSO+FTBS", 10, 0.9120, 2.7443, #
  "Scenario II", 20,  "FTBS",       10, 0.1590, 0.7440, #
  "Scenario II", 20,  "LASSO",      20, 1.0000, 4.2570, #
  "Scenario II", 20,  "LASSO+FTBS", 20, 0.9990, 2.8027, #
  "Scenario II", 20,  "FTBS",       20, 0.7860, 0.9770, #
  "Scenario II", 20,  "LASSO",      30, 1.0000, 4.1407, #
  "Scenario II", 20,  "LASSO+FTBS", 30, 1.0000, 2.7580, #
  "Scenario II", 20,  "FTBS",       30, 0.9880, 1.0107, #
  
  # K = 50
  "Scenario II", 50,  "LASSO",      10, 0.7850, 2.8023, #
  "Scenario II", 50,  "LASSO+FTBS", 10, 0.7550, 2.4890, #
  "Scenario II", 50,  "FTBS",       10, 0.0200, 0.6230, #
  "Scenario II", 50,  "LASSO",      20, 0.9990, 3.2817, #
  "Scenario II", 50,  "LASSO+FTBS", 20, 0.9950, 2.8897, #
  "Scenario II", 50,  "FTBS",       20, 0.3820, 0.8280, #
  "Scenario II", 50,  "LASSO",      30, 1.0000, 3.2467, #
  "Scenario II", 50,  "LASSO+FTBS", 30, 1.0000, 2.8500, #
  "Scenario II", 50,  "FTBS",       30, 0.8800, 0.9870  #
)

# =====================================================================
# Scenario III: True Model Recovery & Expected Size
# =====================================================================
df_scen3 <- tribble(
  ~Scenario,      ~K,  ~Method,      ~m,  ~Prob,  ~Size,
  
  # K = 6
  "Scenario III",  6,  "LASSO",      10, 0.9720, 3.6300, #
  "Scenario III",  6,  "LASSO+FTBS", 10, 0.8030, 0.9547, #
  "Scenario III",  6,  "FTBS",       10, 0.6740, 0.7363, #
  "Scenario III",  6,  "LASSO",      20, 1.0000, 3.8193, #
  "Scenario III",  6,  "LASSO+FTBS", 20, 0.9930, 1.0953, #
  "Scenario III",  6,  "FTBS",       20, 0.9940, 0.9970, #
  "Scenario III",  6,  "LASSO",      30, 1.0000, 3.8090, #
  "Scenario III",  6,  "LASSO+FTBS", 30, 1.0000, 1.1273, #
  "Scenario III",  6,  "FTBS",       30, 1.0000, 1.0010, #
  
  # K = 10
  "Scenario III", 10,  "LASSO",      10, 0.6950, 4.8283, #
  "Scenario III", 10,  "LASSO+FTBS", 10, 0.6290, 1.6247, #
  "Scenario III", 10,  "FTBS",       10, 0.4670, 0.8060, #
  "Scenario III", 10,  "LASSO",      20, 0.9850, 5.6343, #
  "Scenario III", 10,  "LASSO+FTBS", 20, 0.9750, 1.8147, #
  "Scenario III", 10,  "FTBS",       20, 0.9810, 1.0007, #
  "Scenario III", 10,  "LASSO",      30, 1.0000, 5.8233, #
  "Scenario III", 10,  "LASSO+FTBS", 30, 1.0000, 1.7760, #
  "Scenario III", 10,  "FTBS",       30, 1.0000, 1.0047, #
  
  # K = 20
  "Scenario III", 20,  "LASSO",      10, 0.1560, 4.1067, #
  "Scenario III", 20,  "LASSO+FTBS", 10, 0.1450, 2.4987, #
  "Scenario III", 20,  "FTBS",       10, 0.1830, 0.6500, #
  "Scenario III", 20,  "LASSO",      20, 0.8300, 6.0707, #
  "Scenario III", 20,  "LASSO+FTBS", 20, 0.7970, 3.0643, #
  "Scenario III", 20,  "FTBS",       20, 0.8160, 0.9497, #
  "Scenario III", 20,  "LASSO",      30, 0.9890, 6.6463, #
  "Scenario III", 20,  "LASSO+FTBS", 30, 0.9820, 3.2170, #
  "Scenario III", 20,  "FTBS",       30, 0.9900, 1.0030, #
  
  # K = 50
  "Scenario III", 50,  "LASSO",      10, 0.0080, 2.1983, #
  "Scenario III", 50,  "LASSO+FTBS", 10, 0.0050, 1.8577, #
  "Scenario III", 50,  "FTBS",       10, 0.0630, 0.5307, #
  "Scenario III", 50,  "LASSO",      20, 0.3420, 4.3337, #
  "Scenario III", 50,  "LASSO+FTBS", 20, 0.3120, 3.0937, #
  "Scenario III", 50,  "FTBS",       20, 0.4530, 0.8213, #
  "Scenario III", 50,  "LASSO",      30, 0.8310, 5.4850, #
  "Scenario III", 50,  "LASSO+FTBS", 30, 0.7820, 3.3223, #
  "Scenario III", 50,  "FTBS",       30, 0.8960, 0.9660  #
)

# ==========================================
# 2. STRICT WIDTH-CONSTRAINED THEME
# ==========================================

custom_lines <- c("LASSO" = "solid", "LASSO+FTBS" = "dashed", "FTBS" = "dotdash")
custom_shapes <- c("LASSO" = 16, "LASSO+FTBS" = 17, "FTBS" = 15) 
custom_colors <- c("LASSO" = "#E41A1C", "LASSO+FTBS" = "#377EB8", "FTBS" = "#4DAF4A") # RColorBrewer Set1 Colors

# Reduced base_size to 8 to fit narrow columns
academic_theme <- theme_bw(base_size = 8) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 8, margin = margin(2,0,2,0)),
    legend.position = "bottom",
    legend.title = element_blank(), 
    # Severely shortened legend lines to save horizontal space
    legend.key.width = unit(1.2, "cm"),
    legend.text = element_text(size = 8),
    legend.margin = margin(t = -5),
    # Tighten gap between columns so it fits in \textwidth
    panel.spacing.x = unit(0.2, "cm"), 
    panel.spacing.y = unit(0.2, "cm"),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 7)
  )

# ==========================================
# 3. PLOT A: SCENARIOS I & II (Combined)
# ==========================================

df_12 <- bind_rows(df_scen1, df_scen2) %>%
  mutate(
    Method = factor(Method, levels = c("LASSO", "LASSO+FTBS", "FTBS")),
    K_label = factor(K, levels = c(6, 10, 20, 50), labels = paste("K =", c(6, 10, 20, 50))),
    Scenario = factor(Scenario, levels = c("Scenario I", "Scenario II"))
  )

p_prob_12 <- ggplot(df_12, aes(x = m, y = Prob, shape = Method, linetype = Method, color = Method, group = Method)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  facet_grid(Scenario ~ K_label) +
  scale_linetype_manual(name = NULL, values = custom_lines) +
  scale_shape_manual(name = NULL, values = custom_shapes) +
  scale_color_manual(name = NULL, values = custom_colors) +
  scale_x_continuous(breaks = c(10, 20, 30)) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
  labs(
    x = "$m$",
    y = "$\\widehat{\\mathbb{P}}(\\mathcal{S} \\subset \\widehat{\\mathcal{S}})$"
  ) +
  academic_theme

p_size_12 <- ggplot(df_12, aes(x = m, y = Size, shape = Method, linetype = Method, color = Method, group = Method)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  facet_grid(Scenario ~ K_label) +
  scale_linetype_manual(name = NULL, values = custom_lines) +
  scale_shape_manual(name = NULL, values = custom_shapes) +
  scale_color_manual(name = NULL, values = custom_colors) +
  scale_x_continuous(breaks = c(10, 20, 30)) +
  labs(
    x = "$m$",
    y = "$\\widehat{\\mathbb{E}}(|\\widehat{\\mathcal{S}}|/|\\mathcal{S}|)$"
  ) +
  academic_theme

plot_1_2 <- p_prob_12 / p_size_12 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# ==========================================
# 4. PLOT B: SCENARIO III (Isolated)
# ==========================================

df_3 <- df_scen3 %>%
  mutate(
    Method = factor(Method, levels = c("LASSO", "LASSO+FTBS", "FTBS")),
    K_label = factor(K, levels = c(6, 10, 20, 50), labels = paste("K =", c(6, 10, 20, 50)))
  )

p_prob_3 <- ggplot(df_3, aes(x = m, y = Prob, shape = Method, linetype = Method, color = Method, group = Method)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  facet_wrap(~ K_label, nrow = 1) +
  scale_linetype_manual(name = NULL, values = custom_lines) +
  scale_shape_manual(name = NULL, values = custom_shapes) +
  scale_color_manual(name = NULL, values = custom_colors) +
  scale_x_continuous(breaks = c(10, 20, 30)) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
  labs(
    x = "$m$",
    y = "$\\widehat{\\mathbb{P}}(\\boldsymbol{\\nu}_{\\rm cyclic} \\in {\\rm span}(\\widehat{\\mathcal{S}}))$"
  ) +
  academic_theme

p_size_3 <- ggplot(df_3, aes(x = m, y = Size, shape = Method, linetype = Method, color = Method, group = Method)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  facet_wrap(~ K_label, nrow = 1) +
  scale_linetype_manual(name = NULL, values = custom_lines) +
  scale_shape_manual(name = NULL, values = custom_shapes) +
  scale_color_manual(name = NULL, values = custom_colors) +
  scale_x_continuous(breaks = c(10, 20, 30)) +
  labs(
    #title = "B) Expected Relative Model Size (Scenario III)",
    x = "$m$",
    y = "$\\widehat{\\mathbb{E}}(|\\widehat{\\mathcal{S}}|/|\\mathcal{S}|)$"
  ) +
  academic_theme

plot_3 <- p_prob_3 / p_size_3 + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# ==========================================
# 5. TIKZ DEVICE EXPORT (Strict LaTeX Widths)
# ==========================================
# Adjusted width to 6.2 inches to fit within a standard \textwidth
# Adjusted height to maintain a clean aspect ratio

cat("Generating Plot_Scenarios_1_and_2.tex...\n")
tikz(file = "Plot_Scenarios_1_and_2.tex", width = 6.2, height = 6.5)
print(plot_1_2)
invisible(dev.off())

cat("Generating Plot_Scenario_3.tex...\n")
tikz(file = "Plot_Scenario_3.tex", width = 6.2, height = 4.0)
print(plot_3)
invisible(dev.off())

cat("Done! The plots are now constrained to standard LaTeX text widths.\n")