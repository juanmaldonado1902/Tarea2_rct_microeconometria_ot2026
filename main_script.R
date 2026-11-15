# Tarea 2 - Microeconometria Aplicada (ITAM): RCT
# Creador : Juan Pablo Maldonado
# Revisor:

library(tidyverse)
library(fixest)        
library(sandwich)      
library(car)           
library(knitr)

# Correr desde la carpeta del proyecto (donde esta ggl2008.csv).
walk(c("figuras", "tablas"), dir.create, showWarnings = FALSE)
theme_set(theme_minimal(base_size = 11))

# Los grupos se ordenan de menor a mayor presion social
ggl <- read_csv("ggl2008.csv", show_col_types = FALSE) %>%
  mutate(arm = fct_reorder(arm, intensity))
rct <- ggl %>% filter(!is.na(treat))    # Neighbors (1) vs. control (0)

glimpse(ggl)
head(ggl)

covs   <- c("female", "age", "g2000", "g2002", "p2000", "p2002", "p2004", "hh_size")
grupos <- c("Control", "Civic Duty", "Hawthorne", "Self", "Neighbors")

# Tabla de LaTeX con caption y etiqueta 
options(knitr.kable.NA = "")
a_latex <- function(d, archivo, caption, label, digitos = 4, ...)
  writeLines(kable(d, "latex", booktabs = TRUE, digits = digitos, linesep = "",
                   caption = caption, label = label, format.args = list(big.mark = ","), ...),
             archivo)

# Tabla 1:asignacion por hogar y prueba F de igualdad entre grupos.
hogares <- ggl %>%
  summarise(arm = first(arm), across(all_of(covs), mean), .by = hh_id)

tabla1 <- hogares %>%
  pivot_longer(all_of(covs), names_to = "variable") %>%
  summarise(media = mean(value), .by = c(variable, arm)) %>%
  pivot_wider(names_from = arm, values_from = media, names_sort = TRUE) %>%
  mutate(p_value = map_dbl(variable, ~ anova(lm(reformulate("arm", .x), hogares))[["Pr(>F)"]][1]))

# Renglon postratamiento: participacion individual en 2006 con IC (EE por hogar).
# Su p-value es la Wald de que las cuatro cartas tienen efecto cero.
m_grupos <- feols(voted ~ 0 + arm, ggl, cluster = ~hh_id)
voto <- tibble(arm = levels(ggl$arm), voto = coef(m_grupos), ee = se(m_grupos))
p_voto <- wald(feols(voted ~ arm, ggl, cluster = ~hh_id), print = FALSE)$p

print(tabla1, width = Inf)
print(voto)
print(c(p_voto = p_voto))

tabla1 %>%
  bind_rows(voto %>% select(arm, voto) %>% pivot_wider(names_from = arm, values_from = voto) %>%
              mutate(variable = "voted (2006)", p_value = p_voto)) %>%
  a_latex("tablas/tabla1_balance.tex", "Medias por grupo experimental (Tabla 1 del paper)",
          "balance", 3, col.names = c("Variable", grupos, "p-value"))

ggsave("figuras/fig1_participacion.pdf", width = 6.5, height = 4,
       ggplot(voto, aes(fct_inorder(arm), voto)) +
         geom_pointrange(aes(ymin = voto - 1.96 * ee, ymax = voto + 1.96 * ee)) +
         scale_x_discrete(labels = grupos) +
         labs(title = "Participación en la primaria de 2006 por grupo (IC al 95%)",
              x = NULL, y = "Tasa de participación"))

# treat sobre las covariables 
m1c <- feols(reformulate(covs, "treat"), rct, cluster = ~hh_id)
summary(m1c)
wald(m1c)    # Aquí la hipótesis nula es que todos los coefs son cero

# Neyman = MCO
y1 <- rct$voted[rct$treat == 1]
y0 <- rct$voted[rct$treat == 0]

neyman <- mean(y1) - mean(y0)
ee_ney <- sqrt(var(y1) / length(y1) + var(y0) / length(y0))
m2b    <- lm(voted ~ treat, rct)

# El EE de Neyman es identico al robusto HC2 de la regresion
print(c(media_trat = mean(y1), media_control = mean(y0), neyman = neyman,
        mco = coef(m2b)[["treat"]], ee_neyman = ee_ney, ee_hc2 = sqrt(vcovHC(m2b, "HC2")[2, 2]),
        ee_cluster = sqrt(vcovCL(m2b, cluster = ~hh_id)[2, 2]), n1 = length(y1), n0 = length(y0)))

# Tabla 3 del paper, columnas (a) y (b)
m2c <- list(
  "(a)" = feols(voted ~ civicduty + hawthorne + self + neighbors, ggl, cluster = ~hh_id),
  "(b)" = feols(voted ~ civicduty + hawthorne + self + neighbors | cluster, ggl, cluster = ~hh_id))
etable(m2c, digits = 4)
etable(m2c, digits = 4, digits.stats = 4, tex = TRUE, replace = TRUE,
       file = "tablas/tabla2_cartas.tex", fitstat = ~ n + r2,
       title = "Efecto de cada carta sobre la participación (Tabla 3 del paper)",
       label = "tab:cartas", placement = "H",
       notes = "Errores estándar agrupados por hogar entre paréntesis. $^{*}$ $p<0.1$, $^{**}$ $p<0.05$, $^{***}$ $p<0.01$.",
       style.tex = style.tex("aer", model.format = "(a)", yesNo = c("Sí", "No"),
                             fixef.title = "", fixef.suffix = "", stats.title = " ", tablefoot = FALSE),
       dict = c(voted = "voted", civicduty = "Civic Duty", hawthorne = "Hawthorne",
                self = "Self", neighbors = "Neighbors", cluster = "EF por celda de asignación",
                `(Intercept)` = "Constante"))
read_lines("tablas/tabla2_cartas.tex") %>% str_replace("Observations", "Observaciones") %>%
  write_lines("tablas/tabla2_cartas.tex")

# Pruebas del paquete car de la pregunta 2d)
linearHypothesis(m2c$`(a)`, "hawthorne = civicduty")
linearHypothesis(m2c$`(a)`, "neighbors = self")

#presion social como "dosis" 
m_int <- feols(reformulate(covs, "intensity"), ggl, cluster = ~hh_id)
wald(m_int)                       

m_lin <- feols(voted ~ intensity, ggl, cluster = ~hh_id)
summary(m_lin)

# Linealidad, la idea es que si fuera lineal, las dummies de 1, 2 y 3 sobrarian
m_sat <- feols(voted ~ intensity + i(intensity, ref = c(0, 4)), ggl, cluster = ~hh_id)
wald(m_sat, keep = "::")

# diferencia de medias por estrato de historial (prior_share)
ggsave("figuras/fig2_estratos.pdf", width = 6.5, height = 4,
       rct %>%
         summarise(voto = mean(voted), .by = c(prior_share, treat)) %>%
         mutate(grupo = if_else(treat == 1, "Neighbors", "Control")) %>%
         ggplot(aes(prior_share, voto, colour = grupo)) +
         geom_line() + geom_point() +
         scale_x_continuous(labels = scales::percent) +
         labs(title = "Participación en 2006 según el historial de voto",
              x = "prior_share", y = "Tasa de participación", colour = NULL))

# ATE estratificado: diferencias por estrato ponderadas por n_s / N
estratos <- rct %>%
  summarise(n = n(), y1 = mean(voted[treat == 1]), y0 = mean(voted[treat == 0]),
            ee = sqrt(var(voted[treat == 1]) / sum(treat) + var(voted[treat == 0]) / sum(1 - treat)),
            de = sd(voted), .by = prior_share) %>%
  arrange(prior_share) %>%
  mutate(dif = y1 - y0, peso = n / sum(n))

ate_estratos <- sum(estratos$peso * estratos$dif)
ee_estratos  <- sqrt(sum(estratos$peso^2 * estratos$ee^2))
print(estratos)
print(c(ate_estratificado = ate_estratos, ee = ee_estratos, neyman = neyman))

estratos %>%
  transmute(prior_share = str_c(100 * prior_share, "\\%"), n, y1, y0, dif, ee, peso) %>%
  bind_rows(tibble(prior_share = "ATE estratificado", n = sum(estratos$n),
                   dif = ate_estratos, ee = ee_estratos, peso = 1)) %>%
  a_latex("tablas/tabla3_estratos.tex", "Diferencia de medias por estrato de historial de voto",
          "estratos", col.names = c("prior\\_share", "$n$", "$\\bar Y_1$", "$\\bar Y_0$",
                                    "Diferencia", "EE", "Peso"), escape = FALSE)


# Genero la misma estimación en una sola regresión
# sin constante  ( epresentado con el cero). El coeficiente treat: estrato es la diferencia de medias de ese estrato. 
m3b <- feols(voted ~ 0 + i(prior_share) + i(prior_share, treat), rct, vcov = "hetero")
coeftable(m3b)

tibble(Coeficiente = c(str_glue("$\\alpha_{{{estratos$prior_share * 100}\\%}}$"),
                       str_glue("$\\tau_{{{estratos$prior_share * 100}\\%}}$")),
       Estimador = coef(m3b), EE = se(m3b)) %>%
  a_latex("tablas/tabla4_estratos_reg.tex", "Regresión con interacciones estrato $\\times$ treat",
          "estratos_reg", escape = FALSE)

# Con las dummies de estrato centradas en su media, el coeficiente de treat es
# directamente el ATE estratificado y su EE robusto HC2 es el de Neyman estratificado.
Dc  <- scale(model.matrix(~ factor(prior_share), rct)[, -1], scale = FALSE)
m3c <- lm(voted ~ treat * Dc, rct)
print(c(ate_regresion = coef(m3c)[["treat"]], ee_hc2 = sqrt(vcovHC(m3c, "HC2")["treat", "treat"]),
        ate_a_mano = ate_estratos))

# Poder estadístico: en total, p = proporcion de tratados, es = tau / sigma.
# Se usa la sigma del control, segun yo así lo habiamos trabajado la otra vez.
# De igual forma genero tabla de sensibilidad con las demás desviaciones.
poder <- function(n, es, alpha = .06, p = .5) {
  z <- qnorm(1 - alpha / 2)
  pnorm(es * sqrt(n * p * (1 - p)) - z) + pnorm(-es * sqrt(n * p * (1 - p)) - z)
}
n_req <- function(es, alpha = .06, psi = .85, p = .5)
  (qnorm(1 - alpha / 2) + qnorm(psi))^2 / (p * (1 - p) * es^2)

# 4.b: estrato de 3.a con mayor efecto
mayor <- slice_max(estratos, dif)
e_may <- filter(rct, prior_share == mayor$prior_share)
sigmas <- function(d) c(control = sd(d$voted[d$treat == 0]), general = sd(d$voted),
                        tratamiento = sd(d$voted[d$treat == 1]))

poder_tab <- tibble(sigma = names(sigmas(rct)),
                    s_prom = sigmas(rct),  es_prom = neyman / s_prom,    n_prom = n_req(es_prom),
                    s_may  = sigmas(e_may), es_may = mayor$dif / s_may,  n_may  = n_req(es_may))
print(mayor)
print(poder_tab)

poder_tab %>%
  mutate(sigma = c("Control (recomendada)", "General (Neighbors + control)", "Tratamiento"),
         across(starts_with("n_"), ceiling)) %>%
  a_latex("tablas/tabla5_poder.tex", "Tamaño de muestra para 85\\% de poder según la $\\sigma$ utilizada",
          "poder", col.names = c("$\\sigma$ utilizada", "$\\sigma$", "$\\tau/\\sigma$", "$n$",
                                 "$\\sigma$", "$\\tau/\\sigma$", "$n$"), escape = FALSE)
read_lines("tablas/tabla5_poder.tex") %>%
  str_replace("^\\\\toprule$", str_c("\\\\toprule\n & \\\\multicolumn{3}{c}{Efecto promedio (2.b)}",
                                       " & \\\\multicolumn{3}{c}{Estrato de mayor efecto (80\\\\%)} \\\\\\\\",
                                       "\n\\\\cmidrule(lr){2-4} \\\\cmidrule(lr){5-7}")) %>%
  write_lines("tablas/tabla5_poder.tex")

es  <- c("Efecto promedio" = poder_tab$es_prom[[1]], "Estrato de mayor efecto" = poder_tab$es_may[[1]])
n85 <- n_req(es)

ggsave("figuras/fig3_poder.pdf", width = 6.5, height = 4,
       crossing(n = 10:2000, escenario = names(es)) %>%
         mutate(psi = poder(n, es[escenario])) %>%
         ggplot(aes(n, psi, colour = escenario)) +
         geom_line() +
         geom_hline(yintercept = .85, linetype = "dashed") +
         geom_vline(xintercept = n85, linetype = "dotted") +
         labs(title = bquote(paste("Curva de poder: ", alpha, " = 0.06, ", delta, " = ", .(round(es[[1]], 3)),
                                 " y ", .(round(es[[2]], 3)), ", p = 0.5")),
              x = "Tamaño muestral total (n)", y = quote("Poder estadístico " * psi(n)), colour = NULL) +
         theme(legend.position = "bottom"))

# efecto minimo detectable (5%, poder 80%), en puntos porcentuales
mde <- function(n, p, sigma = poder_tab$s_prom[[1]])
  100 * (qnorm(.975) + qnorm(.8)) * sigma / sqrt(n * p * (1 - p))
c(real = mde(nrow(rct), mean(rct$treat)), mitad_y_mitad = mde(nrow(rct), .5))

# tasas de atricion por grupo
rct <- rct %>% mutate(observado = as.integer(!is.na(voted_obs)))
rct %>% summarise(n = n(), observados = sum(observado), atricion = 1 - mean(observado), .by = treat)

# con monotonicidad se recorta, en el grupo que observa mas, p = (q1 - q0) / q1.
# Con Y binaria se recortan unos (cota inferior) o ceros (cota superior).
lee <- function(d) {
  q <- tapply(d$observado, d$treat, mean)
  y <- tapply(d$voted_obs, d$treat, mean, na.rm = TRUE)
  p <- (q[["1"]] - q[["0"]]) / q[["1"]]     # aqui el tratamiento observa mas
  c(q0 = q[["0"]], q1 = q[["1"]], p = p, y1 = y[["1"]], y0 = y[["0"]], ingenuo = y[["1"]] - y[["0"]],
    inferior = max(0, (y[["1"]] - p) / (1 - p)) - y[["0"]],
    superior = min(1, y[["1"]] / (1 - p)) - y[["0"]])
}

(lee_global <- lee(rct))
tibble(Estimador = c("Sin corrección (muestra observada)", "Cota inferior de Lee",
                     "Cota superior de Lee", "Neyman sin atrición (2.b)"),
       Valor = c(lee_global[c("ingenuo", "inferior", "superior")], neyman)) %>%
  a_latex("tablas/tabla6_lee.tex", "Efecto de Neighbors con atrición: cotas de Lee", "lee")

# cotas por celda edad x tamaño del hogar, ponderadas por los
# siempre observados de cada celda (n_c * q0_c)
celdas <- rct %>%
  mutate(hh3 = pmin(hh_size, 3)) %>%
  nest(.by = c(age_group, hh3)) %>%
  mutate(n = map_int(data, nrow), cotas = map(data, lee)) %>%
  unnest_wider(cotas) %>%
  mutate(peso = n * q0 / sum(n * q0)) %>%
  arrange(fct_relevel(age_group, "jovenes", "ideal"), hh3)
print(select(celdas, -data))

(lee_cov <- celdas %>% summarise(across(c(inferior, superior), ~ sum(peso * .x))))

celdas %>%
  transmute(age_group, hh3 = if_else(hh3 == 3, "3+", as.character(hh3)), n, q0, q1, p,
            inferior, superior, peso) %>%
  bind_rows(tibble(age_group = "Agregado", inferior = lee_cov$inferior,
                   superior = lee_cov$superior, peso = 1)) %>%
  a_latex("tablas/tabla7_lee_celdas.tex", "Cotas de Lee por celda de edad y tamaño del hogar",
          "lee_celdas", 3, col.names = c("Edad", "Hogar", "$n$", "$q_0$", "$q_1$", "$p$",
                                         "Inferior", "Superior", "Peso"), escape = FALSE)
