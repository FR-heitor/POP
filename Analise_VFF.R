############################# VFF (Análise presente na TESE)##################################### 
VFF_AOV <- DADOS_TESE <- read_excel("G:/DADOS_TESE.xlsx", sheet = "VFF")

VFF_AOV$LAGARTO <- as.factor(VFF_AOV$LAGARTO)
VFF_AOV$GRUPO <- as.factor(VFF_AOV$GRUPO)
VFF_AOV$CONDIÇÃO <- as.factor(VFF_AOV$CONDIÇÃO)
VFF_AOV$DIAS <- as.factor(VFF_AOV$DIAS)
VFF_AOV$MEMORIA <- as.factor(VFF_AOV$MEMORIA)


######### ANÁLISES D-B

# Passo 1: construir modelo linear
modelo <- aov_car(DBPVL ~ GRUPO*CONDIÇÃO*DIAS*MEMORIA  +Error(LAGARTO/DIAS*CONDIÇÃO), 
                  data = VFF_AOV)
# Passo 2: ANOVA mista


Resultado <- aov_ez("LAGARTO", dv = "DBPVL", data = VFF_AOV, 
                    between = c("GRUPO","MEMORIA"), within = c("DIAS","CONDIÇÃO"),
                    observed = "GRUPO", include_aov = T,
                    anova_table = list(p_adjust_method = "holm"))
Resultado
afex_plot(Resultado, "CONDIÇÃO", "GRUPO", c("DIAS"), error = "TUKEY")
effectsize(model = modelo)
post_hoc_GRUPO_DIAS_CONDIÇÃO <- emmeans(modelo, pairwise ~ GRUPO * DIAS * CONDIÇÃO * MEMORIA, adjust = "tukey")
resumo <- summary(post_hoc_GRUPO_DIAS_CONDIÇÃO)
p_values <- resumo$contrasts[, "p.value"]

# Filter the results for p-values < 0.05
significant_results <- resumo$contrasts[p_values < 0.05,]
######### ANÁLISES T-B

# Passo 1: construir modelo linear
modelo <- aov_car(TBPVL ~ GRUPO*CONDIÇÃO*DIAS*MEMORIA  +Error(LAGARTO/DIAS*CONDIÇÃO), 
                  data = VFF_AOV)
# Passo 2: ANOVA mista


Resultado <- aov_ez("LAGARTO", dv = "TBPVL", data = VFF_AOV, 
                    between = c("GRUPO","MEMORIA"), within = c("DIAS","CONDIÇÃO"),
                    observed = "GRUPO", include_aov = T,
                    anova_table = list(p_adjust_method = "holm"))
Resultado
afex_plot(Resultado, "CONDIÇÃO", "GRUPO", c("DIAS","MEMORIA"), error = "TUKEY")
effectsize(model = modelo)
post_hoc_GRUPO_DIAS_CONDIÇÃO <- emmeans(modelo, pairwise ~ GRUPO * DIAS * CONDIÇÃO * MEMORIA, adjust = "TUKEY")
resumo <- summary(post_hoc_GRUPO_DIAS_CONDIÇÃO)
p_values <- resumo$contrasts[, "p.value"]

# Filter the results for p-values < 0.05
significant_results <- resumo$contrasts[p_values < 0.05,]

######### ANÁLISES D-Gm

# Passo 1: construir modelo linear
modelo <- aov_car(DGmPVL ~ GRUPO*CONDIÇÃO*DIAS*MEMORIA  +Error(LAGARTO/DIAS*CONDIÇÃO), 
                  data = VFF_AOV)
# Passo 2: ANOVA mista


Resultado <- aov_ez("LAGARTO", dv = "DGmPVL", data = VFF_AOV, 
                    between = c("GRUPO","MEMORIA"), within = c("DIAS","CONDIÇÃO"),
                    observed = "GRUPO", include_aov = T,
                    anova_table = list(p_adjust_method = "holm"))
Resultado
afex_plot(Resultado, "CONDIÇÃO", "GRUPO", c("DIAS","MEMORIA"), error = "TUKEY")
effectsize(model = modelo)
post_hoc_GRUPO_DIAS_CONDIÇÃO <- emmeans(modelo, pairwise ~ GRUPO * DIAS * CONDIÇÃO * MEMORIA, adjust = "tukey")
resumo <- summary(post_hoc_GRUPO_DIAS_CONDIÇÃO)
p_values <- resumo$contrasts[, "p.value"]

# Filter the results for p-values < 0.05
significant_results <- resumo$contrasts[p_values < 0.05,]

######### ANÁLISES T-Gm

# Passo 1: construir modelo linear
modelo <- aov_car(TGmPVL ~ GRUPO*CONDIÇÃO*DIAS*MEMORIA  +Error(LAGARTO/DIAS*CONDIÇÃO), 
                  data = VFF_AOV)
# Passo 2: ANOVA mista


Resultado <- aov_ez("LAGARTO", dv = "TGmPVL", data = VFF_AOV, 
                    between = c("GRUPO","MEMORIA"), within = c("DIAS","CONDIÇÃO"),
                    observed = "GRUPO", include_aov = T,
                    anova_table = list(p_adjust_method = "holm"))
Resultado
afex_plot(Resultado, "CONDIÇÃO", "GRUPO", c("DIAS","MEMORIA"), error = "TUKEY")
effectsize(model = modelo)

post_hoc_GRUPO_DIAS_CONDIÇÃO <- emmeans(modelo, pairwise ~ GRUPO * DIAS * CONDIÇÃO * MEMORIA, adjust = "TUKEY")
resumo <- summary(post_hoc_GRUPO_DIAS_CONDIÇÃO)
p_values <- resumo$contrasts[, "p.value"]

# Filter the results for p-values < 0.05
significant_results <- resumo$contrasts[p_values < 0.05,]

boxplot( formula = TGmPVL ~ GRUPO * DIAS * CONDIÇÃO * MEMORIA, data = VFF_AOV)

######### ANÁLISES D-GM

# Passo 1: construir modelo linear
modelo <- aov_car(DGMPVL ~ GRUPO*CONDIÇÃO*DIAS*MEMORIA  +Error(LAGARTO/DIAS*CONDIÇÃO), 
                  data = VFF_AOV)
# Passo 2: ANOVA mista


Resultado <- aov_ez("LAGARTO", dv = "DGMPVL", data = VFF_AOV, 
                    between = c("GRUPO","MEMORIA"), within = c("DIAS","CONDIÇÃO"),
                    observed = "GRUPO", include_aov = T,
                    anova_table = list(p_adjust_method = "holm"))
Resultado
afex_plot(Resultado, "CONDIÇÃO", "GRUPO", c("DIAS"), error = "TUKEY")
effectsize(model = modelo)

post_hoc_GRUPO_DIAS_CONDIÇÃO <- emmeans(modelo, pairwise ~ GRUPO * DIAS * CONDIÇÃO * MEMORIA, adjust = "TUKEY")
resumo <- summary(post_hoc_GRUPO_DIAS_CONDIÇÃO)
p_values <- resumo$contrasts[, "p.value"]

# Filter the results for p-values < 0.05
significant_results <- resumo$contrasts[p_values < 0.05,]

######### ANÁLISES T-GM

# Passo 1: construir modelo linear
modelo <- aov_car(TGMPVL ~ GRUPO * DIAS * CONDIÇÃO * MEMORIA  +Error(LAGARTO/DIAS*CONDIÇÃO), 
                  data = VFF_AOV)
# Passo 2: ANOVA mista


Resultado <- aov_ez("LAGARTO", dv = "TGMPVL", data = VFF_AOV, 
                    between = c("GRUPO","MEMORIA"), within = c("DIAS","CONDIÇÃO"),
                    observed = "GRUPO", include_aov = T)
Resultado
afex_plot(Resultado, "CONDIÇÃO", "GRUPO", c("DIAS"), error = "TUKEY")
effectsize(model = modelo)

post_hoc_GRUPO_DIAS_CONDIÇÃO <- emmeans(modelo, pairwise ~ GRUPO * DIAS * CONDIÇÃO * MEMORIA, adjust = "TUKEY")
resumo <- summary(post_hoc_GRUPO_DIAS_CONDIÇÃO)
p_values <- resumo$contrasts[, "p.value"]

# Filter the results for p-values < 0.05
significant_results <- resumo$contrasts[p_values < 0.05,]

########## CORRELAÇÃO
library(psych)
library(cocor)
R <- corr.test(PLV_CFC[1:40,6:15])
p <- R$p
r <- R$r
corPlot(r,numbers=TRUE,diag=FALSE,stars=TRUE, pval = p, main="Correlation plot
 with Holm corrected 'significance'")

