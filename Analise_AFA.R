################################### CFC (Análise presente na TESE) #################################
AFA_AOV <- DADOS_TESE <- read_excel("G:/DADOS_TESE.xlsx", sheet = "AFA")

AFA_AOV$LAGARTO <- as.factor(AFA_AOV$LAGARTO)
AFA_AOV$GRUPO <- as.factor(AFA_AOV$GRUPO)
AFA_AOV$CONDIÇÃO <- as.factor(AFA_AOV$CONDIÇÃO)
AFA_AOV$DIAS <- as.factor(AFA_AOV$DIAS)
AFA_AOV$MEMORIA <- as.factor(AFA_AOV$MEMORIA)

######### ANÁLISES DBCDM

# Passo 1: construir modelo linear
modelo <- aov_car(DBCDM ~ GRUPO * Dias * Condicao * MEMORIA+Error(LAGARTO/Dias * Condicao), 
                  data = AFA_AOV)
# Passo 2: ANOVA mista


Resultado <- aov_ez("LAGARTO", dv = "DBCDM", data = AFA_AOV, 
                    between = c("GRUPO","MEMORIA"), within = c("Dias","Condicao"),
                    observed = "GRUPO", include_aov = T,
                    anova_table = list(p_adjust_method = "none"))
Resultado
afex_plot(Resultado, "Condicao", "MEMORIA", c("GRUPO","Dias"), error = "mean")

effectsize(model = modelo)

post_hoc_GRUPO_Dias_Condicao <- emmeans(modelo, pairwise ~ Dias*Condicao, adjust = "BONFERRONI")
resumo <- summary(post_hoc_GRUPO_Dias_Condicao)
p_values <- resumo$contrasts[, "p.value"]


# Filter the results for p-values < 0.05
significant_results <- resumo$contrasts[p_values < 0.05,]

######### ANÁLISES DBCL
PerformanceAnalytics::chart.Correlation(AFA_AOV[,8:27], method = )

# Passo 1: construir modelo linear
modelo <- aov_car(DBCL ~ GRUPO * MEMORIA * Dias * Condicao +Error(LAGARTO/Dias * Condicao), 
                  data = AFA_AOV, contr.sum(c("GRUPO","MEMORIA"), contrasts = T))
# Passo 2: ANOVA mista
modelo

Resultado <- aov_ez("LAGARTO", dv = "DBCL", data = AFA_AOV, 
                    between = c("GRUPO","MEMORIA"), within = c("Dias","Condicao"),
                    observed = "GRUPO", include_aov = T,
                    anova_table = list(p_adjust_method = "none"))
Resultado
afex_plot(Resultado, "Condicao", "MEMORIA", c("GRUPO","Dias"), error = "mean")

effectsize(model = modelo)

post_hoc_GRUPO_Dias_Condicao <- emmeans(modelo, pairwise ~ GRUPO*Dias*Condicao*MEMORIA, adjust = "BONFERRONI")
resumo <- summary(post_hoc_GRUPO_Dias_Condicao)
p_values <- resumo$contrasts[, "p.value"]

# Filter the results for p-values < 0.05
significant_results <- resumo$contrasts[p_values < 0.05,]

######### ANÁLISES D-Gm

# Passo 1: construir modelo linear
modelo <- aov_car(DGCDM ~ GRUPO * Dias * Condicao * MEMORIA+Error(LAGARTO/Dias * Condicao), 
                  data = AFA_AOV)
# Passo 2: ANOVA mista


Resultado <- aov_ez("LAGARTO", dv = "DGCDM", data = AFA_AOV, 
                    between = c("GRUPO","MEMORIA"), within = c("Dias","Condicao"),
                    observed = "GRUPO", include_aov = T,
                    anova_table = list(p_adjust_method = "none"))
Resultado
afex_plot(Resultado, "Condicao", "MEMORIA", c("GRUPO","Dias"), error = "mean")

effectsize(model = modelo)

post_hoc_GRUPO_Dias_Condicao <- emmeans(modelo, pairwise ~ Condicao, adjust = "BONFERRONI")
resumo <- summary(post_hoc_GRUPO_Dias_Condicao)
p_values <- resumo$contrasts[, "p.value"]

# Filter the results for p-values < 0.05
significant_results <- resumo$contrasts[p_values < 0.05,]

######### ANÁLISES DGCL

# Passo 1: construir modelo linear
modelo <- aov_car(DGCL ~ GRUPO * Dias * Condicao * MEMORIA+Error(LAGARTO/Dias * Condicao), 
                  data = AFA_AOV)
# Passo 2: ANOVA mista


Resultado <- aov_ez("LAGARTO", dv = "DGCL", data = AFA_AOV, 
                    between = c("GRUPO","MEMORIA"), within = c("Dias","Condicao"),
                    observed = "GRUPO", include_aov = T,
                    anova_table = list(p_adjust_method = "none"))
Resultado
afex_plot(Resultado, "Condicao", "MEMORIA", c("GRUPO","Dias"), error = "mean")

effectsize(model = modelo)

post_hoc_GRUPO_Dias_Condicao <- emmeans(modelo, pairwise ~ GRUPO*MEMORIA*Condicao*Dias, adjust = "bonferroni")
resumo <- summary(post_hoc_GRUPO_Dias_Condicao)
p_values <- resumo$contrasts[, "p.value"]

# Filter the results for p-values < 0.05
significant_results <- resumo$contrasts[p_values < 0.05,]

#################### Potencia
DADOS_ANALISE_AFA_LAGARTOS2 <- read_excel("G:/Meu Drive/Doutorado/Tese/DADOS_ANALISE_AFA_LAGARTOS.xlsx", 
                                          sheet = "POTENCIA", range = "A1:K241")
View(DADOS_ANALISE_AFA_LAGARTOS2)
POT <- DADOS_ANALISE_AFA_LAGARTOS2
POT <- as.data.frame(POT)
names(POT) <- c("Lagarto", "Grupo", "POTA", "POTB", "POTC",
                "POTD", "POTE", "NUCLEO", "Condicao", "Dias", "MEMORIA")
#FATORES
POT$Grupo <- as.factor(POT$Grupo)
POT$Lagarto <- as.factor(POT$Lagarto)
POT$NUCLEO <- as.factor(POT$NUCLEO)
POT$Condicao <- as.factor(POT$Condicao)
POT$Dias <- as.factor(POT$Dias)
POT$MEMORIA <- as.factor(POT$MEMORIA)

# Passo 1: construir modelo linear Para Condicao mnemônica
POT %>% group_by(NUCLEO) %>% filter(NUCLEO == "CDM", MEMORIA == "AQUISICAO") -> CDM

aov_car(POTA ~ Grupo+Condicao+Dias + Error(Lagarto/Condicao+Dias), 
        fun_aggregate = mean, data = CDM)

POT %>% group_by(NUCLEO) %>% filter(NUCLEO == "CL", MEMORIA == "AQUISICAO") -> CL

aov_car(POTA ~ Grupo*Condicao*Dias + Error(Lagarto/Condicao*Dias), 
        fun_aggregate = mean, data = CL)  

POT %>% group_by(NUCLEO) %>% filter(NUCLEO == "CDM", MEMORIA == "CONSOLIDACAO") -> CDM

aov_car(POTA ~ Grupo*Condicao*Dias + Error(Lagarto/Condicao*Dias), 
        fun_aggregate = mean, data = CDM)

POT %>% group_by(NUCLEO) %>% filter(NUCLEO == "CL", MEMORIA == "CONSOLIDACAO") -> CL

aov_car(POTA ~ Grupo*Condicao*Dias + Error(Lagarto/Condicao*Dias), 
        fun_aggregate = mean, data = CL)

POT %>% group_by(NUCLEO) %>% filter(NUCLEO == "CDM", MEMORIA == "EVOCACAO") -> CDM

aov_car(POTA ~ Grupo*Condicao*Dias + Error(Lagarto/Condicao*Dias), 
        fun_aggregate = mean, data = CDM)

POT %>% group_by(NUCLEO) %>% filter(NUCLEO == "CL", MEMORIA == "EVOCACAO") -> CL

aov_car(POTA ~ Grupo*Condicao*Dias + Error(Lagarto/Condicao*Dias), 
        fun_aggregate = mean, data = CL)

# Passo 1: construir modelo linear Para Condicao mnemônica

POT %>% group_by(NUCLEO) %>% filter(NUCLEO == "CDM") -> CDM

aov_car(POTA ~ Grupo*MEMORIA + Error(Lagarto/Condicao*Dias), 
        fun_aggregate = mean, data = CDM)

POT %>% group_by(NUCLEO) %>% filter(NUCLEO == "CL") -> CL

aov_car(POTA ~ Grupo*MEMORIA + Error(Lagarto/Condicao*Dias), 
        fun_aggregate = mean, data = CL)


# Passo 2: ANOVA mista

Resultado <- aov_ez("Lagarto", dv = "POTA", data = CDM, 
                    between = c("Grupo","MEMORIA"), within = c("Condicao","Dias"),
                    observed = "Dias",
                    fun_aggregate = mean,
                    include_aov = T,
                    anova_table = list(p_adjust_method = "bonferroni"))
Resultado
afex_plot(Resultado, "Teste", "Grupo", c("Dia"), error = "mean")

Posteste <- PLV_CFC %>% group_by(Grupo, Dia) %>%
  emmeans_test(DB ~ Teste, p.adjust.method = "bonf")
View(Posteste)


