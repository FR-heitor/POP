################################### FDA (Comparação de anova fatorial) ####################################
#bibliotecas
library(readxl)
library(tidyverse)
library(effectsize)
library(emmeans)
library(rstatix)
library(stats)

#banco de dados
Exemplo <- read_excel("G:/FDA.xlsx")
View(Exemplo)

FDA <- Exemplo

FDA$LAGARTO <- as.factor(FDA$LAGARTO)
FDA$Grupo   <- as.factor(FDA$Grupo)     # CTL, MUS_1, MUS_2, MUS_3
FDA$Dias    <- as.factor(FDA$Dias)      # TREINO e TESTE

FDA$AF_1    <- as.numeric(FDA$AF_1)   # Exposição
FDA$AF_2    <- as.numeric(FDA$AF_2)   # Ambientação

# Filtra dados do treino
FDA %>% filter (Dias == "TREINO") -> FDA_Tr
# Passo 1: Modelo
modelo_fda <- aov(AF1_CL ~ Grupo, data = FDA_Tr) # Adaptar para Exposição ou Ambientação - aqui está em Exposição (AF_1)
# Passo 2: ANOVA mista
Resultado_fda <- anova_test(modelo_fda)
Resultado_fda
# Tamanho do efeito
effectsize(model = modelo_fda)
# Pós-teste
post_hoc_freeze <- emmeans(modelo_fda, pairwise ~ Grupo, adjust = "tukey")
resumo <- summary(post_hoc_freeze)
p_values <- resumo$contrasts[, "p.value"]

# Filtro dos resultados p-valor < 0.05
significant_results <- resumo$contrasts[p_values < 0.05,]

# Filtra dados do teste
FDA %>% filter (Dias == "TESTE") -> FDA_Tt
# Passo 1: Modelo
modelo_fda <- aov(AF1_CL ~ Grupo, data = FDA_Tt)  # Adaptar para Exposição ou Ambientação - aqui está em Exposição (AF_1)
# Passo 2: ANOVA mista
Resultado_fda <- anova_test(modelo_fda)
Resultado_fda
# Tamanho do efeito
effectsize(model = modelo_fda)
# Pós-teste
post_hoc_freeze <- emmeans(modelo_fda, pairwise ~ Grupo, adjust = "tukey")
resumo <- summary(post_hoc_freeze)
p_values <- resumo$contrasts[, "p.value"]
# Filtra valores de p < 0.05
significant_results <- resumo$contrasts[p_values < 0.05,]
#
