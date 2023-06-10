################################### FDA (Comparação de anova fatorial) ####################################
#bibliotecas
library(readxl)
library(afex)
library(tidyverse)
library(effectsize)
library(emmeans)
library(rstatix)
library(ez)

#banco de dados
Exemplo <- read_excel("G:/FDA.xlsx")
View(Exemplo)

FDA <- Exemplo

FDA$LAGARTO <- as.factor(FDA$LAGARTO)
FDA$COMPORTAMENTO <- as.factor(FDA$COMPORTAMENTO) # Treino e Teste
FDA$MEMORIA <- as.factor(FDA$MEMORIA)             # processo mnemônico: aquisição|consolidação|evocação
FDA$Grupo   <- as.factor(FDA$Grupo)               # CTL ou MUS

# Passo 1: Modelo

modelo_fda <- aov(AF1_CL ~ Grupo, data = FDA)
# Passo 2: ANOVA mista


Resultado_fda <- anova_test(modelo_fda)
Resultado_fda

#AF1_CL_GHT <- games_howell_test(FDA, AF1_CL ~ GRUPO+COMPORTAMENTO, detailed = T)
effectsize(model = modelo_fda)

post_hoc_freeze <- emmeans(modelo_fda, pairwise ~ Grupo, adjust = "tukey")
resumo <- summary(post_hoc_freeze)
p_values <- resumo$contrasts[, "p.value"]

# Filter the results for p-values < 0.05
significant_results <- resumo$contrasts[p_values < 0.05,]

# Passo 1: Modelo

modelo_fda <- aov(AF2_CL ~ Grupo, data = FDA)
# Passo 2: ANOVA mista


Resultado_fda <- anova_test(modelo_fda)
Resultado_fda


effectsize(model = modelo_fda)

post_hoc_freeze <- emmeans(modelo_fda, pairwise ~ Grupo, adjust = "tukey")
resumo <- summary(post_hoc_freeze)
p_values <- resumo$contrasts[, "p.value"]

# Filter the results for p-values < 0.05
significant_results <- resumo$contrasts[p_values < 0.05,]


#