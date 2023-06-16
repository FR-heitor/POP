################################## Comportamento (Apenas exposição ou ambientação) ############################

# bibliotecas
library(readxl)
library(tidyverse)
library(effectsize)
library(rstatix)
library(emmeans)
library(WRS2)
library(afex)


#bancao de dados
CONGELAMENTO <- read_excel("G:/Meu Drive/Doutorado/Tese/DADOS_ANALISE_PCA_LAGARTOS.xlsx", 
                           sheet = "CONGELAMENTO_QUALI", range = "a1:h67")
View(CONGELAMENTO)


Freeze <- CONGELAMENTO

Freeze$LAGARTO <- as.factor(Freeze$ANIMAL)
Freeze$Dias <- as.factor(Freeze$Dias)
Freeze$GRUPO <- as.factor(Freeze$GRUPO)
Freeze$MEMORIA <- as.factor(Freeze$MEMORIA)

Freeze <- remove_missing(Freeze, na.rm = T)
# Exposição
# Passo 1: Modelo
modelo_freeze <- aov_car(`C - EXP` ~ GRUPO*MEMORIA
                         +Error(LAGARTO/(Dias)), 
                         data = Freeze)
# Passo 2: ANOVA Robusta
selvagem1 <- bwtrim(formula = `C - EXP` ~ GRUPO+MEMORIA, id = LAGARTO,
         data = Freeze, tr = 0.1, nboot = 10000, MDIS = F, est = MOM) 
selvagem1

effectsize(modelo_freeze)
post_hoc_freeze <- emmeans(modelo_freeze, pairwise ~ GRUPO*MEMORIA, adjust = "tukey")
resumo <- summary(post_hoc_freeze)
p_values <- resumo$contrasts[, "p.value"]

# Filter the results for p-values < 0.05
significant_results <- resumo$contrasts[p_values < 0.05,]

# Ambientação
# Passo 1: Modelo
modelo_freeze <- aov_car(`C - AMB` ~ GRUPO*MEMORIA
                         +Error(LAGARTO/(Dias)), 
                         data = Freeze)
# Passo 2: ANOVA Robusta
selvagem2 <- bwtrim(formula = `C - AMB` ~ GRUPO+MEMORIA,
        data = Freeze, tr = 0.1,nboot = 5000, pro.dis = F, est = MOM)
selvagem2


effectsize(modelo_freeze)
post_hoc_freeze <- emmeans(modelo_freeze, pairwise ~ GRUPO*MEMORIA, adjust = "tukey")
resumo <- summary(post_hoc_freeze)
p_values <- resumo$contrasts[, "p.value"]

# Filter the results for p-values < 0.05
significant_results <- resumo$contrasts[p_values < 0.05,]

esquisse::esquisser(Freeze)


Freeze$Comportamento <- factor(Freeze$Dias, levels = c("TREINO","TESTE"))

library(ggplot2)
library(ggbeeswarm)
ggplot(Freeze) +
  aes(y = MEMORIA, x = `C - EXP`, fill = Dias) +
  geom_boxplot() +
  geom_beeswarm(dodge.width = 0.75) +
  scale_fill_manual(values = c(TREINO = "#FFFFFF", TESTE = "#464546")) +
  labs(x = "Tempo de congelamento (s)", 
       y = "Grupos") +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "top") +
  facet_wrap(vars(GRUPO), ncol = 2L)

ggplot(Freeze) +
  aes(y = MEMORIA, x = `C - AMB`, fill = Dias) +
  geom_boxplot() +
  geom_beeswarm(dodge.width = 0.75) +
  scale_fill_manual(values = c(TREINO = "#FFFFFF", TESTE = "#464546")) +
  labs(x = "Tempo de congelamento (s)", 
       y = "Grupos") +
  coord_flip() +
  theme_classic() +
  theme(legend.position = "top") +
  facet_wrap(vars(GRUPO), ncol = 2L)


#MOVIMENTOS
Freeze <- remove_missing(Freeze, na.rm = T)
# Ambientação
# Passo 1: Modelo
modelo_freeze <- aov_car( `MOVI - AMB` ~ GRUPO*MEMORIA
                         +Error(ID/(Dias)), 
                         data = Freeze)
# Passo 2: ANOVA Robusta
selvagem1 <- bwtrim(formula = `MOVI - AMB` ~ MEMORIA*GRUPO, id = ID,
                    data = Freeze, tr = 0.1, nboot = 5000, MDIS = F, est = MOM) 
selvagem1

effectsize(modelo_freeze)
post_hoc_freeze <- emmeans(modelo_freeze, pairwise ~ GRUPO*MEMORIA, adjust = "tukey")
resumo <- summary(post_hoc_freeze)
p_values <- resumo$contrasts[, "p.value"]

# Filter the results for p-values < 0.05
significant_results <- resumo$contrasts[p_values < 0.05,]

# Exposição
# Passo 1: Modelo
modelo_freeze <- aov_car(`MOVI - EXP` ~ GRUPO*MEMORIA
                         +Error(ID/(Dias)), 
                         data = Freeze)
# Passo 2: ANOVA Robusta
selvagem2 <- bwtrim(formula = `MOVI - EXP` ~ GRUPO+MEMORIA, id = ID,
                    data = Freeze, tr = 0.1,nboot = 5000, pro.dis = F, est = MOM)
selvagem2


effectsize(modelo_freeze)
post_hoc_freeze <- emmeans(modelo_freeze, pairwise ~ GRUPO*MEMORIA, adjust = "tukey")
resumo <- summary(post_hoc_freeze)
p_values <- resumo$contrasts[, "p.value"]

# Filter the results for p-values < 0.05
significant_results <- resumo$contrasts[p_values < 0.05,]

esquisse::esquisser(Freeze)


Freeze$Dias <- factor(Freeze$Dias, levels = c("TREINO","TESTE"))

library(ggplot2)
library(ggbeeswarm)
ggplot(Freeze) +
  aes(y = MEMORIA, x = `MOVI - AMB`, fill = Dias) +
  geom_boxplot() +
  geom_beeswarm(dodge.width = 0.75) +
  scale_fill_manual(values = c(TREINO = "#FFFFFF", TESTE = "#464546")) +
  labs(x = "Número de movimentos", 
       y = "Grupos") +
  xlim(-5,5)+
  coord_flip() +
  theme_classic() +
  theme(legend.position = "top") +
  facet_wrap(vars(GRUPO), ncol = 2L)


ggplot(Freeze) +
  aes(y = MEMORIA, x = `MOVI - EXP`, fill = Dias) +
  geom_boxplot() +
  geom_beeswarm(dodge.width = 0.75) +
  scale_fill_manual(values = c(TREINO = "#FFFFFF", TESTE = "#464546")) +
  labs(x = "Número de movimentos", 
       y = "Grupos") +
  xlim(-5,5)+
  coord_flip() +
  theme_classic() +
  theme(legend.position = "top") +
  facet_wrap(vars(GRUPO), ncol = 2L)
