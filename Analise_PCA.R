####### final- PCA (Resumir variáveis em torno da medida) #################
# Carregar bibliotecas
library(readxl)
library(FactoMineR)
library(factoextra)
library(tidyverse)

# Carregar banco de dados
PCA  <- read_excel("G://DADOS_ANALISE_PCA_LAGARTOS.xlsx", 
                   sheet = "PCA")
# Definir variáveis independentes
PCA$LAGARTO <- as.factor(PCA$LAGARTO)
PCA$GRUPO <- as.factor(PCA$GRUPO)
PCA$CODICAO <- as.factor(PCA$CODICAO)
PCA$DIAS <- as.factor(PCA$DIAS)
PCA$MEMORIA <- as.factor(PCA$MEMORIA)
PCA$INTERACAO <- as.factor(PCA$INTERACAO)

PCA <- as.data.frame(PCA)

# Realizando a Análise de Componentes Principais (PCA) nos dados das colunas 3 a 29. 
# Os dados são escalonados para ter média 0 e desvio padrão 1.
# Três dimensões principais são mantidas na análise. 
# As colunas 1 a 5 são mantidas como variáveis suplementares.
res.pca <- PCA(PCA[,c(3:29)], scale.unit = T, ncp = 3, quali.sup = c(1:5))

# Extraindo e visualizando os autovalores (variações) associados a cada componente principal.
eig.val <- get_eigenvalue(res.pca)
fviz_eig(res.pca)
head(eig.val)

# Extraindo os resultados para as variáveis na análise PCA.
var <- get_pca_var(res.pca)
head(var$coord)
head(var$contrib)
head(var$cor)
head(var$cos2)

# Visualizando as variáveis da análise PCA.
# As cores são determinadas pelo quadrado do cosseno das variáveis.
fviz_pca_var(res.pca, col.var = "cos2", 
             axes = c(1,2), 
             palette = c("ggplot2"), 
             legend.title = "cos2", 
             fill.var = "white")

# Agrupando as variáveis em 3 grupos usando o algoritmo k-means.
set.seed(123)
res.km <- kmeans(var$coord, centers = 3, nstart = 1)
grp <- as.factor(res.km$cluster)

# Extraindo os resultados para os indivíduos na análise PCA.
ind <- get_pca_ind(res.pca)
head(ind$coord)
head(ind$contrib)
head(ind$cos2)

# Descrevendo as dimensões na análise PCA.
res.desc <- dimdesc(res.pca, axes = c(1:3), proba = 0.05)

# Visualizando a matriz de correlação entre as dimensões e os constructos.
corrplot(var$cor, is.corr=T, 
         title = "Corelação das dimensões com os constructos", 
         col = col4(150), 
         tl.col = "black")

# Visualizando o quadrado do cosseno das variáveis para as duas primeiras dimensões.
fviz_cos2(res.pca, choice = "var", axes = 1)
fviz_cos2(res.pca, choice = "var", axes = 2)

# Visualizando a contribuição de cada variável para as três primeiras dimensões.
fviz_contrib(res.pca, "var", axes = 1)
fviz_contrib(res.pca, "var", axes = 2)
fviz_contrib(res.pca, "var", axes = 3)

############################## PLOT - PCA #################
# Create the PCA biplot
fviz_pca_biplot(res.pca, 
                axes = c(1,2),
                # Individuals
                geom.ind = "point",
                fill.ind = PCA$GRUPO, col.ind = 3, col.ind.sup = "black",
                pointshape = 21, pointsize = 2,
                palette = c("darkred","red","orange","yellow"),
                addEllipses = T,
                ellipse.type = "confidence",
                # Variables
                alpha.var ="contrib", col.var = "cos2",
                gradient.cols = c("black", "darkred", "darkorange"),
                repel = T,
                
                legend.title = list(fill = "CODICAO", color = "Correlation",
                                    alpha = "Quality")
)


summary.PCA(res.pca, nbelements = 50, )
esquisser(PCA)              


library(ggplot2)

ggplot(PCA) +
  aes(x = CODICAO, fill = GRUPO, y = DGMPVL) +
  geom_boxplot(position = "dodge") +
  scale_fill_hue(direction = 1) +
  labs(x = "TAREFA", y = "VFF", title = "INTERAÇÃO ENTRE CL --> CDM") +
  ggthemes::theme_base() +
  theme(legend.position = "top") +
  facet_grid(MEMORIA~DIAS, scales = "free")


