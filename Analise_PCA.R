####### final- PCA (Resumir variáveis em torno da medida) #################
library(readxl)
PCA  <- read_excel("G://DADOS_ANALISE_PCA_LAGARTOS.xlsx", 
                   sheet = "PCA")

PCA$LAGARTO <- as.factor(PCA$LAGARTO)
PCA$GRUPO <- as.factor(PCA$GRUPO)
PCA$CODICAO <- as.factor(PCA$CODICAO)
PCA$DIAS <- as.factor(PCA$DIAS)
PCA$MEMORIA <- as.factor(PCA$MEMORIA)
PCA$INTERACAO <- as.factor(PCA$INTERACAO)

PCA <- as.data.frame(PCA)

res.pca <- PCA(PCA[,c(3:29)], scale.unit = T, ncp = 3,
               quali.sup = c(1:5))

grupo <- cbind(PCA[,1:5])
grupo <- as.data.frame(grupo)

#The proportion of variances retained by the different dimensions
eig.val <- get_eigenvalue(res.pca)
fviz_eig(res.pca)
head(eig.val)

var <- get_pca_var(res.pca)
var
head(var$coord)
head(var$contrib)

head(var$cor)
head(var$cos2)


# Color variables by groups
fviz_pca_var(res.pca, col.var = "cos2", 
             axes = c(1,2),
             palette = c("ggplot2"),
             legend.title = "cos2",
             fill.var = "white",
)

# Create a grouping variable using kmeans
# Create 3 groups of variables (centers = 3)
set.seed(123)
res.km <- kmeans(var$coord, centers = 3, nstart = 1)
grp <- as.factor(res.km$cluster)


ind <- get_pca_ind(res.pca)
ind
head(ind$coord)
head(ind$contrib)
head(ind$cos2)

res.desc <- dimdesc(res.pca, axes = c(1:3), proba = 0.05)
# Description of dimension 1
res.desc$Dim.1
# Description of dimension 2
res.desc$Dim.2
# Description of dimension 3
res.desc$Dim.3

library("corrplot")
col4 <- colorRampPalette(c("red", "white","blue"))
corrplot(var$cor, is.corr=T, title = "Corelação das dimensões com os constructos", col = col4(150), tl.col = "black")
corrplot(res.desc$Dim.1$quanti, is.corr=T, title = "Correlation DIM. 1")
corrplot(res.desc$Dim.2$quanti, is.corr=T, title = "Correlation DIM. 2")


# Total cos2 of variables on Dim.1 until Dim.2
fviz_cos2(res.pca, choice = "var", axes = 1)
fviz_cos2(res.pca, choice = "var", axes = 2)

fviz_pca_var(res.pca)
# Contribution to the first dimension
fviz_contrib(res.pca, "var", axes = 1)
# Contribution to the second dimension
fviz_contrib(res.pca, "var", axes = 2)
# Contribution to the third dimension
fviz_contrib(res.pca, "var", axes = 3)



############################## PLOT - PCA #################
# Create the PCA biplot
fviz_pca_biplot(res.pca, 
                axes = c(1,2),
                # Individuals
                geom.ind = "point",
                fill.ind = PCA$GRUPO, col.ind = 3, col.ind.sup = "black",
                pointshape = 21, pointsize = 2,
                palette = "ggplot2",#c("darkred","red","orange","yellow"),
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


