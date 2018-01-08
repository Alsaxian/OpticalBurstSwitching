#' ---
#' title: Attaque réseaux BHP
#' author: Jean-Baptiste AUJOGUE <jean-baptiste.aujogue@etu.univ-lyon1.fr>, Xian YANG <xian.yang@etu.univ-lyon1.fr>
#' date: 20 déc. 2017
#' output:
#'    html_document:
#'      theme: spacelab
#'      toc: true
#'      highlight: kate
#'      
#' 
#' ---
#+ results = "hide", echo = FALSE
set.seed(11715490)

#' &ensp; &ensp;   
#' &ensp; &ensp; 
#' 
#' ## 0. Introduction
#' Nous allons traiter un jeu de données concernant la sécurité du réseau. Il s'agit 
#' apparamment, lors d'un certain genre d'attaque du réseau, des états du réseau, de la transmission des paquets informatiques
#' etc. Le jeu de données est téléchargeable 
#' [sous ce lien](https://archive.ics.uci.edu/ml/machine-learning-databases/00404/OBS-Network-DataSet_2_Aug27.arff)
#' et la description des variables se trouve 
#' [ici](https://archive.ics.uci.edu/ml/datasets/Burst+Header+Packet+%28BHP%29+flooding+attack+on+Optical+Burst+Switching+%28OBS%29+Network#).
#' Notre but principal dans cet article est de faire une analyse en composantes principales sur le jeu de données. Néanmoins,
#' on verra que les données sont assez brutes qui nécessiteront un pré-traitement avant de faire l'ACP. 
#' 
#' L'article est organisé come suit :   
#' Dans les parties I et II on importe les données et traite les valeurs manquantes. Dans la partie III on étudie la 
#' colinéarité deux à deux des variables et discute la garde ou non des couples de variables extrèmement corrélés. On
#' argumente pourquoi dans notre analyse nous optons pour l'enlève des variables redondantes. La partie IV est une 
#' analyse complète des données en composantes principales, qui étudie et illustre entre autres la variance expliquée
#' par les dimensions principales,
#' les relations entres les variables numériques, la qualité de représentation et la contribution des variables par
#' rapport aux dimensions principales. Cette partie finit par 
#' une analyse factorielle des données mixtes sans détails mathématiques
#' pour révéler les relations entre les variables quantitatives et les variables qualitatives.
#' 
#' 
#' 
#' 
#' &ensp; &ensp;   
#' &ensp; &ensp; 
#' 
#' ## I. Importation du jeu de données et traitement des valeurs manquantes
#' Le package `foreign` nous permet de lire le jeu de données au format `.arff`. 
#+ results = "hide", warning = FALSE
library(foreign)
#+
setwd("C:/Users/ellie/Documents/R/OpticalBurstSwitching")
OBS <- read.arff("OBS-Network-DataSet_2_Aug27.arff")
# head(OBS) 
#' Il nous convient de connaître la taille du jeu de données et d'en faire une statistique descriptive.
dim(OBS)
summary(OBS) 
#' On voit que la variable `Packet Size_Byte` est complètement inutile car elle ne change jamais (=1440). On va la virer.
OBS <- OBS[names(OBS) != "Packet Size_Byte"]
#' Pour rendre les noms de variables faciles à traiter par `R`, on va remplacer des espaces par des points.
names(OBS) <- make.names(names(OBS), unique = TRUE)
#' Maintenant le jeu de données est prêt pour `R`.
#' 
#' &ensp; &ensp;   
#' &ensp; &ensp; 
#' 
#' ## II. Traitement des valeurs manquantes
#' On demande d'abord un aperçu des valeurs manquantes dans le jeu de données pour décider quel traitement prendre par 
#' la suite.
valMan <- which(is.na(OBS), arr.ind = TRUE, useNames = TRUE)
dim(valMan) 
#' Il s'avère que dans ce jeu de données de taille `r dim(OBS)[1]` * `r dim(OBS)[2]``,  il n'y a que `r dim(valMan)[1]` 
#' valeurs manquantes. On n'a donc pas besoin d'en faire une visualisation et le traitement sera relativement libre et simple.
#' La fonction suivante peut illustrer à l'aide d'un tableau la structure des valeurs manquantes.
#+ warning = FALSE
library(mice)
mdp <- md.pattern(OBS)
mdp
#' On voit que toutes les 15 valeurs manquantes appartiennent à la variable `Packet_lost`. Comme elles sont peu par rapport
#' au nombre d'observation, on pourrait bien entendu les ignorer toutes. Mais ici, on applique une technique d'imputation
#' qui porte le nom k plus proches voisins, remplaçant la valeur manquante par la moyenne des observations voisinée. Le nombre
#' des voisins sera fixé avec le reste des données.
#+ echo = FALSE, eval = FALSE
library(ipred)
# preproc <- preProcess(OBS, method = "bagImpute")

# mice(OBS, m=5, maxit=50, meth='pmm', seed=500)
# OBScomplet <- mice(OBS)
# aggr(mdp, prop = FALSE, numbers = TRUE) 
#+ warning = FALSE
library(DMwR)
OBScomplet <- knnImputation(OBS[! names(OBS) %in% c("Node.Status", "Class")])
# Quelles sont les valeurs imputées ?
OBScomplet[valMan]
#' 
#+ echo = FALSE, eval = FALSE
# install.packages("Amelia")
# library(Amelia)
# missmap(OBS)
#'
#' &ensp; &ensp; 
#' &ensp; &ensp;
#' 
#' ## III. Analyse des corrélations deux à deux 
#' Dans une analyse de données brutes il peut être intéressant de regarder d'abord s'il y a des corrélations élevées 
#' deux à deux parmi les variables numériques, afin de détecter celles qui sont susceptibles d'être redondantes. Une 
#' carte thermique (heatmmap) peut faciliter la visualisation.
#+ results="hide", warning = FALSE
library(plotly)
#+ 
corOBS <- cor(OBScomplet[! names(OBScomplet) %in% c("Node.Status", "Class")])
corHmp <- plot_ly(x = names(OBScomplet[! names(OBScomplet) %in% c("Node.Status", "Class")]), 
                  y = names(OBScomplet[! names(OBScomplet) %in% c("Node.Status", "Class")]), z = corOBS, type = "heatmap")
corHmp
#' In the correlation heatmap we can repeatly spot extremely high pairwise correlation close or even equal to 1 or -1.  
#'    
#' Highly to perfectly correlated variables could bias the PCA result in a way that PCA will overemphasize the common contribution 
#' of the (nearly) redundant variables. Therefore, it might make sense to find out and remove such variables before doing a PCA,
#' __Especially when they describe actually (nearly) the same aspect of an issue.__ 
#' For instance, `Tansmitted_Byte` and `Packet_Transmitted` have a correlation of 1, because they reflect the same quantity 
#' up to a ratio `Packet Size_Byte`, which, as being mentioned here above, never changes. Thus, we will remove one of both.
#' Which one to remove is of our free choice. Here, we consider that quantitative variables in packets may be more reader 
#' friendly than those in bytes in terms of unit, so that we keep the variable `Packet_Transmitted` and drop the other one.
#' But what about the high correlation between other variables (e.g. higher than 97% but not perfectly equal to 1) ? 
#' Shall we remove the nearly redundant variables too before doing our PCA?   
#'   
#' The statistical community doesn't have a 
#' straightforward anwer to it. As a matter of fact, it hugely depends on the nature of data and the purpose to do the PCA.
#' On one hand, like we said, highly correlated variables would be possible to strongly influence the result of the PCA
#' and, as a result, the real contributions to the principal components of the underlying variables that are truely meaningful.
#' If our PCA is meant to give such information, then high correlation should better be avoided prior to the PCA.
#' On the other hand, however, a PCA with redundant variables can still faithfully reveal the high correlation between them, though
#' principal components would be probably established otherwise. In that sense, if it is an exploratory PCA that we are doing,
#' which only aims to find a broad outline of the relationships between variables disregarding how principal components are built,
#' then including some redundant ones may be fine.  
#'   
#' In the light of this, in the next section 
#' we decide to go an onerous but careful way, in which we do firstly a PCA with almost all the variables.
#' With both the correlation circle of PCA and the correlation heatmap, we then kick out the redundant variables and do a second
#' PCA with only variables that we consider to be meaningful. Further analysis (variable and individual relationship, 
#' quantitative and qualitative variable relationship) shall also be based on the latter PCA.  
#'   
#' More technically, 
#' in a rough and rapid preprocessing procedure, we could
#' limit our focus on high pairwise correlations (the reader should know that in a more rigorous treatment,
#' high correlations within a multuple tuple of variables should also be considered) and
#' refer to the heatmap of variable correlations. We may set a threshold, e.g. 0.97, and find out all couples whose absolute
#' value of correlation exceeds this threshold, before we remove one of both variables by verifying that they are indeed telling (nearly) the same story.
#' 
#' &ensp; &ensp;   
#' &ensp; &ensp; 
#' 
#' ## IV. Analyse en composantes principales
#' ### 1. Analyse des variables numériques
#+warning = FALSE
library("factoextra")
pcaOBS <- prcomp(OBScomplet[! names(OBScomplet) %in% c("Node.Status", "Class")], scale. = TRUE, rank. = 3)
fviz_eig(pcaOBS, ncp = 3, addlabels = TRUE, ylim = c(0, 60))
summary(pcaOBS)$importance[3, c(1, 2)]
fviz_pca_var(pcaOBS, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, title =
               "PCA of all numeric variables")
#' We see that in the correlation circle there are four subgroups of very closely situated variables which are potentially
#' redundant variables (visibly there are five, but the one at the very left is highly negatively correlated with the one at the
#' very right). For example, `Tansmitted_Byte`, `Packet_Transmitted` and `Full_Bandwidth` are even perfectly correlated,
#' which can be also confirmed by the previous correlation heatmap. No doubt those subgroups consist of the _most contributing
#' variables (the longest arrows that touch almost the circle)
#' _, because they are redundant! This may be due to the fact that raw data are collected by non-statisticians 
#' and a primary treatment of unusual variables is missing. Including them in the PCA might not be of interest.   
#'    
#' We now decide to remove extremely high correlations among variables (greater than
#' 97%) for our next PCA so that only
#' one variable of each of the four subgroups should stay in the game. Another reason for doing this, in a point of view
#' of field knowledge, is that all the variables
#' within the same subgroup mean in fact the same thing just in some different way. At the end, we choose to keep 
#' `Packet.Received..Rate`, `Packet_Received`, `Packet_Transmitted` and `Packet_Lost` as representatives of their belonging subgroup along with 
#' other well defined variables. In the more balanced coming PCA with uniquely the non redundant variables, we could spot
#' a change in variance contribution of the variables, as well as a diminuation of explained variance by the principal 
#' components because data are now distributed in a more balanced way.
#+ warning = FALSE
pcaOBS2 <- prcomp(OBScomplet[c("Packet.Received..Rate", "Packet_Received", "Packet_Transmitted", "Packet_lost",
                             "Average_Delay_Time_Per_Sec", "X10.Run.Delay", "Node", "Received_Byte", "Flood.Status")], scale. = TRUE, rank. = 3)
summary(pcaOBS2)$importance[3, c(1, 2, 3)]
fviz_eig(pcaOBS2, ncp = 3, addlabels = TRUE, ylim = c(0, 60))
fviz_pca_var(pcaOBS2, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, title =
               "PCA of meaningful numeric variables")
#' Following information can be extracted from the correlation circle of the PCA graph above:   
#'   
#'   `Flood.Status`, `Packet_lost`, `Packet_Transmitted`, `Packet_Received` and `Packet.Received..Rate` are the five most well
#' represented or in other words most contributing variables to the first both principal dimensions.
#'   `Packet_Received` has no contribution to the first dimension while `Node` has no contribution to the second one.
#'   `Received_Byte` has little contribution to the first dimension while `Packe_lost` has litte contribution to the second one.
#'   `X10.Run.Delay` and Àverage_Delay_Time_Per_Sec` seem to be correlated in the projected dataset onto the first two dimensions.
#'   `Packet.Received..Rate` and `Flood.Status` are highly negatively correlated.  
#'  
#' To extract all information about the remaining variables, we can do
#+ warning = FALSE
variables <- get_pca_var(pcaOBS2)
#' What are the scores of them?
#+ warning = FALSE
variables$coord[, 1:3]
#' We may want to illustrate them in a `corrplot`:
#+ warning = FALSE
variables$cos2[, 1:2]
library("corrplot")
corrplot(variables$cos2[, 1:2], is.corr=FALSE, title = "Quality of representation of the variables by the first two components")
#' Or sum it up: 
#+ warning = FALSE
fviz_cos2(pcaOBS2, choice = "var", axes = 1:2, title = "Total representation quality on dimension 1 and 2")
#' How are the contributions of the variables to the first two components?
#+ warning = FALSE
variables$contrib[, 1:2]
corrplot(variables$contrib[, 1:2], is.corr=FALSE, title = "Contributions of the variables to the first two components")
#' Or let's present them ranked:
#+ warning = FALSE
fviz_contrib(pcaOBS2, choice = "var", axes = 1, title = "Contribution of variables to PC1")
fviz_contrib(pcaOBS2, choice = "var", axes = 2, title = "Contribution of variables to PC2")

#' ### 2. Analyse des individus
#' Now we would like to plot some individuals. We can see that here, using redundant variables in the PCA or not will
#' lead to a considerable difference.  
#'   
#' Firstly let's see a graph of 20 randomly selected individuals drawn in the PCA including redundant variables, coloured by their quality of 
#' representation.
fviz_pca_ind(pcaOBS, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, 
             select.ind = list(name = sample.int(n = nrow(OBS), size = 20)))
#' Their quality of representation is always the best, equal to 1! This is due to the fact that the first dimension
#' is biased on the redundant information, so that an individual can be easily well represented by 
#' only these variables because they repeatly occur.  
#'  
#' Then we show a graph of 20 randomly selected individuals drawn in the PCA _without_ redundant variables, coloured by their quality of 
#' representation
#+ warning = FALSE
fviz_pca_ind(pcaOBS2, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, 
             select.ind = list(name = sample.int(n = nrow(OBS), size = 20)))
#' All kinds of quality can be observed this time, which is more meaningful than in the previous case. Hence, this PCA pragh is what 
#' we are going to keep using in the rest of our analysis.  
#'   
#' We now try to include the categorical information: `Node.Status` and `Class` in our visualisation.
#' Firstly a graph of some randomly selected individuals coloured by subgroups of `Node.Status`
#+ warning = FALSE
fviz_pca_ind(pcaOBS2, geom = "point", habillage = OBS$Node.Status,  palette = "jco", 
             addEllipses = FALSE, legend.title = "Node Status")
#' Then a graph of randomly selected individuals coloured by subgroups of `Node.Status`
#+ warning = FALSE
fviz_pca_ind(pcaOBS2, geom = "point", habillage = OBS$Class,  palette = "RdBu",
            select.ind = list(name = 1:1075), legend.title = "Class")
#' ### 3. Présentation graphique intégrée des variables et des individus            
#' Now we draw some biplots with randomly selected individuals and variables
#+ warning = FALSE
library(gridExtra)
plot1 <- fviz_pca_biplot(pcaOBS2, 
                habillage = OBS$Node.Status, palette = "jco", 
                addEllipses = FALSE, label = "var",
                col.var = "darkgreen", repel = TRUE, title = NULL,
                legend.title = "Node Status", 
                subtitle = "Coloured by Node Status", xlab = "PC1", ylab = "PC2") 

plot2 <- fviz_pca_biplot(pcaOBS2, 
                         habillage = OBS$Class, palette = "RdBu", 
                         addEllipses = FALSE, label = "var",
                         col.var = "darkgreen", repel = TRUE,
                         legend.title = "Class", title = NULL,
                         subtitle = "Coloured by Class", xlab = "PC1", ylab = "PC2") 
grid.arrange(plot1, plot2, ncol=2, top = "Randomly selected individuals with meaningful variables")
#' ### 4. Analyse Factorielle des Données Mixtes
#' Finalement, on veut faire l'AFDM pour résumer les ralations entre les variables quantitatives et les variables
#' qualitatives, dans un seul graphe de façon
#' plus intuitive. On vient de visualiser ces relations dans le paragraphe précédent mais à travers la coloration des individus.
#+ warning = FALSE 
library(FactoMineR)
#+ echo = FALSE, results = "hide"
#### Attention: On ne peut pas faire le graphe de quanti et quali ensemble! Car elles ne partagent pas les mêmes coordonnees.
#### Par contre quali et individus ensemble c'est possible mais ca revient à sous-groupes avec ellipse.
#### Pour tracer quali et quanti ensemble, il faut absolument passer par FAMD.
#+ warning = FALSE
pcaOBS3 <- PCA(cbind(OBScomplet[c("Packet.Received..Rate", "Packet_Received", "Packet_Transmitted", "Packet_lost",
                         "Average_Delay_Time_Per_Sec", "X10.Run.Delay", "Node", "Received_Byte", "Flood.Status")],
                         OBS[c("Node.Status", "Class")]), quali.sup = 10:11, graph = FALSE)
p <- fviz_pca_var(pcaOBS3)
fviz_add(p, pcaOBS3$quali.sup$coord, color = "red")


famdOBS <- FAMD(cbind(OBScomplet[c("Packet.Received..Rate", "Packet_Received", "Packet_Transmitted", "Packet_lost",
                              "Average_Delay_Time_Per_Sec", "X10.Run.Delay", "Node", "Received_Byte", "Flood.Status")],
                 OBS[c("Node.Status", "Class")]), graph = FALSE)

cumsum(famdOBS$eig[, 2])
plot(famdOBS, choix = "var")

