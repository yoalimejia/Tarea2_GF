######EJERCICIO 5######
#Selecciona alg ??un experimento que te parezca interesante de micro-arreglosde expresi ??o
#v??ia la herramienta Geo2R de NCBI. Realiza un an ??alisis de exprei ??on diferencial ent
#re al menos dos grupos experimentales que se reporten en el art ??iculo correspondiente.
#Genera el c ??odigo para realizar todos los analisis como las gr ??aficas y almac??enalo en tu cuenta de Github.
#Comenta el codigo de cade secci??on para que me quede claro que entiendes cada paso par ainferir la red.
#Manda la liga de ese c ??odigo. Genera las tablas en formato csv de genes diferencialmente expresados con las siguientes
#condiciones: logFC de al menos 2 y p-value inferior a 0.1.


# Análisis de expresión diferencial con limma
#   Differential expression analysis with limma
#primero se cargaron las libreríias a utilizar 
library(GEOquery) #se instala de BioConductor
library(limma) #previamente instalada 
library(umap) #esta la instalé a parte en R de CRAN
BiocManager :: install ("GEOquery")
# load series and platform data from GEO
#se deben de cargar la serie del GEO con el que se buscó en GEO2R que comiena con GSE
gset <- getGEO("GSE18388", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6246", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]] #se establecen los parámetros con ayuda de un if else, sus condicionales. 

fvarLabels(gset) <- make.names(fvarLabels(gset)) #se generan columnas para asignar los datos 

gsms <- "11110000" #grupo al que pertenecen cad auno de los miembros de genes 
sml <- strsplit(gsms, split="")[[1]]

#se lleva a cabo un logFC 
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

#se agregan ya genes o experimentos a cada grupo realidado, en este caso, control y muted 
gs <- factor(sml)
groups <- make.names(c("control","muted"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # establecer el modelo lineal de los grupos 

#Aquí se seleccionan los grupos de interés para el contraste de los modelos
#Se calculan los coeficientes con base en esa selección 
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

#modelo estadístico con el modelo de Bayes y se establecen parámetros específicos 
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
#se fijan en una tabla con difernetes características para filtrar los datos 
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

#Aquí se generó una tabla de calidad con base en el control y se realiza un histograma para visuallizar los datos
#el histograma se realiza con base en el p valor de 0.1 para ver los genes diferencialmente expresados.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")
#observar genes que están sobre regulados, subregulados y no expresados con base en el p valor 
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.1)

#Diagrama de venn vara observar resultados obtenidos hasta ahora 
vennDiagram(dT, circle.col=palette())

# QQ plot estadístico, se quitan todos aquellos grupos que no cumplen con los parámetros--> T
t.good <- which(!is.na(fit2$F))
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (P valor y logFC)
colnames(fit2) #lista de datos
ct <- 1        # grupo contraste 
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2))) #parámetros para el volcano plot 

# MD plot (logFC vs mean log expression)
#Se destaca estadísticamente a aquellos que tienen un pvalor <-0.1
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
#Análisis de datos, se lo aisganmos a un vector
ex <- exprs(gset)

# box-and-whisker plot
ord <- order(gs)  #se ordenan la smuestras con abse en los grupos establecidos 
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE18388", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topleft", groups, fill=palette(), bty="n")

# valor de la distribución d elos grupos 
par(mar=c(4,4,2,1))
title <- paste ("GSE18388", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")


ex <- na.omit(ex) # 
ex <- ex[!duplicated(ex), ]  # 
ump <- umap(t(ex), n_neighbors = 4, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=4", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)
library("maptools")  # cargar librería 
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

#promedio de las varianzas para ver precisión de los pesos de los grupos.
plotSA(fit2, main="Mean variance trend, GSE18388")


