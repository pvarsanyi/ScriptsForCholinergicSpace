
library("RPostgreSQL")
library(DBI)
library(tidyr)
library(ComplexHeatmap)
library(gplots)
library(colorRamp2)
library(pvclust)
library(reshape)
library(cluster)
library(ggplot2)
library(ggrepel)
library(NMF)
library(factoextra) 
set.seed(123)



idvectorspace <-99#103;#99 #87 #94 #87;#94#87;#92;#90#89;#70;#67;#87;


tryCatch({
  drv <- dbDriver("PostgreSQL")
  print("Connecting to Database…")
  connec <- dbConnect(drv, 
                      dbname = 'anatomy',
                      host = 'localhost', 
                      user = 'anatomy', 
                      password = 'anatomy')
  
  print("Database Connected!")
},
error=function(cond) {
  print("Unable to connect to Database.")
  
})

#(case when parent is not null then parent||'/'|| nm else nm end) as
#order by  nmiddelineation2d,namedstructure.structureorder
query<-paste("select iddelineation2d , nm as name,value from (select iddelineation2d, (select  larger1  from ontology where Larger11  =namedstructure.name and namedstructure.name <> larger1 limit 1) parent , namedstructure.name nm, value from vectorcomponent, namedstructure where namedstructure.idnamedstructure  = vectorcomponent.idnamedstructure  and namedstructure.name not like '%unid%' and idvectorspace =",idvectorspace," order by  iddelineation2d,namedstructure.structureorder ) st")               

#order by original order 
ord <- " order by (case when parent is not null then parent||'/'|| nm else nm end)";

query<-paste(query, ord)
               
print(query)

df <- dbGetQuery(connec, query)


df<-df %>% pivot_wider(names_from = 'name', values_from = 'value')

orignames<-names(df)
preClusteringColumns<-orignames[!orignames %in% c('iddelineation2d')]

#print(df[preClusteringColumns])

cortable<-cor(df[preClusteringColumns])



col_fun = colorRamp2(c(-1, 0, 1), c("blue","white","red"))



hm <- Heatmap(cortable, 
        name = "mat",
        column_dend_height = unit(4, "cm"),
        row_dend_width = unit(4, "cm"),
        row_names_gp = gpar(fontsize =11),
        column_names_gp = gpar(fontsize =11),
        clustering_method_rows = "complete",
        clustering_method_columns = "complete",
        clustering_distance_columns = "euclidean",
        clustering_distance_rows = "euclidean",
        cluster_columns = TRUE, 
        cluster_rows = TRUE, 
        row_names_side = c("right"),
        column_dend_side = c("top"),
        row_dend_side = c( "right"),
        column_names_side = c("top"),
        show_column_dend =FALSE,
        show_row_dend = TRUE,
        #row_dend_reorder =FALSE,
        #column_dend_reorder =FALSE,
 
        
        #clustering_distance_rows = "pearson",
        #clustering_method_rows = "single",
        #clustering_method_columns = "single",
       # column_title = "Spacial Correlation of BF neurons projecting to Cortical Areas",
         column_title_gp = gpar(fontsize = 14, fontface = "bold"),
        
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", cortable[i, j]), x, y, gp = gpar(fontsize =5))
        }
        )
hm

cophenetic <- cophcor(cortable, linkage =  "complete")

#par(mar=c(11,15,4,1))
par(cex.axis=0.7)

df_ct <- data.frame(cortable)

rownames(df_ct) <- c("A1","AAF","PDAF","SRAF","AIP","Amygdala","LEC","MEC","Hipp","IDIGI","M1Wh","M1FL","M1HL","M2","Cg1","IL","PrL","LO","VO","PER","POR","RS","S1Or","S1Wh","S1FL","S1HL","S1Tr","S2","V1","V2")
colnames(df_ct) <-rownames(df_ct)


dis = dist(df_ct);

frame <- data.frame(dis);

cl <- hclust(dist(cortable), "complete")

cut = cutree(cl, k=18)

print (cut)

sil_cl <- silhouette(cut ,dis)

#sil_df <- data.frame(sil_cl )

print(sil_cl)

rownames(sil_cl) <- rownames(df_ct )

#rownames(sil_cl) <- rownames(cortable)

#rownames(sil_cl) <- c("A1","AAF","PDAF","SRAF","AIP","Amygd"  ,"LEC", "MEC", "Hipp","IDIGI","M1Wh","M1FL","M1HL","M2","Cg1","IL","PrL","LO","VO","PER","POR","RS", "S1OR","S1WH","S1FL","S1HL" ,"S1Tr","S2",      "V1",      "V2")     ;


plot(sil_cl)


#reordered_corr_matrix <- data.frame(hm@matrix)
#print(hm@matrix)
#hm@row_dend_reorder

heatmap(cortable, 
          name = "mat",
          column_dend_height = unit(4, "cm"),
          row_dend_width = unit(4, "cm"),
          row_names_gp = gpar(fontsize =8),
          column_names_gp = gpar(fontsize =8),
          hclustfun = function (x)hclust(x, method="complete"),
            #clustering_distance_rows = "pearson",
            #clustering_method_rows = "single",
            #clustering_method_columns = "single",
            column_title = "Spacial Correlation of BF neurons projecting to Cortical Areas",
          column_title_gp = gpar(fontsize = 14, fontface = "bold"),
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", cortable[i, j]), x, y, gp = gpar(fontsize =5))
          }
)




############################################### PCA of correlation matrix ####################################################

pca <- prcomp(cortable, scale = TRUE)

print(pca)

plot(pca$x[,1],pca$x[,2])

pca.var <- pca$sdev^2
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)

barplot(pca.var.per, main="Percentage Variances of the Principal Components", xlab="Principal Component", ylab="Percent Variation")

pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2])
pca.data

ggplot(data=pca.data, aes(x=X, y=Y, label=Sample)) +

  xlab(paste("PC1 - ", pca.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pca.var.per[2], "%", sep="")) +

  ggtitle("")+
  geom_label_repel(aes(label = Sample),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme_classic()+
  geom_point(color = "blue", size = 3)



############################################### PCA of original data ####################################################

origdata<-t(df[preClusteringColumns])

pcaorig <- prcomp(origdata, scale = TRUE)

print(pcaorig )

plot(pcaorig$x[,1],pcaorig$x[,2])

pcaorig.var <- pcaorig$sdev^2
pcaorig.var.per <- round(pcaorig.var/sum(pcaorig.var)*100, 1)

barplot(pcaorig.var.per, main="Scree Plot", xlab="Principal Component", ylab="Percent Variation")

pcaorig.data <- data.frame(Sample=rownames(origdata),
                       X=pcaorig$x[,1],
                       Y=pcaorig$x[,2])

pcaorig.data
ggplot(data=pcaorig.data, aes(x=X, y=Y, label=Sample)) +
  geom_text() +
  xlab(paste("PC1 - ", pcaorig.var.per[1], "%", sep="")) +
  ylab(paste("PC2 - ", pcaorig.var.per[2], "%", sep="")) +
  theme_bw() +
  ggtitle("My PCA Graph")

##############################################  Correlation between dendogram distance andd real distance #################################


hc <- hclust(dist(cortable), "complete") 



query <- " SELECT name1 as X1 , name2 as X2, median as value FROM public.stat_matview_median_distance_bf_cell_populations where not (name1  like '%unide%' or name2  like '%unide%')"


phisical_distance <- dbGetQuery(connec, query)

colnames(phisical_distance ) <- c("X1","X2","value");

phisical_distance <- subset(phisical_distance , X1!=X2);

print(phisical_distance["X1"])

phisical_distance <-phisical_distance[order(phisical_distance$X1,phisical_distance$X2),]
rownames(phisical_distance) <- 1:nrow(phisical_distance)

cophenetic_dist <- cophenetic(hc) 

cophenetic <- data.frame(melt(as.matrix(cophenetic_dist))); #

cophenetic  <- subset(cophenetic  , X1!=X2);
rownames(cophenetic ) <- 1:nrow(cophenetic )


cophenetic <-cophenetic[order(as.character(cophenetic$X1),as.character(cophenetic$X2)),]



cor(cophenetic$value, phisical_distance$value)



plot(cophenetic$value, phisical_distance$value, pch = 16, cex = 1.3, col = "blue", main = "", xlab = "Cophenetic Distance", ylab = "Median cell distance between projection cell populations")

abline(lm( phisical_distance$value ~ cophenetic$value ),col="black")

############################################# Kmeans reflect to review by HANGYA####################################################



km.res <- kmeans(cortable, 3, nstart = 25)

print(km.res["betweenss"])

fviz_cluster(km.res, data = cortable,labelsize = 12)


fviz_nbclust(cortable, kmeans, method = "wss")

fviz_nbclust(cortable, kmeans, method = "silhouette")


for (k in 1: 20)
{
  
  km.res <- kmeans(cortable, k, nstart = 25)
  print(cbind(k,':',km.res["tot.withinss"], '   ', km.res["withinss"]))
}


##############################################  Correlation between correlation vector and real distance #################################


corr_vector_distance <- data.frame(melt(as.matrix(dist(cortable))));

corr_vector_distance <- subset(corr_vector_distance , X1!=X2);

corr_vector_distance <-corr_vector_distance[order(as.character(corr_vector_distance$X1),as.character(corr_vector_distance$X2)),]

rownames(corr_vector_distance) <- 1:nrow(corr_vector_distance)


cor(corr_vector_distance$value, phisical_distance$value)

cor.test(corr_vector_distance$value,  phisical_distance$value)$p.value


plot(corr_vector_distance$value, phisical_distance$value, pch = 16, cex = 1.3, col = "blue", main = "", xlab = "Correlation vector distance (Euclidean)", ylab = "Median cell distance between projection cell populations")

abline(lm( phisical_distance$value ~ corr_vector_distance$value ),col="black")



##############################################  Correlation between correlation vector distance  and dendogram distance  #################################


cor(corr_vector_distance$value, cophenetic$value)

plot(corr_vector_distance$value, cophenetic$value, pch = 16, cex = 1.3, col = "blue", main = "", xlab = "Correlation vector distance (Euclidean)", ylab = "Cophenetic Distance")

abline(lm( cophenetic$value ~ corr_vector_distance$value ),col="black")

write.csv(corr_vector_distance, file='corr_vector_distance.csv')

write.csv(cophenetic, file='cophenetic.csv')


############################################## Correlation between vector distance and injection site distance ##########################################



query <- "select name1 as X1 , name2 as X2, avg as value from injection_site_centroids_distance_larger11 where name1 not like '%unid%' and name2 not like '%unid%' and name1 <> name2 order by name1, name2"



injection_site_distance <- dbGetQuery(connec, query)


colnames(injection_site_distance) <- c("X1","X2","value");


#(X1 =='A1' | X1 =='SRAF (A2)' | X1 =='PDAF (A2)' | X1 =='PDAF (A2)' | X1 =='AAF (A2)' | X1 =='IDIGI' | X1 =='Amygdala' | X1 =='PER' | X1 =='POR' | X1 =='AIP') &
#(X2 =='A1' | X2 =='SRAF (A2)' | X2 =='PDAF (A2)' | X2 =='PDAF (A2)' | X2 =='AAF (A2)' | X2 =='IDIGI' | X2 =='Amygdala' | X2 =='PER' | X2 =='POR' | X2 =='AIP')


#(X1 =='V1' | X1 =='Hippocampus' | X1 =='V2' | X1 =='Cg1' | X1 =='IL' | X1 =='LEC' | X1 =='MEC' | X1 =='RS' | X1 =='PrL' ) &
#(X2 =='V1' | X2 =='Hippocampus' | X2 =='V2' | X2 =='Cg1' | X2 =='IL' | X2 =='LEC' | X2 =='MEC' | X2 =='RS' | X2 =='PrL' )


#(X1 =='LO' | X1 =='S1Tr (trunk) lower trunk' | X1 =='M2' | X1 =='VO' | X1 =='M1 (Agl,PrCl)hindlimb' | X1 =='S1-whisker' | X1 =='S1HL(hindlimb)' | X1 =='S1-orofacial' | X1 =='M1 (Agl,PrCl) whisker' | X1=='M1 (Agl,PrCl)forelimb' | X1=='S1FL (forelimb)'  ) &
#(X2 =='LO' | X2 =='S1Tr (trunk) lower trunk' | X2 =='M2' | X2 =='VO' | X2 =='M1 (Agl,PrCl)hindlimb' | X2 =='S1-whisker' | X2 =='S1HL(hindlimb)' | X2 =='S1-orofacial' | X2 =='M1 (Agl,PrCl) whisker' | X2=='M1 (Agl,PrCl)forelimb' | X2=='S1FL (forelimb)'  ) 


corr_vector_distance_filtered <- subset(corr_vector_distance, ( X1 != 'S2' & X2 !='S2' & 
                                                                  (X1 =='LO' | X1 =='S1Tr (trunk) lower trunk' | X1 =='M2' | X1 =='VO' | X1 =='M1 (Agl,PrCl)hindlimb' | X1 =='S1-whisker' | X1 =='S1HL(hindlimb)' | X1 =='S1-orofacial' | X1 =='M1 (Agl,PrCl) whisker' | X1=='M1 (Agl,PrCl)forelimb' | X1=='S1FL (forelimb)'  ) &
                                                                  (X2 =='LO' | X2 =='S1Tr (trunk) lower trunk' | X2 =='M2' | X2 =='VO' | X2 =='M1 (Agl,PrCl)hindlimb' | X2 =='S1-whisker' | X2 =='S1HL(hindlimb)' | X2 =='S1-orofacial' | X2 =='M1 (Agl,PrCl) whisker' | X2=='M1 (Agl,PrCl)forelimb' | X2=='S1FL (forelimb)'  ) ))

                                                                
injection_site_distance_filtered <-  subset(injection_site_distance, ( X1 != 'S2' & X2 !='S2' & 
                                                                         
                                                                         (X1 =='LO' | X1 =='S1Tr (trunk) lower trunk' | X1 =='M2' | X1 =='VO' | X1 =='M1 (Agl,PrCl)hindlimb' | X1 =='S1-whisker' | X1 =='S1HL(hindlimb)' | X1 =='S1-orofacial' | X1 =='M1 (Agl,PrCl) whisker' | X1=='M1 (Agl,PrCl)forelimb' | X1=='S1FL (forelimb)'  ) &
                                                                         (X2 =='LO' | X2 =='S1Tr (trunk) lower trunk' | X2 =='M2' | X2 =='VO' | X2 =='M1 (Agl,PrCl)hindlimb' | X2 =='S1-whisker' | X2 =='S1HL(hindlimb)' | X2 =='S1-orofacial' | X2 =='M1 (Agl,PrCl) whisker' | X2=='M1 (Agl,PrCl)forelimb' | X2=='S1FL (forelimb)'  ) ))


correlation <- cor(corr_vector_distance_filtered$value,injection_site_distance_filtered$value)

pvalue <- cor.test(corr_vector_distance_filtered$value,injection_site_distance_filtered$value)

pvalue



plot(corr_vector_distance_filtered$value, injection_site_distance_filtered$value, pch = 16, cex = 1.3, col = "blue", main = "", xlab = "Correlation vector distance (Euclidean)", ylab = "Injection Site Distance")

abline(lm( injection_site_distance_filtered$value ~ corr_vector_distance_filtered$value ),col="black")



############################################## Finding  P value ####################################################


result <- pvclust(cortable, method.dist="euclidean", method.hclust="complete", nboot=1000, parallel=TRUE)          
          
 plot (result)         
 pvrect(result, alpha=0.90)
          
          
          
  ################################## Complex heatmap parameters (Heatmap function with a capital "H") ******************************************        
          
 na_col = "grey",
 color_space = "LAB",
 rect_gp = gpar(col = NA),
 border = NA,
 border_gp = gpar(col = "black"),
 cell_fun = NULL,
 layer_fun = NULL,
 jitter = FALSE,
 
 row_title = character(0),
 row_title_side = c("left", "right"),
 row_title_gp = gpar(fontsize = 13.2),
 row_title_rot = switch(row_title_side[1], "left" = 90, "right" = 270),
 column_title = character(0),
 column_title_side = c("top", "bottom"),
 column_title_gp = gpar(fontsize = 13.2),
 column_title_rot = 0,
 
 cluster_rows = TRUE,
 cluster_row_slices = TRUE,
 clustering_distance_rows = "euclidean",
 clustering_method_rows = "complete",
 row_dend_side = c("left", "right"),
 row_dend_width = unit(10, "mm"),
 show_row_dend = TRUE,
 row_dend_reorder = is.logical(cluster_rows) || is.function(cluster_rows),
 row_dend_gp = gpar(),
 cluster_columns = TRUE,
 cluster_column_slices = TRUE,
 clustering_distance_columns = "euclidean",
 clustering_method_columns = "complete",
 column_dend_side = c("top", "bottom"),
 column_dend_height = unit(10, "mm"),
 show_column_dend = TRUE,
 column_dend_gp = gpar(),
 column_dend_reorder = is.logical(cluster_columns) || is.function(cluster_columns),
 
 row_order = NULL,
 column_order = NULL,
 
 row_labels = rownames(matrix),
 row_names_side = c("right", "left"),
 show_row_names = TRUE,
 row_names_max_width = unit(6, "cm"),
 row_names_gp = gpar(fontsize = 12),
 row_names_rot = 0,
 row_names_centered = FALSE,
 column_labels = colnames(matrix),
 column_names_side = c("bottom", "top"),
 show_column_names = TRUE,
 column_names_max_height = unit(6, "cm"),
 column_names_gp = gpar(fontsize = 12),
 column_names_rot = 90,
 column_names_centered = FALSE,
 
 top_annotation = NULL,
 bottom_annotation = NULL,
 left_annotation = NULL,
 right_annotation = NULL,
 
 km = 1,
 split = NULL,
 row_km = km,
 row_km_repeats = 1,
 row_split = split,
 column_km = 1,
 column_km_repeats = 1,
 column_split = NULL,
 gap = unit(1, "mm"),
 row_gap = unit(1, "mm"),
 column_gap = unit(1, "mm"),
 show_parent_dend_line = ht_opt$show_parent_dend_line,
 
 heatmap_width = unit(1, "npc"),
 width = NULL,
 heatmap_height = unit(1, "npc"),
 height = NULL,
 
 show_heatmap_legend = TRUE,
 heatmap_legend_param = list(title = name),
 
 use_raster = NULL,
 raster_device = c("png", "jpeg", "tiff", "CairoPNG", "CairoJPEG", "CairoTIFF", "agg_png"),
 raster_quality = 1,
 raster_device_param = list(),
 raster_resize_mat = FALSE,
 raster_by_magick = requireNamespace("magick", quietly = TRUE),
 raster_magick_filter = NULL,
 
 post_fun = NULL)          
          
          
          
          
          
          
          
          