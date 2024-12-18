library("RPostgreSQL")
library(DBI)
library(tidyr)
library(tidyverse)
library(ComplexHeatmap)
library(colorRamp2)


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



query<-"select bfstruct,cortexstruct, cellcount/injectioncount  cellcount_normalized  from  bf_region_larger11_stat where  cortexstruct  not like '%uniden%' ";

print(query)

df <- dbGetQuery(connec, query)


distanceMat<-df %>% pivot_wider(names_from = 'cortexstruct', values_from = 'cellcount_normalized')


distanceMat[is.na(distanceMat)] <- 0

distanceMat %>% column_to_rownames(., var = "bfstruct")

rnames <- distanceMat$bfstruct

distanceMat$bfstruct <- NULL

rownames(distanceMat) <-rnames


#distanceMat <- distanceMat[,-1]



#orignames<-names(df)
#preClusteringColumns<-orignames[!orignames %in% c('iddelineation2d')]

#print(df[preClusteringColumns])

#cortable<-cor(df[preClusteringColumns])




col_fun = colorRamp2(c(0,60), c("white", "red"))

Heatmap(data.matrix(distanceMat), 
        name = "mat",
        column_dend_height = unit(4, "cm"),
        row_dend_width = unit(4, "cm"),
        row_names_gp = gpar(fontsize =8),
        column_names_gp = gpar(fontsize =8),
        #clustering_distance_rows = "pearson",
        #clustering_method_rows = "single",
        #clustering_method_columns = "single",
        col_fun,
        column_title = "Cortex projection cell numbers in BF regions normalized by injection number",
        column_title_gp = gpar(fontsize = 14, fontface = "bold"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.0f", distanceMat[i, j]), x, y, gp = gpar(fontsize =5))
        }
        
)
