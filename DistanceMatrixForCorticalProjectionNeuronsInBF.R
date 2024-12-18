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



#query<-"select * from average_cell_distance_cell_populations where name1 not like '%uniden%' and name2 not like '%uniden%'  ";

query <-"select name1, name2, median from stat_matview_median_distance_bf_cell_populations where name1 not like '%uniden%' and name2 not like '%uniden%' "

print(query)

df <- dbGetQuery(connec, query)

#df <- round(df, digits = 0)

#distanceMat<-df %>% pivot_wider(names_from = 'name2', values_from = 'average_distance')

distanceMat<-df %>% pivot_wider(names_from = 'name2', values_from = 'median')

#set first column as row names

distanceMat %>% column_to_rownames(., var = "name1")

rnames <- distanceMat$name1

distanceMat$name1 <- NULL

rownames(distanceMat) <-rnames




#distanceMat <- distanceMat[,-1]



#orignames<-names(df)
#preClusteringColumns<-orignames[!orignames %in% c('iddelineation2d')]

#print(df[preClusteringColumns])

#cortable<-cor(df[preClusteringColumns])

jet.colors <-
  colorRamp2(c(0,400,800,1200,1600,2000,2400,2800,3200),c("#7F0000",
                                                          "red",
                                                          "#FF7F00", 
                                                          "yellow",
                                                          "#7FFF7F", 
                                                          "cyan",
                                                          "#007FFF", 
                                                          "blue", 
                                                          "#00007F"
  )
  )




col_fun = jet.colors; #jet.colors(5) #colorRamp2(c(0,3000), c("red","white"))

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
        column_title = "Median Distance between Cortical Projection Populations in BF",
        column_title_gp = gpar(fontsize = 14, fontface = "bold"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.0f", distanceMat[i, j]), x, y, gp = gpar(fontsize =5))
        }
        
)
