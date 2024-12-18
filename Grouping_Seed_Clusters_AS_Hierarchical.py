
import inline as inline
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import psycopg2
from pandas import DataFrame
import math

baseclustering =  173;
idvectorspace = 103;#58;

df = pd.read_excel('C:/Users/pvars/Dropbox/PVARSANYI/PROCESSING/2022-2023/Manuscript+statistics+figures/DRAFT/Tables-Excels-Suppl/HierarchicalClusterCompositionsnnNoUnidentified.xlsx');
conn = psycopg2.connect("dbname=anatomy user=anatomy password=anatomy");


cursor = conn.cursor()

cursor.execute(
    "INSERT INTO clustering  (name,idvectorspace,idvectorspacedisplay) values ('ClusteringNoUnidentified( Hierarchical  ) '," + str(
        idvectorspace) + "," + str(idvectorspace) + ") RETURNING idclustering")
idclustering = cursor.fetchone()[0]



for index, row in df.iterrows():

    cluster_label = row[0];

    cluster_in_query="";
    cluster_description="";

    for i in range(1,len(row)):

        if str(row[i]) != 'nan':
            cluster_in_query+="'"+str(row[i])+"'"
            cluster_description +=  str(row[i]) + ","


    cluster_in_query= "insert into clusterpoint ( idclustering,iddelineation2d,clusterlabel,description)  select distinct "+str(idclustering)+",iddelineation2d,"+str(cluster_label)+",'" + cluster_description[:60]+"' from clusterpoint where idclustering = "+str(baseclustering) +" and description in ("+cluster_in_query.replace("''","','")+")";

    cursor.execute(cluster_in_query);

    print(cluster_in_query)

conn.commit();


