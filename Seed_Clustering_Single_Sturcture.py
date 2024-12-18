import inline as inline
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import psycopg2
from pandas import DataFrame

pd.options.mode.chained_assignment = None  # default='warn'

idvectorspace = 103;#94;#58;
correlation_threshold=0.2;

seed_cluster_point_mean_threshold =1;

conn = psycopg2.connect("dbname=anatomy user=anatomy password=anatomy");

cursor = conn.cursor()


def getDBData(query) :

    cursor = conn.cursor()

    cursor.execute(query);

    names = [ x[0] for x in cursor.description]

    data = DataFrame(cursor.fetchall(),columns=names );

    return data;



vsq = ("select name from vectorspace where idvectorspace = "+str(idvectorspace));
vectorspacename= getDBData(vsq)['name'].iloc[0]

print(vectorspacename);




def saveClustrersToDB(data):

    cursor = conn.cursor()

    cursor.execute(
        "INSERT INTO clustering  (name,idvectorspace,idvectorspacedisplay) values ('SeedClusteringSingleStructure(" + vectorspacename + ") ',"+str(idvectorspace)+","+str(idvectorspace)+") RETURNING idclustering")
    idclustering = cursor.fetchone()[0]

    for vector_idx, vector in data.iterrows():

        #print(vector)

        cursor.execute("Insert into clusterpoint (idclustering,iddelineation2d, clusterlabel, description) values ( "+str(idclustering)+"  , "+str(vector_idx) +", " + str(vector["cluster_label"])+", '" +vector["cluster_description"]  +"')");



    conn.commit();


# query = ("select * from crosstab ('select iddelineation2d, namedstructure.name, value "
# 		"from vectorcomponent, namedstructure "
# 		"where namedstructure.idnamedstructure  = vectorcomponent.idnamedstructure "
# 		"and idvectorspace =48 "
# 		"order by 1,2'::text) as vectorcomponent(id bigint, auditory double precision, entorhinal double precision, insular double precision, motor double precision, mpfc double precision, ofc double precision, "
#          "perirhinal double precision, somatosensory double precision, visual double precision)"
# 		);



query  = ("select iddelineation2d, namedstructure.name, value from vectorcomponent, namedstructure where namedstructure.idnamedstructure  = vectorcomponent.idnamedstructure and idvectorspace ="+str(idvectorspace)+" order by 1,2");

query = ( "select iddelineation2d , (case when parent is not null then parent||'/'|| nm else nm end) as name,value "
            "from "
            "(select iddelineation2d, (select  larger1  from ontology where Larger11  =namedstructure.name and namedstructure.name <> larger1 limit 1) parent , namedstructure.name nm, value "
            "from vectorcomponent, namedstructure where namedstructure.idnamedstructure  = vectorcomponent.idnamedstructure and idvectorspace ="+str(idvectorspace)+") st" )

print(query)


originalData = getDBData(query)




df = originalData.pivot(index='iddelineation2d', columns='name', values='value')

mn=df;


regions =mn.columns.tolist();
#Zero out all below the median
for region in regions:

    mn_region_column = mn.loc[:,region];
    regionmedian = mn_region_column .median();
    mn_region_column .where(mn_region_column  > regionmedian, 0, inplace=True);

    print(mn.loc[:,region])



meanNonZero = mn.replace(0, np.NaN).mean();
print(meanNonZero)


cols_ordered = df.columns.tolist();

cols_ordered.sort(key=str.lower);

#print(cols_ordered)

df = df.reindex(columns=cols_ordered)

#print(df);

localInjectionDensities = df[df.columns[~df.columns.isin(['iddelineation2d'])]]


#print(list(localInjectionDensities.columns))

injDensityCorr= localInjectionDensities.corr();

cortical_regions1 =df.columns.tolist();
cortical_regions2 =df.columns.tolist();

df["cluster_label" ]=0;
df["cluster_description" ]="";

print(type(injDensityCorr))



print(type(cortical_regions1))

clusterdata = pd.DataFrame(columns = df.columns.tolist())

cluster_label =1;


for region1 in  cortical_regions1:

    for region2 in cortical_regions2:

        if region1 == region2:

            correlation = injDensityCorr.at[region1,region2];

            if correlation >correlation_threshold:


                print(region1 + " : " + region2 + "  corr:" + str(correlation)+ "  "+str(meanNonZero[region1]) + "  " +str(meanNonZero[region2]))

                seed_cluster_points = df.loc[(df[region1] >meanNonZero[region1]*seed_cluster_point_mean_threshold) & (df[region2] >meanNonZero[region2]*seed_cluster_point_mean_threshold), :];

                seed_cluster_points.loc[:,"cluster_label"] = cluster_label;

                seed_cluster_points.loc[:,"cluster_description"] = region1 ;

                clusterdata =pd.concat([clusterdata,seed_cluster_points], ignore_index=False)

                cluster_label=cluster_label+1;

                #seed_cluster_points = df.where((df[region1] >meanNonZero[region1]) & (df[region2] >meanNonZero[region2]));

    cortical_regions2.pop(0);


print("check")

saveClustrersToDB(clusterdata)










