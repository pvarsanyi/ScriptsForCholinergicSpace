import inline as inline
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import psycopg2
from pandas import DataFrame
from scipy.stats import pearsonr
import io
from sqlalchemy import create_engine

idvectorspace = 99;#87;#58;
Pvalue_treshold=0.05;
isPvalue=False

def getDBData(query) :

    cursor = conn.cursor()

    cursor.execute(query);

    names = [ x[0] for x in cursor.description]

    data = DataFrame(cursor.fetchall(),columns=names );

    return data;

def calculate_pvalues(df):
    dfcols = pd.DataFrame(columns=df.columns)
    pvalues = dfcols.transpose().join(dfcols, how='outer')
    pvalues = pvalues.apply(pd.to_numeric)
    for r in df.columns:
        for c in df.columns:
            tmp = df[df[r].notnull() & df[c].notnull()]
            pvalues[r][c] = pearsonr(tmp[r], tmp[c])[1]
    return pvalues



def saveDfToDB (df, tablename):
    engine = create_engine('postgresql+psycopg2://anatomy:anatomy@127.0.0.1:5432/anatomy')

    df.head(0).to_sql(tablename, engine, if_exists='replace',
                      index=False)  # drops old table and creates new empty table

    conn = engine.raw_connection()
    cur = conn.cursor()
    output = io.StringIO()
    df.to_csv(output, sep='\t', header=False, index=False)
    output.seek(0)
    contents = output.getvalue()
    cur.copy_from(output, tablename, null="")  # null values become ''
    conn.commit()


conn = psycopg2.connect("dbname=anatomy user=anatomy password=anatomy");

cursor = conn.cursor()

# query = ("select * from crosstab ('select iddelineation2d, namedstructure.name, value "
# 		"from vectorcomponent, namedstructure "
# 		"where namedstructure.idnamedstructure  = vectorcomponent.idnamedstructure "
# 		"and idvectorspace =48 "
# 		"order by 1,2'::text) as vectorcomponent(id bigint, auditory double precision, entorhinal double precision, insular double precision, motor double precision, mpfc double precision, ofc double precision, "
#          "perirhinal double precision, somatosensory double precision, visual double precision)"
# 		);


vsq = ("select name from vectorspace where idvectorspace = "+str(idvectorspace));

vectorspacename= getDBData(vsq)['name'].iloc[0]

print(vectorspacename);

query  = ("select iddelineation2d, namedstructure.name, value from vectorcomponent, namedstructure where namedstructure.idnamedstructure  = vectorcomponent.idnamedstructure and idvectorspace ="+str(idvectorspace)+" order by 1,2");

query = ( "select iddelineation2d , (case when parent is not null then parent||'/'|| nm else nm end) as name,value "
            "from "
            "(select iddelineation2d, (select  larger1  from ontology where Larger11  =namedstructure.name and namedstructure.name <> larger1 limit 1) parent , namedstructure.name nm, value "
            "from vectorcomponent, namedstructure where namedstructure.idnamedstructure  = vectorcomponent.idnamedstructure and idvectorspace ="+str(idvectorspace)+") st" )

print(query)


df = getDBData(query)






df = df.pivot(index='iddelineation2d', columns='name', values='value')


cols_ordered = df.columns.tolist();

cols_ordered.sort(key=str.lower);

print(cols_ordered)

df = df.reindex(columns=cols_ordered)

print(df);

localInjectionDensities = df[df.columns[~df.columns.isin(['iddelineation2d'])]]



#print(list(localInjectionDensities.columns))
injDensityCorr=pd.DataFrame();



if isPvalue:



    injDensityCorr= calculate_pvalues(localInjectionDensities);



    saveDfToDB(injDensityCorr.reset_index(), "Pvalue_of_Density_Corr_Overlap_larger11")

    #saveDfToDB (injDensityCorr, "InjectionCorrelationPValue")

    injDensityCorr.mask(injDensityCorr >=Pvalue_treshold ,100, inplace=True)
    injDensityCorr.mask(injDensityCorr <Pvalue_treshold ,float(1), inplace=True)
    injDensityCorr.mask(injDensityCorr ==100 ,float(-1), inplace=True)
else:

    injDensityCorr = localInjectionDensities.corr();



print(injDensityCorr )


##Creating De-pivoted table to save to db
renamedCorr = injDensityCorr.rename_axis(None).rename_axis(None, axis=1)
corrTable =  renamedCorr.stack().reset_index()
#saveDfToDB (corrTable, "Corr_"+ vectorspacename );

f, ax = plt.subplots(figsize=(len(injDensityCorr.index),len(injDensityCorr.columns) ))

sns.set(font_scale =3)
heatmap = sns.heatmap(injDensityCorr,
                      square = True ,
                      linewidths = 2,
                      cmap ='coolwarm',
                      #cbar_kws={"shrink": .4, 'ticks' : [-1, -.5, 0, 0.5, 1]},
                      vmin = -1,
                      vmax = 1,

                      annot = True,
                      annot_kws = {"size": 10}
                      )

heatmap.set_xlabel("Cortical Region", fontsize = 16)
heatmap.set_ylabel("Cortical Region", fontsize = 16)
heatmap.set_title("Spatial correlation between cell populations in BF ", fontsize = 16)

ax.set_yticklabels(injDensityCorr.index, rotation = 0,fontsize =15)
ax.set_xticklabels(injDensityCorr.columns,fontsize = 15)
sns.set_style({'xtick.bottom': True}, {'ytick.left': True})

if isPvalue:

    heatmap.get_figure().savefig('target_density_pvalue_' + vectorspacename + '.png', bbox_inches='tight')

else:
    heatmap.get_figure().savefig('target_density_correlation_'+vectorspacename+'.png', bbox_inches='tight')
#ax.set_yticklabels(corr_matrix.columns, rotation = 0)
#ax.set_xticklabels(corr_matrix.columns) #add the column names as labels










