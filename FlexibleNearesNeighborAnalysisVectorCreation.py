from multiprocessing import Pool

import inline as inline
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
import psycopg2
from pandas import DataFrame
import warnings
from sqlalchemy import create_engine
import psycopg2
import io



conn = psycopg2.connect("dbname=anatomy user=anatomy password=anatomy");

def getDBData(query) :

    cursor = conn.cursor()

    cursor.execute(query);

    names = [ x[0] for x in cursor.description]

    data = DataFrame(cursor.fetchall(),columns=names );

    return data;

cortical_region_partition_base ="larger11";
exception_regions = ["S1BF (barrel field uniden)","M1 (Agl,PrCl) (uniden)"]

warnings.simplefilter(action='ignore', category=FutureWarning)

max_context_radius = 300;

normalization_by_innervation_densities= False;

normalization_by_overlapping_injections= True;

normalization_power=1;

cortical_regions = list(getDBData("select distinct "+cortical_region_partition_base+" from injectionsummary where "+cortical_region_partition_base + " not in "+str(exception_regions).replace("[","(").replace("]",")") )[cortical_region_partition_base]);

cortical_innervation_densities_query = "select target_parent,cellnumber_in_mm3 from cortical_innervation_densities ";


projectionQuery=("select crb.idbrain,crb.idfeature, crb.idfeaturetype, crb.featurename, crb.iddelineation,crb.iddelineation2d,crb.idmapping , crb.paxinos_idsection,"
                 "mpr.*, ST_X(geometry) as x, ST_Y(geometry) as  y,ST_Z(geometry) as z "
                 "from cell_populations_in_reference_brain crb, injectionsummary mpr "
                 "where mpr.idlabeling = crb.idlabeling "
                 "and ")+cortical_region_partition_base +"  in "+str(cortical_regions).replace("[","(").replace("]",")");

print(projectionQuery)

target_parent_map_query = "select distinct target_parent, "+cortical_region_partition_base+" as partition_base from matview_injection_max_participation_in_region ";


referencequery = "select iddelineation2d, ST_X(cell) as x, ST_Y(cell) as y,ST_Z(cell) as z from referenceCholinergicSystemInVRB";



def parallelize_dataframe(df, func, n_cores=24):
    df_split = np.array_split(df, n_cores)
    pool = Pool(n_cores)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df


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

def saveContextVectors(vectors):

    cursor = conn.cursor()

    #comp_list =vectors.columns;

    comp_list = cortical_regions;

    contextspace_name = "NearestNeighbor_"+cortical_region_partition_base;



    cursor.execute("INSERT INTO vectorspace (name) values ('"+contextspace_name+"') RETURNING idvectorspace")
    id_context_space = cursor.fetchone()[0]

    for vector_idx, vector in vectors.iterrows():

        for component in comp_list:

            cursor.execute("INSERT INTO vectorcomponent (idnamedstructure,iddelineation2d,value,idvectorspace) values ((SELECT idnamedstructure from namedstructure where name = '" +component+ "'), "+ str(vector["iddelineation2d"])+","+str(vector[component])+", "+str(id_context_space)+" ) ")





    conn.commit();


def save3DVectors(vectors):

    cursor = conn.cursor()

    comp_list =vectors.columns;

    comp_list = list(set(comp_list) - set(["iddelineation2d"]))

    contextspace_name = "3DReferenceXYZ";




    cursor.execute("INSERT INTO vectorspace (name) values ('"+contextspace_name+"') RETURNING idvectorspace")
    id_context_space = cursor.fetchone()[0]

    for vector_idx, vector in vectors.iterrows():

        for component in comp_list:

            cursor.execute("INSERT INTO vectorcomponent (idnamedstructure,iddelineation2d,value,idvectorspace) values ((SELECT idnamedstructure from namedstructure where name = '" +component+ "'), "+ str(vector["iddelineation2d"])+","+str(vector[component])+", "+str(id_context_space)+" ) ")





    conn.commit();


def proximity_measure (distance):

    return 1/math.pow(distance,1);






def createContextVector(points_of_measure):

    i=0;
    #print(len(referencecells.index))

    referencecells = getDBData(referencequery);
    projectioncells = getDBData(projectionQuery);
    injection_overlap = getDBData("select * from injection_overlap");

    projection_reference_cell_ratio=len(projectioncells.index)/len(referencecells.index)



    ###Adding projection cells  to the reference system

    for measureindex, measuringpoint in points_of_measure.iterrows():

        contextProjectionCells = pd.DataFrame(data=None, columns=projectioncells.columns);
        #print(len(contextProjectionCells.index))
        contextProjectionCells["proximity"] =0;
        contextProjectionCells["distance"] = 0;

        ############################# computing projectioncell proximity ##########################################
        number_of_projectioncells_in_radius = 0;
        for prindex, projectioncell in projectioncells.iterrows():

           distance   = math.sqrt( math.pow(measuringpoint["x"] - projectioncell["x"],2) + math.pow(measuringpoint["y"] - projectioncell["y"],2) + math.pow(measuringpoint["z"] - projectioncell["z"],2));

           if   (distance < max_context_radius):
               #print("Reference (" + str(measuringpoint["x"]) + "," + str(measuringpoint["y"]) + "," + str(
               #    measuringpoint["z"]) + ")" + "     Proj (" + str(projectioncell["x"]) + "," + str(
               #    projectioncell["y"]) + "," + str(projectioncell["z"]) + ")  Distance: " + str(distance))

               proximity=  1;

               projectioncell["proximity"] = proximity;

               projectioncell["distance"] = distance;

               #contextProjectionCells= contextProjectionCells.append(projectioncell);

               contextProjectionCells.loc[len(contextProjectionCells)] = projectioncell;

               #number_of_projectioncells_in_radius=number_of_projectioncells_in_radius +1;


        ############################# counting referencecells in proximity ##########################################
        number_of_referencecells_in_radius = 0;

        sum_reference_cell_proximity=0;
        for refindex2, referencecell in referencecells.iterrows():

           distance   = math.sqrt( math.pow(measuringpoint["x"] - referencecell["x"],2) + math.pow(measuringpoint["y"] - referencecell["y"],2) + math.pow(measuringpoint["z"] - referencecell["z"],2));

           if  (0 < distance) & (distance < max_context_radius):

               refproximity = proximity_measure(distance);

               sum_reference_cell_proximity=sum_reference_cell_proximity+refproximity;


        #if number_of_referencecells_in_radius != 0:
         #   local_projection_reference_cell_ratio = number_of_projectioncells_in_radius/ number_of_referencecells_in_radius;
        #else:
        #    local_projection_reference_cell_ratio =1000;


        ###################### removing measuring points where projection density is ot sufficient to makedensity measurments #################################


        # if local_projection_reference_cell_ratio < projection_reference_cell_ratio :
        #     print("delete "+ str(number_of_referencecells_in_radius)+ "    "+str(number_of_projectioncells_in_radius)+" "+ str(local_projection_reference_cell_ratio)+"  "+str(projection_reference_cell_ratio) );
        #
        #     points_of_measure.drop(index = measureindex,inplace=True);
        #     continue;
        #
        # else:
        #     print("keep "+ str(number_of_referencecells_in_radius)+ "    "+str(number_of_projectioncells_in_radius)+" "+ str(local_projection_reference_cell_ratio)+"  "+str(projection_reference_cell_ratio) );



        ############################################# Creating normalizing contants by overlapping injections #############################################################

        if normalization_by_overlapping_injections:

            injections_in_context = contextProjectionCells["idlabeling"].unique();

            normalizing_dict={}

            for inj1 in injections_in_context:



                try:
                    inj1_area  = (injection_overlap.loc[(injection_overlap.idlabeling_1 == inj1) & (injection_overlap.idlabeling_2 == inj1),"sumarea_1"]).values[0]

                    sum_overlap_area = 0;

                    for inj2 in injections_in_context:

                        overlap_area=0;

                        inj1_inj2_overlap = (injection_overlap.loc[(injection_overlap.idlabeling_1 == inj1) & (injection_overlap.idlabeling_2 == inj2),"intersection_area"])

                        if len(inj1_inj2_overlap)>0:
                            overlap_area = inj1_inj2_overlap.values[0]

                        sum_overlap_area = sum_overlap_area+overlap_area;

                    if sum_overlap_area != 0 :   #This can be 0 only if the injection does not have injection site added to the db

                        normalizing_constant = inj1_area/sum_overlap_area;

                    else:

                        normalizing_constant=1;


                    normalizing_dict[inj1] = normalizing_constant;

                    #print("before: "+str(contextProjectionCells.loc[contextProjectionCells['idlabeling'] == inj1, 'proximity']));
                    #before = contextProjectionCells.loc[contextProjectionCells['idlabeling'] == inj1, 'proximity'];
                    #print(contextProjectionCells['idlabeling']);

                    act_proximity = contextProjectionCells.loc[contextProjectionCells['idlabeling'] == inj1, 'proximity'];


                    act_proximity= act_proximity * pow(normalizing_constant,normalization_power);

                    contextProjectionCells.loc[contextProjectionCells['idlabeling'] == inj1, 'proximity'] = act_proximity;

                    print("Act proximity:" + str(proximity))

                except IndexError:

                    print("Injection IndexError:" + str(inj1))


                # if normalizing_constant < 1:
                #
                #     print("");
                #     print("");
                #     print( "before:" +str(before) + "   after: " + str(contextProjectionCells.loc[contextProjectionCells['idlabeling'] == inj1, 'proximity']) + "  normalizing_constant:" +str(normalizing_constant));

                #print(normalizing_constant);

            #print(contextProjectionCells)
            #print(injections_in_context);

        ################################## creating context vectors #################################################

        #innervated_cortical_regions_in_context = contextProjectionCells[cortical_region_partition_base].unique();

        #sum_proximities_by_cortical_region = contextProjectionCells.groupby(cortical_region_partition_base)['proximity'].sum();
        #projectioncell_count_by_cortical_region  = contextProjectionCells.groupby(cortical_region_partition_base)['iddelineation2d'].count();
        projectioncell_count_by_cortical_region = contextProjectionCells.groupby(cortical_region_partition_base)['proximity'].sum();  ##### adding overlapping normalization



        #overall_sum_proximities =projectioncell_count_by_cortical_region.sum(axis = 0, skipna = True);

        #print(overall_sum_proximities)

        points_of_measure.at[measureindex,"sum_reference_proximity"] =sum_reference_cell_proximity;

        for region in projectioncell_count_by_cortical_region.index:


            region_sum_proximities = projectioncell_count_by_cortical_region.at[region];


            if region_sum_proximities != np.NaN:

                context_vector_region_member = region_sum_proximities; #adding only the cortical region proximity without diviing by:  /overall_sum_proximities;

                points_of_measure.at[measureindex,region] = context_vector_region_member;



        i=i+1;

        print(__name__ + ":NN"+str(i));

    return points_of_measure;

#if __name__ == '__main__':

#points_of_measure = getDBData( "select iddelineation2d, ST_X(cell) as x, ST_Y(cell) as y,ST_Z(cell) as z from referenceCholinergicSystemInVRB");

##########save3DVectors(refcells);

projectioncells = getDBData(projectionQuery);

#points_of_measure = points_of_measure.append(projectioncells[['iddelineation2d', 'x', 'y', 'z']])
points_of_measure = projectioncells[['iddelineation2d', 'x', 'y', 'z']];

#points_of_measure = points_of_measure.head(20);

points_of_measure.reset_index(drop=True, inplace=True);


for colname in cortical_regions:

    points_of_measure[colname]=0;



#points_of_measure =parallelize_dataframe(points_of_measure, createContextVector);

points_of_measure = createContextVector(points_of_measure);


#saveDfToDB(points_of_measure, "points_of_measure");

#points_of_measure= getDBData( "select * from points_of_measure where iddelineation2d  not in (select iddelineation2d from referenceCholinergicSystemInVRB)");











################################################# creating the reference vectors ############################









        #   local_region_reference_proximity_ratio = regionMean / measuringpoint["sum_reference_proximity"];
        #
        #   median_region_proximities[rp_index] = median_region_proximities[rp_index] - innervation[innervation["target_parent"]==rp_index]["zscore"] *std_of_median_region_proximities;
        #
        # regionMean = measuringpoint[cortical_regions].replace(0, np.NaN).mean();

        #
        #
        #
        # if local_region_reference_proximity_ratio < global_region_reference_proximity_ratio:
        #
        #     points_of_measure.drop(index = measureindex,inplace=True);
        #
        #     dropped_points=dropped_points+1;





points_of_measure =points_of_measure.replace(np.NaN,0 );

saveContextVectors(points_of_measure);

    #saveDfToDB(points_of_measure,"ovelapping_injections_in_bf")










