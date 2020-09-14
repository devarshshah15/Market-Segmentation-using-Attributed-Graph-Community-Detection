import sys
import numpy as np
import pandas as pd
from igraph import *
from scipy import spatial
import math



def SAC1_clustering_ALGO(all_communities):
        global graph
        cnt = 0
        
        for v in range(leng):
                cur_vertex = []
                
                found=False   
                for comm in all_communities:
                        if v in comm:
                                cur_vertex = comm
                                found = True
                        if found==True:
                                break

                max_g_val = -math.inf   
                max_comm  = []
                for comm in all_communities: 
                        wt = 0
                        comm1 = list(set(comm.copy()))
                        for points in comm1:
                                if graph.get_eid(int(v),int(points),error=False)!=-1:
                                        index = graph.get_eid(v,points)
                                        t = graph.es["weight"]
                                        wt += t[index]

                        
                        di = sum(graph.degree(comm1))
                        dj = graph.degree(v)

                        
                        delta_q_neuman = (wt - di * dj / (2*len(graph.es)))/(2*len(graph.es))

                        
                        delta_q_attri = 0.0
                        c = 0

                        ### calcuating delta_q_attri from section 4
                        while(c<len(comm1)):
                                ind_1 = comm1[c]
                                delta_q_attri += similarity_matrix[ind_1][v]
                                c+=1
                                
                        if len(comm1)!=0:
                                delta_q_attri = delta_q_attri/len(comm1)
                                delta_q_attri /= len(comm1)

                        
                        delta_q = alpha*delta_q_neuman
                        delta_q += (1-alpha)*delta_q_attri
                       
                        if delta_q>0:
                                if delta_q > max_g_val:
                                        max_g_val = delta_q
                                        max_comm = comm

                set_cur_vertex = set(cur_vertex)   ## getting unique
                set_max_comm = set(max_comm)

                if set_cur_vertex != set_max_comm:  
                        if max_g_val <= 0:
                                pass
                        else: 
                                if v in cur_vertex:    
                                        cur_vertex.remove(v)
                        max_comm.append(v)  

                        cnt=cnt+1

                        if(len(cur_vertex)==0):
                                all_communities.remove([])  

        return cnt 

def phase1(attri,all_communities,alpha):##,leng):
        global leng,similarity_matrix,graph
        iterations = 0
        similarity_matrix = []
        for i in range(leng):
                temp = []
                for j in range(leng):
                        temp.append(0)
                similarity_matrix.append(temp)
                        
        st = 0
        inside=st
        while(st<leng):
                inside = st
                while(inside<leng):
                        similarity_matrix[st][inside] = spatial.distance.cosine(list(graph.vs.select(st)[0].attributes().values()),list(graph.vs.select(inside)[0].attributes().values()))
                        similarity_matrix[st][inside] = 1/(1+similarity_matrix[st][inside])
                        similarity_matrix[inside][st] = similarity_matrix[st][inside]
                        inside+=1
                st+=1
                
        cnt = SAC1_clustering_ALGO(all_communities)
        while(cnt>0 and iterations<15):    
                cnt = SAC1_clustering_ALGO(all_communities)   
                iterations+=1
                

def phase2(similarity_matrix,all_communities,grp_of_vertices,alpha):
        global leng,graph
        pointer = 0;
        st=0
        while(st<len(all_communities)):
                inside=0
                while(inside<len(all_communities[st])):
                        grp_of_vertices[all_communities[st][inside]]= pointer
                        inside+=1
                pointer +=1
                st+=1

        
        graph.contract_vertices(grp_of_vertices,combine_attrs="mean")
        leng = len(all_communities)#pointer

       
        graph.simplify()
        
        
        all_communities = []
        for i in range(leng):
                all_communities.append([i])

        
        graph.es["weight"] = [0]
        for i in range(len(graph.es)-1):
                graph.es["weight"].append(0) 

        for e in edgs:
                if grp_of_vertices[e[0]] != grp_of_vertices[e[1]]:  ### updating weight here
                        graph.es["weight"][graph.get_eid(grp_of_vertices[e[0]],grp_of_vertices[e[1]])]+=1 #()

        similarity_matrix = []
        for i in range(leng):
                temp = []
                for j in range(leng):
                        temp.append(0)
                similarity_matrix.append(temp)
        st = 0
        inside=st
        while(st<leng):
                inside = st
                while(inside<leng):
                        similarity_matrix[st][inside] = spatial.distance.cosine(list(graph.vs.select(st)[0].attributes().values()),list(graph.vs.select(inside)[0].attributes().values()))
                        similarity_matrix[st][inside] = 1/(1+similarity_matrix[st][inside])
                        similarity_matrix[inside][st] = similarity_matrix[st][inside]
                        inside+=1
                st+=1

        
        phase1(attri,all_communities,alpha)


        
alpha = 0.0
if len(sys.argv)!= 2:
	print("Enter alpha value")
	sys.exit()
else:
	alpha = float(sys.argv[1])


attri = pd.read_csv("./data/fb_caltech_small_attrlist.csv")

graph = Graph()
graph.add_vertices(len(attri))


edges = open('./data/fb_caltech_small_edgelist.txt')
list_of_edges = edges.read().split("\n")


edgs = []
for edg in list_of_edges:
        if edg.split(' ')[0] != '' and edg.split(' ')[1] != '':
                edgs.append((int(edg.split(' ')[0]),int(edg.split(' ')[1])))

graph.add_edges(edgs)

graph.es["weight"] = [1]
for i in range(len(graph.es)-1):
        graph.es["weight"].append(1)


similarity_matrix = np.zeros((len(attri),len(attri)))

leng = len(attri)

for column in attri.keys():
        graph.vs[column] = attri[column]
	

all_communities = []
for i in range(leng):
        all_communities.append([i])
        

st = 0
inside=st
while(st<leng):
        inside = st
        while(inside<leng):
                ### using cosine function is spatial to calculate similarity
                similarity_matrix[st][inside] = spatial.distance.cosine(list(graph.vs.select(st)[0].attributes().values()),list(graph.vs.select(inside)[0].attributes().values()))
                similarity_matrix[st][inside] = 1/(1+similarity_matrix[st][inside])
                similarity_matrix[inside][st] = similarity_matrix[st][inside]
                inside+=1
        st+=1

phase1(attri,all_communities,alpha)#,similarity_matrix)#,len(attri))

grp_of_vertices = []
for x in range(len(attri)):
        grp_of_vertices.append(0)
        
phase2(similarity_matrix,all_communities,grp_of_vertices,alpha)#attri,list_of_edges,alpha)


if alpha==0.5:
        val = 5    
else:
        val = int(alpha)
output_file = open("./communities_"+str(val)+".txt","w")
for c in all_communities:
        output_file.write(",".join([str(i) for i in c]))
        output_file.write("\n")
output_file.close()

