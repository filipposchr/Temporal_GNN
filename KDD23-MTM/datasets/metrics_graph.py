import networkx as nx
import itertools
import re
import math
import numpy as np
from scipy import stats
from scipy.stats import kstest

def count_motif_events():
       #Experiment 4: MSRE Results, counting of temporal motifs
       
       motif_2 = 0
       motif_3 = 0
       motif_4 = 0
       
       file_original = "out_CollegeMsg-orig"
       #file_original = "out_SMS-orig"
       #file_original = "out_superuser-orig"
       #file_original = "out_facebook-wall-orig"
       
       with open(file_original, "r") as f:
              for line in itertools.islice(f, 3, None):
                     fields = re.split("\s+",line)
                     if (len(fields[0]) == 4):
                            #print(fields[0], 2)
                            motif_2 = motif_2 + int(fields[1])
                     elif len(fields[0]) == 6:
                            #print(fields[0], 3)
                            motif_3 = motif_3 + int(fields[1])
                     elif (len(fields[0]) == 8):
                            #print(fields[0], 4)
                            motif_4 = motif_4 + int(fields[1])   
       
       total_motifs_original = motif_2 + motif_3 + motif_4
       
       motif_2_g= 0
       motif_3_g = 0
       motif_4_g = 0
       
       file_gen = "out_CollegeMsg-gen"
       #file_gen = "out_SMS-gen"
       #file_gen = "out_superuser-gen"
       #file_gen = "out_facebook-wall-gen"
       
       with open(file_gen, "r") as f:
              for line in itertools.islice(f, 3, None):
                     fields = re.split("\s+",line)
                     if (len(fields[0]) == 4):
                            motif_2_g = motif_2_g + int(fields[1])
                     elif len(fields[0]) == 6:
                            motif_3_g = motif_3_g + int(fields[1])
                     elif (len(fields[0]) == 8):
                            motif_4_g = motif_4_g + int(fields[1])  
                             
       total_motifs_gen = motif_2_g + motif_3_g + motif_4_g
       
       print("MSRE: ",  ((motif_2_g-motif_2) / motif_2_g)**2, ((motif_3_g-motif_3) / motif_3_g)**2, ((motif_4_g-motif_4) / motif_4_g)**2 )
       print("Ratio Number of events: ", total_motifs_original / total_motifs_gen)
       
def global_graph_statistics(): 
       #Experiment 2: Evaluation global graph metrics

       #ORIGINAL GRAPH
       
       file = "CollegeMsg.txt"
       #file = "facebook-wall.txt"
       #file = "stackoverflow.txt"
       #file = "SMS-A.txt"
       #file = "superuser.txt"
       
       with open(file, "r") as f:
           G_orig = nx.read_edgelist(f.readlines(),
                                create_using=nx.MultiDiGraph,
                                data=[("Timestamp", str)],
                                delimiter=" ",
                                nodetype=int)
       
       timestamps_orig = []
       timestamps_per_edge = []
       for u, v, t in G_orig.edges(data=True):
              timestamps_per_edge.append(len(G_orig.get_edge_data(u, v)))
              timestamps_orig.append(int(t['Timestamp']))
       
       max_events_on_edge_orig = max(timestamps_per_edge)
       #print("MAX ORIG : " , max_events_on_edge_orig)
       sorted_timestamps_orig = sorted(timestamps_orig)
       iet_list = [t - s for s, t in zip(sorted_timestamps_orig, sorted_timestamps_orig[1:])]
       mean_iet_orig = sum(iet_list) / len(iet_list)
       #print("Mean IET: ", mean_iet_orig)
       
       timespan_orig = max(timestamps_orig) - min(timestamps_orig)
       print("Timespan (days) : " , int(timespan_orig /  86400))
              
       
       #GENERATED GRAPH
       
       file = "CollegeMsg_gen.txt"
       #file = "facebook-wall-gen.txt"
       #file = "SMS-A-gen.txt"
       #file = "superuser-gen.txt"
       
       with open(file, "r") as f:
              for line in itertools.islice(f, 4, None):
                     G_gen = nx.read_edgelist(f.readlines(),
                     create_using=nx.MultiDiGraph, data=[("Timestamp", str)], delimiter=" ", nodetype=int)

       timestamps_gen = []
       timestamps_per_edge = []
       for u, v, t in G_gen.edges(data=True):
              timestamps_per_edge.append(len(G_gen.get_edge_data(u, v)))
              timestamps_gen.append(int(t['Timestamp']))
       
       max_events_on_edge_gen = max(timestamps_per_edge)
       print("Max events on edge GEN : " , max_events_on_edge_gen)
       sorted_timestamps_gen = sorted(timestamps_gen)
       iet_list = [t - s for s, t in zip(sorted_timestamps_gen, sorted_timestamps_gen[1:])]
       mean_iet_gen = sum(iet_list) / len(iet_list)
       
                     
       timespan_gen = max(timestamps_gen) - min(timestamps_gen)
       
       print("Ratio Num of edges: " , G_orig.number_of_edges() / G_gen.number_of_edges())
       print("Ratio Mean Degree: ", (2*G_orig.size()/G_orig.order()) / (2*G_gen.size()/G_gen.order())) 
       print("Ratio N-Components: " , nx.number_strongly_connected_components(G_orig) / nx.number_strongly_connected_components(G_gen))

       largest_orig = len(max(nx.strongly_connected_components(G_orig), key=len))
       largest_gen = len(max(nx.strongly_connected_components(G_gen), key=len))
       print("Ratio LCC: ",  largest_orig / largest_gen)     
       print("Ratio Timespan: ", timespan_orig / timespan_gen)
       print("Ratio Mean IET: ", mean_iet_orig / mean_iet_gen)
       print("Ratio Maximum events on edge: " , max_events_on_edge_gen / max_events_on_edge_orig)
       

def ks_statistics():
       #Experiment 3: Kolmogorov-Smirnov (KS) test
       
       #ORIGINAL GRAPH

       file = "CollegeMsg.txt"
       #file = "facebook-wall.txt"
       #file = "SMS-A.txt"
       #file = "superuser.txt"
       
       with open(file, "r") as f:
           G_orig = nx.read_edgelist(f.readlines()[0:],
                                create_using=nx.MultiDiGraph,
                                data=[("Timestamp", str)],
                                delimiter=" ",
                                nodetype=int)
                                
       #GENERATED GRAPH
       
       file = "CollegeMsg_gen.txt"
       #file = "superuser-gen.txt"
       #file = "facebook-wall-gen.txt"
       #file = "SMS-A-gen.txt"
      
       with open(file, "r") as f:
              for line in itertools.islice(f, 4, None):
                     G_gen = nx.read_edgelist(f.readlines()[0:],
                                       create_using=nx.MultiDiGraph,
                                       data=[("Timestamp", str)],
                                       delimiter=" ",
                                       nodetype=int)
                                       
       orig_in_degree = []
       gen_in_degree = []
       
       for e,d in G_orig.in_degree:
              orig_in_degree.append(d)
              
       for e,d in G_gen.in_degree:
              gen_in_degree.append(d)              
              
       orig_out_degree = []
       gen_out_degree = []
       
       for e,d in G_orig.out_degree:
              orig_out_degree.append(d)
              
       for e,d in G_gen.out_degree:
              gen_out_degree.append(d)              
       
       timestamps_orig = []
       timestamps_gen = []
       
       for u, v, t in G_orig.edges(data=True):
              timestamps_orig.append(int(t['Timestamp']))
              
       sorted_timestamps_orig = sorted(timestamps_orig)
       iet_list_orig = [t - s for s, t in zip(sorted_timestamps_orig, sorted_timestamps_orig[1:])]
       
       for u, v, t in G_gen.edges(data=True):
              timestamps_gen.append(int(t['Timestamp']))
      
       sorted_timestamps_gen = sorted(timestamps_gen)
       iet_list_gen = [t - s for s, t in zip(sorted_timestamps_gen, sorted_timestamps_gen[1:])]
       
       print("KS in-degree: ", kstest(gen_in_degree, orig_in_degree))
       print("KS out-degree: ", kstest(orig_out_degree, gen_out_degree))
       print("KS IET: ", kstest(iet_list_orig, iet_list_gen))       
       print("KS Timestamp: ", kstest(timestamps_orig, timestamps_gen))

def global_graph_statistics2():
       #Statistics of the temporal network datasets - Generating Table 1

       #ORIGINAL GRAPH

       file = "CollegeMsg.txt"
       #file = "facebook-wall.txt"
       #file = "stackoverflow.txt"
       #file = "SMS-A.txt"
       #file = "superuser.txt"  
       
       with open(file, "r") as f:
           G_orig = nx.read_edgelist(f.readlines(),
                                create_using=nx.MultiDiGraph,
                                data=[("Timestamp", str)],
                                delimiter=" ",
                                nodetype=int)
       
       timestamps_orig = []
       timestamps_per_edge = []
       for u, v, t in G_orig.edges(data=True):
              timestamps_per_edge.append(len(G_orig.get_edge_data(u, v)))
              timestamps_orig.append(int(t['Timestamp']))
       
       max_events_on_edge_orig = max(timestamps_per_edge)
       #print("MAX ORIG : " , max_events_on_edge_orig)
       sorted_timestamps_orig = sorted(timestamps_orig)
       iet_list = [t - s for s, t in zip(sorted_timestamps_orig, sorted_timestamps_orig[1:])]
       mean_iet_orig = sum(iet_list) / len(iet_list)
       print("Nodes : " , G_orig.number_of_nodes())
       
       
       print("Edges : " , G_orig.number_of_edges())
       
       print("Mean IET: ", mean_iet_orig)
       
       timespan_orig = max(timestamps_orig) - min(timestamps_orig)
       print("Timespan (days) : " , int(timespan_orig /  86400))
              

global_graph_statistics2()             
ks_statistics()
count_motif_events()
global_graph_statistics()

