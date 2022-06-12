import networkx as nx
import random
import pandas as pd
import matplotlib.pyplot as plt
import numpy

"""For this small world network, I suppose that the probability can make difference on the spread
of the epdimic."""
#Initialize parameters (we need N >> K >> lnN)
N = 1000 #total number of nodes, lnN = 7
k = 40 #average path length from each nodes
p = 0.01  #from p = 0 to p = 1, the model changes from regular to random

#Create the WS model
G_WS = nx.watts_strogatz_graph(N, k, p)

#check small world property
#sw_sigma = nx.sigma(G_WS, nrand = 3) #>1 then small world
#sw_omega = nx.omega(G_WS) #near 0 then small world
#print(sw_sigma)
#print(sw_omega)

#visualize model
pos = nx.circular_layout(G_WS)
nx.draw_networkx(G_WS, pos)
ax = plt.gca()
ax.margins(0.20)
plt.axis("off")
plt.show()

#verify WS model
d_mean = nx.average_shortest_path_length(G_WS)
print("<d>:", d_mean)

k_mean = 2 * G_WS.number_of_edges() / N
print("lnN/ln<k>:", numpy.log(N) / numpy.log(k_mean))

average_clustering = nx.average_clustering(G_WS)
print("<c>:", average_clustering)

#plot the graph related to c, d and p
def plot_together(N, k):
    p = numpy.logspace(-5, 0, num = 21)
    #Repeat the process 20 times to get the average value of c and d.
    #For clustering coefficient
    c_list_average = []
    d_list_average = []
    for i in p:
        c_list = []
        for l in range(20):
            G_WS_func = nx.watts_strogatz_graph(N,k,i)
            average_clustering = nx.average_clustering(G_WS_func)
            c_list.append(average_clustering)
        m = 0
        for v in range(20):
            m = m + c_list[v]
        c_list_average.append(m / 20)
        
    for i in p:
        d_list = []
        for d in range(20):
            G_WS_func = nx.watts_strogatz_graph(N,k,i)
            distance = nx.average_shortest_path_length(G_WS_func)
            d_list.append(distance)
        n = 0
        for b in range(20):
            n = n + d_list[b]
        d_list_average.append(n / 20)
        
    average_clustering_div = [i / c_list_average[0] for i in c_list_average]
    average_distance_div = [i / d_list_average[0] for i in d_list_average]
    
    plt.plot(p, average_clustering_div, label='average clustering')
    plt.plot(p, average_distance_div, label='average distance')
    plt.xscale('log')
    plt.xlabel('p')
    plt.legend()
    plt.show()

#Since when k = 10, the network is more similar to real network, so I decide to set k = 10 all time.
#plot_together(1000, 10)

"""Now find hubs"""
#by using degree
def degree_hub(G):
    k1 = 0
    node_list_deg=[]
    for i in range(0, 999):
        if G_WS.degree[i] >= k1:
            k1 = G.degree[i]
            node_list_deg.append((i, k1))
    lar_degree = node_list_deg[-1][1]
    lar_deg_node_list = []
    for k in range(1, 5):
        if node_list_deg[-k][1] == lar_degree:
            lar_deg_node_list.append(node_list_deg[-k][0])
    return lar_deg_node_list

#by using betweenness
def betweenness_hub(G):
    k2 = 0
    node_list_bet=[]
    betweenness = nx.betweenness_centrality(G, normalized=False) #betweenness is a dict
    betweenness_list_key = list(betweenness) 
    betweenness_list_value = list(betweenness.values())
    for i in range(0, 999):
        if betweenness_list_value[i] >= k2:
            k2 = betweenness_list_value[i]
            node_list_bet.append((betweenness_list_key[i], betweenness_list_value[i]))
    lar_bet = node_list_bet[-1][1]
    lar_bet_node_list = []
    for k in range(1, 5):
        if node_list_bet[-k][1] == lar_bet:
            lar_bet_node_list.append(node_list_bet[-k][0])
    return lar_bet_node_list

#by using eigenvector centrality
def evec_central_hub(G):
    k3 = 0
    node_list_evec = []
    eigenvector = nx.eigenvector_centrality(G)
    eigenvector_list_key = list(eigenvector)
    eigenvector_list_value = list(eigenvector.values())
    for i in range(0, 999):
        if eigenvector_list_value[i] >= k3:
            k3 = eigenvector_list_value[i]
            node_list_evec.append((eigenvector_list_key[i], eigenvector_list_value[i]))
    lar_evec = node_list_evec[-1][1]
    lar_evec_node_list = []
    for k in range(1, 5):
        if node_list_evec[-k][1] == lar_evec:
            lar_evec_node_list.append(node_list_evec[-k][0])
    return lar_evec_node_list

#by using clustering coefficient
def clustering_hub(G):
    k4 = 0
    node_list_clustering = []
    node_list_clustering_lar = []
    for k in range(0, 999):
        clustering = nx.clustering(G, k)
        node_list_clustering.append((k, clustering))
    for i in range(0, 999):
        if node_list_clustering[i][1] >= k4:
            k7 = node_list_clustering[i][1]
            node_list_clustering_lar.append((i, node_list_clustering[i][1]))
    lar_clustering = node_list_clustering[-1][1]
    lar_c_node_list = []
    for k in range(1, 5):
        if node_list_clustering[-k][1] == lar_clustering:
            lar_c_node_list.append(node_list_clustering[-k][0])
    return lar_c_node_list

"""I would like to simulate the SIR model"""
#initialize the state of all nodes
#for node in G_WS.nodes():
#   G_WS.nodes[node]["state"] = "S"

#The function that can update the state of any nodes
def updateNodeState(G, node, alpha, beta):
    if G.nodes[node]["state"] == "I": #infected
        p = random.random() #generate random number between 0 and 1
        if p < beta:   # recover rate
            G.nodes[node]["state"] = "R" #change state to R
    elif G.nodes[node]["state"] == "S": #susceptible
        p = random.random() #generate random number between 0 and 1
        k = 0  #to calculate number of infected neighbours
        for neibor in G.adj[node]: #go through every neighbous
            if G.nodes[neibor]["state"] == "I": #if infected, then k+=1
                k = k + 1
        if p < 1 - (1 - alpha)**k:  #S are infected, then change to I
            G.nodes[node]["state"] = "I" 

#The function updates states of all of nodes
def updateNetworkState(G, alpha, beta):
    for node in G: #go through every node to update their state
        updateNodeState(G, node, alpha, beta)
        
#To count the total number of different states people
def countSIR(G):
    S = 0
    I = 0
    for node in G:
        if G.nodes[node]["state"] == "S":
            S = S + 1
        elif G.nodes[node]["state"] == "I":
            I = I + 1
    return S,I, len(G.nodes) - S - I

#initializing the graph and parameters
#for k in range(1, 5):
#    G_WS.nodes[k]["state"] = "I"  #choose some nodes to be infected firstly
 
days = 50 #time
alpha = 0.05 #infected rate
beta = 0.05  #recover rate

#set color to different states of people
color_dict = {"S":"orange","I":"red","R":"green"}

SIR_list_total = []
for i in range(20):
    SIR_list = []
    for node in G_WS.nodes():
        G_WS.nodes[node]["state"] = "S"
    for k in range(1, 5):
        G_WS.nodes[k]["state"] = "I"
    for t in range(0,days):
        updateNetworkState(G_WS,alpha,beta)
        SIR_list.append(list(countSIR(G_WS)))
    SIR_list_total.append(SIR_list)

print(SIR_list_total)

average_SIR_list = []
for j in range(days):
    total_S = 0
    total_I = 0
    for i in range(len(SIR_list_total)):
        total_S += SIR_list_total[i][j][0]
        total_I += SIR_list_total[i][j][1]
    average_S = total_S / len(SIR_list_total)
    average_I = total_I / len(SIR_list_total)
    average_R = N - average_S - average_I
    average_SIR_list.append([average_S, average_I, average_R])
    
#visualizing the graph
df = pd.DataFrame(average_SIR_list,columns=["S","I","R"])
df.plot(figsize=(9,6),color=[color_dict.get(x) for x in df.columns])
plt.show()
