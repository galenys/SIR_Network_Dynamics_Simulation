import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random
import numpy as np
from enum import Enum


class Status(Enum):
    SUS = 1
    INF = 2
    REC = 3


class Node:
    def __init__(self, index, status) -> None:
        self.index = index
        self.status = status

class Simulation:
    def __init__(self, G, N, inf_rate, rec_rate, vacc_inf_rate):
        self.nodeStatuses = [Status.SUS for _ in range(N)]
        self.G = G
        self.N = N
        self.inf_rate = inf_rate
        self.rec_rate = rec_rate
        self.vacc_inf_rate = vacc_inf_rate
        self.vaccinated = []

    def initialise_population(self, infected):
        for i in range(self.N):
            if i in infected:
                self.nodeStatuses[i] = Status.INF
            else:
                self.nodeStatuses[i] = Status.SUS

    def no_vaccines(self):
        self.vaccinated = []

    def vaccinate_randomly(self):
        edge_nodes = self.edge_nodes()
        self.vaccinated = random.sample(range(self.N), len(edge_nodes))

    def edge_nodes(self):
        G = self.G
        set_nodes = set()
        for e in G.edges():
            if (e[0] <= 249 and e[1] >= 250) or (e[1] <= 249 and e[0] >= 250):
                set_nodes.add(e[0])
                set_nodes.add(e[1])
        return list(set_nodes)

    def vaccinate_edge_nodes(self):
        edge_nodes = self.edge_nodes()
        self.vaccinated = edge_nodes

    def get_degree(self, i):
        return self.G.degree[i]

    def timestep(self):
        for node in self.G:
            if self.nodeStatuses[node] == Status.INF:
                for neighbour in self.G.neighbors(node):
                    rate = self.vacc_inf_rate if neighbour in self.vaccinated else self.inf_rate
                    if (random.random() < rate and 
                            self.nodeStatuses[neighbour] == Status.SUS):
                        self.nodeStatuses[neighbour] = Status.INF
                if random.random() < self.rec_rate:
                    self.nodeStatuses[node] = Status.REC

    def get_col(self, status):
        if status == Status.SUS:
            return "blue"
        if status == Status.INF:
            return "red"
        return "green"

    def get_col_map(self):
        return [self.get_col(status) for status in self.nodeStatuses]

    def get_proportions(self):
        susceptible = 0
        infected = 0
        recovered = 0
        for status in self.nodeStatuses:
            if status == Status.SUS:
                susceptible += 1
            elif status == Status.INF:
                infected += 1
            else:
                recovered += 1
        return [susceptible/self.N, infected/self.N, recovered/self.N]

    def disease_graph(self, infected_proportion, should_display=False):
        self.initialise_population(infected=list(range(round(self.N*infected_proportion))))
        sus, inf, rec = [], [], []
        proportions = self.get_proportions()
        timesteps = 0
        while proportions[1] != 0:
            sus.append(proportions[0])
            inf.append(proportions[1])
            rec.append(proportions[2])
            self.timestep()
            proportions = self.get_proportions()
            timesteps += 1
        if should_display:
            X = np.linspace(0, timesteps, timesteps)
            plt.plot(X, sus, label="Susceptible")
            plt.plot(X, inf, label="Infected")
            plt.plot(X, rec, label="Recovered")
            plt.legend()
            plt.ylim((0, 1))
            plt.show()
        return proportions[2]

    def get_final_proportions(self):
        proportions = self.get_proportions()
        timesteps = 0
        while proportions[1] != 0:
            self.timestep()
            proportions = self.get_proportions()
            timesteps += 1
        return proportions[2]

    def average_disease_graph_recovered_rate(self, infected_proportion, iterations):
        total = 0
        for i in range(iterations):
            total += self.disease_graph(infected_proportion=infected_proportion)
            print(f"Iteration {i}")
        return total / iterations

    def display_network(self):
        random_pos = nx.random_layout(self.G, seed=42)
        pos = nx.spring_layout(self.G, pos=random_pos)
        nx.draw(self.G, node_color=self.get_col_map(), pos=pos, with_labels=True)
        plt.show()

    # def animate(self, _):
    #     random_pos = nx.random_layout(self.G, seed=42)
    #     pos = nx.spring_layout(self.G, pos=random_pos)
    #     nx.draw(self.G, node_color=self.get_col_map(), pos=pos, with_labels=True)
    #     # plt.show()
    #     self.timestep()

    # def continuous_display(self):
    #     self.initialise_population([i for i in range(20)])
    #     anim = animation.FuncAnimation(plt.gcf(), self.animate, frames=20, interval=20)
    #     plt.show()


N = 500
p_er = 0.01
l = N*(N-1)*p_er/2
r = 0.01

def p_2_sbm(N, r):
    n = N / 2
    sizes = [int(n), int(n)]
    p = l / (n*n*r + n*(n-1))
    return sizes, [[p, r*p], [r*p, p]]


n_sbm = [[], [], []]
sizes, p = p_2_sbm(N, r)

for i in range(1000):
    G = nx.stochastic_block_model(sizes, p)
    S = Simulation(G=G, N=N, inf_rate=0.05, rec_rate=0.05, vacc_inf_rate=0.005)
    S.initialise_population([random.randint(0, 249) for i in range(5)])
    S.no_vaccines()
    n = S.get_final_proportions()
    n_sbm[0].append(n)
    print(i)

for i in range(1000):
    G = nx.stochastic_block_model(sizes, p)
    S = Simulation(G=G, N=N, inf_rate=0.05, rec_rate=0.05, vacc_inf_rate=0.005)
    S.initialise_population([random.randint(0, 249) for i in range(5)])
    S.vaccinate_randomly()
    n = S.get_final_proportions()
    n_sbm[1].append(n)
    print(i)

for i in range(1000):
    G = nx.stochastic_block_model(sizes, p)
    S = Simulation(G=G, N=N, inf_rate=0.05, rec_rate=0.05, vacc_inf_rate=0.005)
    S.initialise_population([random.randint(0, 249) for i in range(5)])
    S.vaccinate_edge_nodes()
    n = S.get_final_proportions()
    n_sbm[2].append(n)
    print(i)


print(np.mean(n_sbm[0]), np.mean(n_sbm[1]), np.mean(n_sbm[2]))

fig, ax = plt.subplots()
ax.boxplot(n_sbm)
plt.xticks([1, 2, 3], ['No Vaccination', 'Random Vaccination', 'Targeted Vaccination'])
plt.ylabel('Proportion of population affected by illness')
plt.show()
