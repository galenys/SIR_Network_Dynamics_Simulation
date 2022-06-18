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
    def __init__(self, G, N, inf_rate, rec_rate):
        self.nodeStatuses = [Status.SUS for _ in range(N)]
        self.G = G
        self.N = N
        self.inf_rate = inf_rate
        self.rec_rate = rec_rate

    def initialise_population(self, infected):
        for i in range(self.N):
            if i in infected:
                self.nodeStatuses[i] = Status.INF
            else:
                self.nodeStatuses[i] = Status.SUS

    def timestep(self):
        for node in self.G:
            if self.nodeStatuses[node] == Status.INF:
                for neighbour in self.G.neighbors(node):
                    if (random.random() < self.inf_rate and 
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

    def disease_graph(self, infected_proportion):
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
    #     self.infect_list([i for i in range(20)])
    #     anim = animation.FuncAnimation(plt.gcf(), self.animate, frames=20, interval=20)
    #     plt.show()


N = 500
p_er = 0.01
l = N*(N-1)*p_er/2


def p_2_sbm(N, r):
    n = N / 2
    sizes = [int(n), int(n)]
    p = l / (n*n*r + n*(n-1))
    return sizes, [[p, r*p], [r*p, p]]


def p_3_sbm(N, r):
    n = N / 3
    sizes = [int(n), int(n), int(n)]
    p = l / (3*r*n*n + 3*n*(n-1)/2)
    return sizes, [[p, r*p, r*p], [r*p, p, r*p], [r*p, r*p, p]]


def p_4_sbm(N, r):
    n = N / 4
    sizes = [int(n), int(n), int(n), int(n)]
    p = l / (6*r*n*n + 2*n*(n-1))
    return sizes, [[p, r*p, r*p, r*p], [r*p, p, r*p, r*p], [r*p, r*p, p, r*p], [r*p, r*p, r*p, p]]


list_r = np.logspace(-6, 0, base=10.0, num=51)
n_sbm = []
r_sbm = []
n_m_sbm = []
r_m_sbm = []
n_v_sbm = []

for r in list_r:
    sizes, p = p_4_sbm(N, r)
    list_n = []
    for i in range(50):
        G = nx.stochastic_block_model(sizes, p)
        S = Simulation(G=G, N=N, inf_rate=0.05, rec_rate=0.05)
        S.initialise_population([random.randint(0, 124) for i in range(5)])
        n = S.get_final_proportions()
        list_n.append(n)
        n_sbm.append(n)
        r_sbm.append(np.log10(r))
    r_m_sbm.append(np.log10(r))
    n_m_sbm.append(np.mean(list_n))
    n_v_sbm.append(np.std(list_n))

plt.scatter(r_sbm, n_sbm, alpha=0.25)
plt.plot(r_m_sbm, n_m_sbm, 'b', linestyle=':', label='mean value')
plt.xlabel('log(q/p)')
plt.ylabel('Proportion of population affected by illness')
plt.legend(loc='best')
plt.show()

plt.plot(r_m_sbm, n_v_sbm)
plt.xlabel('log(q/p)')
plt.ylabel('Proportion of population affected by illness')
plt.show()
