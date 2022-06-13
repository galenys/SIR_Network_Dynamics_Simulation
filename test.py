import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random
import numpy as np
from enum import Enum

N = 100
inf_rate = 0.2
rec_rate = 0.2

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

    def vaccinate_randomly(self, p):
        self.vaccinated = random.sample(range(self.N), round(self.N * p))

    def vaccinate_highest_degree_first(self, p):
        self.vaccinated = list(range(self.N))
        self.vaccinated.sort(reverse=True, key=self.get_degree)
        self.vaccinated = self.vaccinated[:round(self.N * p)]

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

    def average_disease_graph_recovered_rate(self, infected_proportion, iterations):
        no_vac_samples = []
        self.no_vaccines()
        for i in range(iterations):
            no_vac_samples.append(self.disease_graph(infected_proportion=infected_proportion))
            print(f"No vaccination: {i}")

        random_vac_samples = []
        self.vaccinate_randomly(0.5)
        for i in range(iterations):
            random_vac_samples.append(self.disease_graph(infected_proportion=infected_proportion))
            print(f"Random vaccination: {i}")

        hub_vac_samples = []
        self.vaccinate_highest_degree_first(0.5)
        for i in range(iterations):
            hub_vac_samples.append(self.disease_graph(infected_proportion=infected_proportion))
            print(f"Hub vaccination: {i}")

        plt.xticks([1,2,3], ["No Vaccines", "Random Vaccination", "Vaccinate Hubs"])
        plt.boxplot([no_vac_samples, random_vac_samples, hub_vac_samples])
        plt.legend()
        plt.show()

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

S = Simulation(
        # G=nx.gnp_random_graph(500, 0.005, seed=42),
        G=nx.scale_free_graph(500, seed=42).to_undirected(),
        N=500,
        inf_rate=0.1,
        rec_rate=0.1,
        vacc_inf_rate=0.001
        )

