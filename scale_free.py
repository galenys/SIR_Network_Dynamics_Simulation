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
            plt.ylabel("Proportion of population affected by illness")
            plt.xlabel("Time")
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

        plt.boxplot([no_vac_samples, random_vac_samples, hub_vac_samples])
        plt.ylabel("Proportion of population affected by illness")
        plt.xticks([1,2,3], ["No Vaccination", "Random Vaccination", "Targeted Vaccination"])
        plt.show()

    def display_network(self):
        random_pos = nx.random_layout(self.G, seed=42)
        pos = nx.spring_layout(self.G, pos=random_pos)
        nx.draw(self.G, node_color=self.get_col_map(), pos=pos, with_labels=True)
        plt.show()

S = Simulation(
        # G=nx.gnp_random_graph(500, 0.005, seed=42),
        # G=nx.scale_free_graph(500, seed=42).to_undirected(),
        G=nx.watts_strogatz_graph(500, 20, 0.0001),
        N=500,
        inf_rate=0.01,
        rec_rate=0.1,
        vacc_inf_rate=0.001
        )




def varied_parameter_simulations(iterations):
    simulationA = Simulation(
                G=power_law_graph(500, 2),
                N=500,
                inf_rate=0.1,
                rec_rate=0.1,
                vacc_inf_rate=0.001
            )
    samplesA = []
    for i in range(iterations):
        samplesA.append(simulationA.disease_graph(0.1))
        print(f"A: {i}")

    simulationB = Simulation(
                G=power_law_graph(500, 2.25),
                N=500,
                inf_rate=0.1,
                rec_rate=0.1,
                vacc_inf_rate=0.001
            )
    samplesB = []
    for i in range(iterations):
        samplesB.append(simulationB.disease_graph(0.1))
        print(f"B: {i}")


    simulationC = Simulation(
                G=power_law_graph(500, 2.5),
                N=500,
                inf_rate=0.1,
                rec_rate=0.1,
                vacc_inf_rate=0.001
            )
    samplesC = []
    for i in range(iterations):
        samplesC.append(simulationC.disease_graph(0.1))
        print(f"C: {i}")

    simulationD = Simulation(
                G=power_law_graph(500, 2.75),
                N=500,
                inf_rate=0.1,
                rec_rate=0.1,
                vacc_inf_rate=0.001
            )
    samplesD = []
    for i in range(iterations):
        samplesD.append(simulationD.disease_graph(0.1))
        print(f"D: {i}")

    simulationE = Simulation(
                G=power_law_graph(500, 3),
                N=500,
                inf_rate=0.1,
                rec_rate=0.1,
                vacc_inf_rate=0.001
            )
    samplesE = []
    for i in range(iterations):
        samplesE.append(simulationE.disease_graph(0.1))
        print(f"E: {i}")

    plt.boxplot([samplesA, samplesB, samplesC, samplesD, samplesE])
    plt.ylabel("Proportion of population affected by illness")
    plt.xticks([1, 2, 3, 4, 5], ["gamma=2", "gamma=2.25", "gamma=2.5", "gamma=2.75", "gamma=3"])
    plt.show()


def power_law_graph(n, exp):
    while True:  
        s = []
        while len(s) < n:
            nextval = int(nx.utils.powerlaw_sequence(1, exp)[0]) #100 nodes, power-law exponent 2.5
            if nextval != 0:
                s.append(nextval)
        if sum(s)%2 == 0:
            break
    G = nx.configuration_model(s)
    G = nx.Graph(G) # remove parallel edges
    G.remove_edges_from(nx.selfloop_edges(G))
    return G

