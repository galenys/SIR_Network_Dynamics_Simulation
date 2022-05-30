import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random
import time
from enum import Enum

N = 100
B = 0.2
K = 0.2

class Status(Enum):
    SUS = 1
    INF = 2
    REC = 3


class Node:
    def __init__(self, index, status) -> None:
        self.index = index
        self.status = status


class Simulation:
    def __init__(self) -> None:
        self.nodeStatuses = [Status.SUS for _ in range(N)]
        self.G = nx.gnp_random_graph(N, 0.05, seed=42)

    def infect(self, index):
        self.nodeStatuses[index] = Status.INF

    def infect_list(self, indices):
        for i in indices:
            self.infect(i)

    def timestep(self):
        for node in self.G:
            if self.nodeStatuses[node] == Status.INF:
                for neighbour in self.G.neighbors(node):
                    if (random.random() < B and 
                            self.nodeStatuses[neighbour] == Status.SUS):
                        self.nodeStatuses[neighbour] = Status.INF
                if random.random() < K:
                    self.nodeStatuses[node] = Status.REC

    def get_col(self, status):
        if status == Status.SUS:
            return "blue"
        if status == Status.INF:
            return "red"
        return "green"

    def get_col_map(self):
        return [self.get_col(status) for status in self.nodeStatuses]

    def display_network(self):
        random_pos = nx.random_layout(self.G, seed=42)
        pos = nx.spring_layout(self.G, pos=random_pos)
        nx.draw(self.G, node_color=self.get_col_map(), pos=pos, with_labels=True)
        plt.show()

    def animate(self, _):
        random_pos = nx.random_layout(self.G, seed=42)
        pos = nx.spring_layout(self.G, pos=random_pos)
        nx.draw(self.G, node_color=self.get_col_map(), pos=pos, with_labels=True)
        # plt.show()
        self.timestep()

    def continuous_display(self):
        self.infect_list([i for i in range(20)])
        anim = animation.FuncAnimation(plt.gcf(), self.animate, frames=20, interval=20)
        plt.show()

S = Simulation()


