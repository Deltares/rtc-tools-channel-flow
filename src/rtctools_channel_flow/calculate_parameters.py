import math
import numpy as np


class GetLinearSVVariables:
    def __init__(
        self,
        n_level_nodes=4,
        length=0,
        h_b_up=0,
        h_b_down=0,
        q_nominal=0,
        width=0,
        y_nominal=0,
        y_nominal_down=0,
        friction_coefficient=0,
    ):
        self.n_level_nodes = n_level_nodes
        self.length = length
        self.h_b_up = h_b_up
        self.h_b_down = h_b_down
        self.q_nominal = q_nominal
        self.width = width
        self.y_nominal = y_nominal
        self.y_nominal_down = y_nominal_down
        self.friction_coefficient = friction_coefficient

    def getVariables(self):
        s_b = (self.h_b_up - self.h_b_down) / self.length
        g_n = 9.80665
        dx = self.length / (self.n_level_nodes - 1)
        self.q0 = [self.q_nominal] * (self.n_level_nodes + 1)
        self.t0 = [self.width] * self.n_level_nodes
        self.y0 = np.linspace(
            self.y_nominal_down, self.y_nominal, self.n_level_nodes
        )
        self.a0 = [None] * self.n_level_nodes
        for node in range(self.n_level_nodes):
            self.a0[node] = self.t0[node] * self.y0[node]

        # Calculate self.v0
        self.v0 = [None] * (2 * self.n_level_nodes - 1)
        self.v0[0] = self.q0[0] / self.a0[0]
        for node in range(1, self.n_level_nodes):  # Q points
            self.v0[node * 2 - 1] = self.q0[node] / (
                (self.a0[node - 1] + self.a0[node]) / 2
            )
        for node in range(1, self.n_level_nodes - 1):  # H points
            self.v0[node * 2] = ((self.q0[node + 1] + self.q0[node]) / 2) / (
                self.a0[node]
            )
        self.v0[2 * self.n_level_nodes - 2] = (
            self.q0[self.n_level_nodes] / (self.a0[self.n_level_nodes - 1])
        )

        # Calculate self.p0
        self.p0 = [None] * self.n_level_nodes
        for node in range(self.n_level_nodes):
            self.p0[node] = self.t0[node] + 2 * self.y0[node]

        self.r = [None] * self.n_level_nodes
        for node in range(self.n_level_nodes):
            self.r[node] = self.a0[node] / self.p0[node]

        self.sf = [None] * (2 * self.n_level_nodes - 1)
        self.sf[0] = (
            math.pow(self.q0[0], 2) * math.pow(self.friction_coefficient, 2)
        ) / (math.pow((self.a0[0]), 2) * math.pow(self.r[0], 4 / 3))
        for node in range(1, self.n_level_nodes):  # Q points
            self.sf[node * 2 - 1] = (
                math.pow(self.q0[node], 2) * math.pow(self.friction_coefficient, 2)
            ) / (
                math.pow(((self.a0[node - 1] + self.a0[node]) / 2), 2)
                * math.pow(((self.r[node - 1] + self.r[node]) / 2), 4 / 3)
            )
        for node in range(1, self.n_level_nodes - 1):  # H points
            self.sf[node * 2] = (
                math.pow((self.q0[node + 1] + self.q0[node]) / 2, 2)
                * math.pow(self.friction_coefficient, 2)
            ) / (math.pow((self.a0[node]), 2) * math.pow(self.r[node], 4 / 3))
        self.sf[2 * self.n_level_nodes - 2] = (
            math.pow(self.q0[self.n_level_nodes], 2)
            * math.pow(self.friction_coefficient, 2)
        ) / (
            math.pow((self.a0[self.n_level_nodes - 1]), 2)
            * math.pow(self.r[self.n_level_nodes - 1], 4 / 3)
        )

        # Calculate self.c0
        self.c0 = [None] * (self.n_level_nodes + 1)
        self.c0[0] = math.sqrt(g_n * self.y0[0])
        for node in range(1, self.n_level_nodes):
            self.c0[node] = math.sqrt(g_n * (self.y0[node - 1] + self.y0[node]) / 2)
        self.c0[self.n_level_nodes] = math.sqrt(g_n * self.y0[self.n_level_nodes - 1])

        self.f0 = [None] * (2 * self.n_level_nodes - 1)
        self.f0[0] = self.v0[0] / self.c0[0]
        for node in range(1, self.n_level_nodes):  # Q points
            self.f0[node * 2 - 1] = self.v0[node * 2 - 1] / self.c0[node]
        for node in range(1, self.n_level_nodes - 1):  # H points
            self.f0[node * 2] = self.v0[node * 2 - 1] / (
                (self.c0[node + 1] + self.c0[node]) / 2
            )
        self.f0[2 * self.n_level_nodes - 2] = (
            self.v0[2 * self.n_level_nodes - 2] / self.c0[self.n_level_nodes]
        )

        self.dydx = [None] * (2 * self.n_level_nodes - 1)
        for node in range(2 * self.n_level_nodes - 1):  # Q points
            self.dydx[node] = (s_b - self.sf[node]) / (1 - math.pow((self.f0[node]), 2))

        # Calculate self.deltas
        self.delta = [None] * (self.n_level_nodes + 1)
        self.delta[0] = (2 * g_n / self.v0[0]) * (s_b - self.dydx[0])
        for node in range(1, self.n_level_nodes):
            self.delta[node] = (2 * g_n / self.v0[node * 2 - 1]) * (
                s_b - self.dydx[node * 2 - 1]
            )
        self.delta[self.n_level_nodes] = (
            2 * g_n / self.v0[2 * self.n_level_nodes - 2]
        ) * (s_b - self.dydx[2 * self.n_level_nodes - 2])

        # Calculate self.kappas
        self.kappa = [None] * self.n_level_nodes
        for node in range(self.n_level_nodes):
            self.kappa[node] = (
                7 / 3
                - ((4 * (self.a0[node])) / (3 * self.t0[node] * self.p0[node])) * 2
            )  # This 2 is the value of the PDE because of rectangular shape

        # Calculate self.f2
        self.f2 = [None] * self.n_level_nodes
        for node in range(self.n_level_nodes):
            self.f2[node] = (math.pow(self.v0[node * 2], 2) * self.t0[node]) / (
                g_n * (self.a0[node])
            )

        # Calculate self.gammas
        self.gamma = [None] * self.n_level_nodes
        self.gamma[0] = math.pow(self.v0[0], 2) * (
            (self.t0[1] - self.t0[0]) / dx
        ) + g_n * self.t0[0] * (
            (1 + self.kappa[0]) * s_b
            - (1 + self.kappa[0] - (self.kappa[0] - 2) * self.f2[0]) * self.dydx[0]
        )
        for node in range(1, self.n_level_nodes - 1):
            self.gamma[node] = math.pow(self.v0[node * 2], 2) * (
                (self.t0[node + 1] - self.t0[node - 1]) / (2 * dx)
            ) + g_n * self.t0[node] * (
                (1 + self.kappa[node]) * s_b
                - (1 + self.kappa[node] - (self.kappa[node] - 2) * self.f2[node])
                * self.dydx[node * 2]
            )
        self.gamma[self.n_level_nodes - 1] = math.pow(
            self.v0[2 * self.n_level_nodes - 2], 2
        ) * (
            (self.t0[self.n_level_nodes - 1] - self.t0[self.n_level_nodes - 2]) / dx
        ) + g_n * self.t0[self.n_level_nodes - 1] * (
            (1 + self.kappa[self.n_level_nodes - 1]) * s_b
            - (
                1
                + self.kappa[self.n_level_nodes - 1]
                - (self.kappa[self.n_level_nodes - 1] - 2)
                * self.f2[self.n_level_nodes - 1]
            )
            * self.dydx[2 * self.n_level_nodes - 2]
        )
