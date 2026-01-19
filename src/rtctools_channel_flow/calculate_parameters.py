import math
import numpy as np


def CircularChannelData(y, r, q):
    """
    Compute hydraulic properties for a partially filled circular channel.

    Parameters
    ----------
    y : float
        Flow depth measured from the channel bottom [m].
    r : float
        Radius of the circular channel [m].
    q : float
        Discharge through the channel [m³/s].

    Returns
    -------
    alpha : float
        Central angle (rad) corresponding to the wetted portion of the channel.
    area : float
        Wetted cross-sectional area [m²].
    perimeter : float
        Wetted perimeter [m].
    hydraulic_radius : float
        Hydraulic radius, defined as area divided by wetted perimeter [m].
    c : float
        Kinematic wave celerity, computed as sqrt(g * hydraulic_radius) [m/s].
    v : float
        Mean flow velocity [m/s].

    Notes
    -----
    The circular cross-section is assumed to be partially filled.
    """
    diameter = 2 * r
    alpha = np.arccos(1 - y / r)
    area = (diameter**2) / 4 * (alpha - np.sin(2 * alpha) / 2)
    perimeter = alpha * diameter
    hydraulic_radius = (diameter / 4) * (1 - np.sin(2 * alpha) / (2 * alpha))
    celerity = np.sqrt(9.81 * hydraulic_radius)
    velocity = q / area

    return alpha, area, perimeter, hydraulic_radius, celerity, velocity


def froudeC(y, q, n, r, m):
    """
    Compute the Froude number for flow in a circular channel.

    Parameters
    ----------
    y : float
        Flow depth [m].
    q : float
        Discharge [m³/s].
    n : float
        Manning roughness coefficient [-].
        (Included for interface consistency; not used in the calculation.)
    r : float
        Radius of the circular channel [m].
    m : float
        Side slope parameter.
        (Not used for circular channels; included for interface consistency.)

    Returns
    -------
    fr : float
        Froude number [-].

    Notes
    -----
    The Froude number is defined as the ratio of mean velocity to
    kinematic wave celerity:

        Fr = v / c

    where c = sqrt(g * hydraulic_radius).
    """
    g = 9.81
    diameter = 2 * r
    alpha = np.arccos(1 - y / r)
    area = (diameter**2) / 4 * (alpha - np.sin(2 * alpha) / 2)
    hydraulic_radius = diameter / 4 * (1 - np.sin(2 * alpha) / (2 * alpha))
    celerity = np.sqrt(g * hydraulic_radius)
    velocity = q / area
    fr = velocity / celerity

    return fr


def froude(y, q, n, b0, m):
    """
    Compute the Froude number for flow in a trapezoidal channel.

    Parameters
    ----------
    y : float
        Flow depth [m].
    q : float
        Discharge [m³/s].
    n : float
        Manning roughness coefficient [-].
        (Included for interface consistency; not used in the calculation.)
    b0 : float
        Bottom width of the channel [m].
    m : float
        Side slope of the channel (horizontal/vertical) [-].

    Returns
    -------
    fr : float
        Froude number [-].

    Notes
    -----
    The kinematic wave celerity is computed as:

        c = sqrt(g * A / T)

    where A is the cross-sectional area and T is the top width.
    """
    g = 9.81
    area = b0 * y + m * y**2
    top_width = b0 + 2 * y * m
    celerity = np.sqrt(g * area / top_width)
    velocity = q / area
    fr = velocity / celerity

    return fr


def sf0C(y, q, n, r, m):
    """
    Compute the friction slope for a circular channel using Manning's equation.

    Parameters
    ----------
    y : float
        Flow depth [m].
    q : float
        Discharge [m³/s].
    n : float
        Manning roughness coefficient [-].
    r : float
        Radius of the circular channel [m].
    m : float
        Side slope parameter.
        (Not used for circular channels; included for interface consistency.)

    Returns
    -------
    sf0_ : float
        Friction slope [-].

    Notes
    -----
    The friction slope is computed using Manning's equation:

        Sf = (Q² n²) / (A² R^(4/3))
    """
    diameter = 2 * r
    alpha = np.arccos(1 - y / r)
    area = (diameter**2) / 4 * (alpha - np.sin(2 * alpha) / 2)
    hydraulic_radius = diameter / 4 * (1 - np.sin(2 * alpha) / (2 * alpha))
    sf0_ = (q**2 * n**2) / (area**2 * hydraulic_radius ** (4 / 3))

    return sf0_


def sf0(y, q, n, b0, m):
    """
    Compute the friction slope for a trapezoidal channel using Manning's equation.

    Parameters
    ----------
    y : float
        Flow depth [m].
    q : float
        Discharge [m³/s].
    n : float
        Manning roughness coefficient [-].
    b0 : float
        Bottom width of the channel [m].
    m : float
        Side slope of the channel (horizontal/vertical) [-].

    Returns
    -------
    sf0_ : float
        Friction slope [-].
    """
    area = b0 * y + m * y**2
    perimeter = b0 + 2 * y * np.sqrt(1 + m**2)
    hydraulic_radius = area / perimeter
    sf0_ = (q**2 * n**2) / (area**2 * hydraulic_radius ** (4 / 3))

    return sf0_


def normal_depth(q, n, b0, m, sb, yx, L, shape):
    """
    Compute the normal flow depth using a bisection method.

    Parameters
    ----------
    q : float
        Discharge [m³/s].
    n : float
        Manning roughness coefficient [-].
    b0 : float
        Bottom width (trapezoidal) or radius (circular) [m].
    m : float
        Side slope (trapezoidal) [-].
    sb : float
        Channel bed slope [-].
    yx : float
        Reference depth used to define the upper search bound [m].
    L : float
        Channel length [m].
        (Not used in the current implementation.)
    shape : int
        Channel shape identifier:
        - 0 : Trapezoidal channel
        - 1 : Circular channel

    Returns
    -------
    yn : float
        Normal flow depth [m].

    Notes
    -----
    The normal depth is determined by solving Manning's equation
    such that the friction slope equals the bed slope.
    A fixed number of bisection iterations (25) is used.
    """
    if shape == 0:
        y1 = 0
        y2 = yx * 5
        for k in range(25):
            y = (y1 + y2) / 2
            a = b0 * y + m * y**2
            p = b0 + 2 * y * (1 + m**2) ** 0.5
            r = a / p
            dif = (q**2 * n**2) / (a**2 * r ** (4 / 3)) - sb
            if dif < 0:
                y2 = (y1 + y2) / 2
            else:
                y1 = (y1 + y2) / 2
    elif shape == 1:
        y1 = 0
        y2 = yx * 5
        for k in range(25):
            y = (y1 + y2) / 2
            alpha, a, p, r, c, v = CircularChannelData(y, b0, q)
            dif = (q**2 * n**2) / (a**2 * r ** (4 / 3)) - sb
            if dif < 0:
                y2 = (y1 + y2) / 2
            else:
                y1 = (y1 + y2) / 2

    yn = y  # Normal depth

    return yn


def calculate_p11_inf(t, c, v, fr, alpha, gamma, x):
    p11_inf = (
        1
        / (t * c * (1 - fr))
        * (
            (1 + ((1 - fr) / (1 + fr)) ** 2 * np.exp(alpha * x))
            / (1 + np.exp(alpha * x))
        )
        ** 0.5
    )
    return p11_inf


def calculate_p12_inf(t, c, v, fr, alpha, gamma, x):
    p12_inf = (
        2
        / (t * c * (1 - fr**2))
        * (
            np.exp(-(gamma) / (2 * t * (c**2 - v**2)) * x)
            / (1 + np.exp(alpha * x)) ** 0.5
        )
    )
    return p12_inf


def calculate_p21_inf(t, c, v, fr, alpha, gamma, x):
    p21_inf = (
        2
        / (t * c * (1 - fr**2))
        * (
            np.exp((gamma) / (2 * t * (c**2 - v**2)) * x)
            / (1 + np.exp(alpha * x)) ** 0.5
        )
    )
    return p21_inf


def calculate_p22_inf(t, c, v, fr, alpha, gamma, x):
    p22_inf = (
        1
        / (t * c * (1 + fr))
        * (
            (1 + ((1 + fr) / (1 - fr)) ** 2 * np.exp(alpha * x))
            / (1 + np.exp(alpha * x))
        )
        ** 0.5
    )
    return p22_inf


def backwater_area_downstream_direction(t, c, v, gamma, x):
    ad = (t**2 * (c**2 - v**2)) / gamma * (1 - np.exp(-gamma / (t * (c**2 - v**2)) * x))
    return ad


def backwater_area_upstream_direction(t, c, v, gamma, x):
    au = (t**2 * (c**2 - v**2)) / gamma * (np.exp(gamma / (t * (c**2 - v**2)) * x) - 1)
    return au


def alpha_kappa(a, hydraulic_radius, p, m, fr, sb):
    kappa = 7 / 3 - 4 * a / (3 * hydraulic_radius * p) * 2 * (1 + m**2) ** 0.5
    alpha = (hydraulic_radius * (2 + (kappa - 1) * fr**2) * sb) / (a * fr * (1 - fr**2))
    return alpha, kappa


def calculate_gamma(hydraulic_radius, kappa, sb, fr, sx, v, m):
    g = 9.81
    gamma = v**2 * 2 * m * sx + (
        g
        * hydraulic_radius
        * ((1 + kappa) * sb - (1 + kappa - fr**2 * (kappa - 2)) * sx)
    )
    return gamma


def IdzFun(q, n, B, m, Sb, Y0, L, shape):
    """
    Compute parameters for an IDZ (Impulse-Delay-Zero) channel flow model.

    This function derives linearized hydraulic parameters for open-channel
    flow based on gradually varied flow theory. The channel is divided into
    upstream and downstream sections, and frequency-domain as well as
    time-delay properties are computed.

    The implementation supports trapezoidal and circular channel geometries
    and is intended for use in reduced-order channel flow models such as
    ``Deltares.ChannelFlow.Hydraulic.Branches.IDZ``.

    Parameters
    ----------
    q : float
        Discharge through the channel [m³/s].
    n : float
        Manning roughness coefficient [-].
    B : float
        Bottom width of the channel [m] (trapezoidal) or
        radius of the channel [m] (circular).
    m : float
        Side slope of the channel (horizontal/vertical) [-].
        For circular channels this parameter is retained for interface
        consistency.
    Sb : float
        Channel bed slope [-].
    Y0 : float
        Reference water depth, typically the downstream water depth [m].
    L : float
        Channel length [m].
    shape : int
        Channel geometry identifier:
        - ``0`` : Trapezoidal channel
        - ``1`` : Circular channel

    Returns
    -------
    p11 : float
        High-frequency upstream-to-upstream transfer coefficient [-].
    p12 : float
        High-frequency downstream-to-upstream transfer coefficient [-].
    p21 : float
        High-frequency upstream-to-downstream transfer coefficient [-].
    p22 : float
        High-frequency downstream-to-downstream transfer coefficient [-].
    au : float
        Low-frequency upstream storage coefficient [s²].
    ad : float
        Low-frequency downstream storage coefficient [s²].
    tu : float
        Upstream wave travel time [s].
    td : float
        Downstream wave travel time [s].
    yn : float
        Normal flow depth [m].
    x2 : float
        Location of the transition point between upstream and downstream
        regions, measured from the upstream end [m].

    Notes
    -----
    - The bottom slope is constrained to a minimum value to avoid numerical
      singularities.
    - Normal depth is computed using Manning’s equation and a bisection method.
    - The channel is split into upstream and downstream regions depending on
      the relationship between the normal depth and downstream depth.
    - Both low-frequency (storage, delay) and high-frequency (transfer matrix)
      components are computed.
    - Gravitational acceleration is assumed to be 9.81 m/s².

    References
    ----------
    The formulation follows standard linearized open-channel flow theory
    used in IDZ-type reduced-order hydraulic models. The detailed description can be found in
    Litrico, X., & Fromion, V. (2004). Analytical approximation of open-channel flow for controller
    design. Applied Mathematical Modelling, 28(7), 677-695.
    """
    # ------------------------------------------------------------------
    # Initialization
    # ------------------------------------------------------------------
    sb = Sb
    b0 = B  # Bottom width or radius
    h = Y0  # Initial / downstream depth
    yx = h  # Downstream depth
    g = 9.81  # Gravitational acceleration [m/s^2]

    yn = normal_depth(q, n, B, m, sb, Y0, L, shape)

    # ------------------------------------------------------------------
    # Water surface slope at downstream end
    # ------------------------------------------------------------------
    if shape == 0:
        sxx = (sb - sf0(yx, q, n, b0, m)) / (1 - froude(yx, q, n, b0, m) ** 2)
    elif shape == 1:
        sxx = (sb - sf0C(yx, q, n, b0, m)) / (1 - froudeC(yx, q, n, b0, m) ** 2)

    # ------------------------------------------------------------------
    # Location of transition point
    # ------------------------------------------------------------------

    x1 = max(0, L - (yx - yn) / sxx)

    if x1 == 0:
        y1 = yx - sxx * L
    else:
        y1 = yn

    x2 = (L + x1) / 2

    if x1 == 0:
        y2 = yx - (L - x2) * sxx
    else:
        y2 = yn + (x2 - x1) * sxx

    sxu = 0  # The upstream slope of the water surface wrt bottom is zero, as it is normal flow

    # ------------------------------------------------------------------
    # Part 2: Upstream section
    # ------------------------------------------------------------------
    if x1 != 0:
        if shape == 0:
            a = b0 * y1 + m * y1**2
            p = b0 + 2 * y1 * (1 + m**2) ** 0.5
            t = b0 + 2 * y1 * m
            c = (g * a / t) ** 0.5
            v = q / a
            fr = froude(y1, q, n, b0, m)
            alpha, kappa = alpha_kappa(a, t, p, m, fr, sb)
            gamma = calculate_gamma(t, kappa, sb, fr, sxu, v, m)
        elif shape == 1:
            alpha, a, p, t, c, v = CircularChannelData(y1, b0, q)
            fr = froudeC(y1, q, n, b0, m)
            alpha, kappa = alpha_kappa(a, t, p, m, fr, sb)
            gamma = calculate_gamma(t, kappa, sb, fr, sxu, v, 0)

        # Low frequencies
        x = x1
        td = x / (c + v)
        tu = x / (c - v)

        ad = backwater_area_downstream_direction(t, c, v, gamma, x)
        au = backwater_area_upstream_direction(t, c, v, gamma, x)

        # High frequencies
        p11_inf = calculate_p11_inf(t, c, v, fr, alpha, gamma, x)
        p12_inf = calculate_p12_inf(t, c, v, fr, alpha, gamma, x)
        p21_inf = calculate_p21_inf(t, c, v, fr, alpha, gamma, x)
        p22_inf = calculate_p22_inf(t, c, v, fr, alpha, gamma, x)

    # Part 3 ----Downstream part------------------------------------------------
    if shape == 0:
        a = b0 * y2 + m * y2**2
        p = b0 + 2 * y2 * (1 + m**2) ** 0.5
        # r = a / p
        t = b0 + 2 * y2 * m
        c = (g * a / t) ** 0.5
        v = q / a
        fr = v / c

        sxd = (sb - sf0(yx, q, n, b0, m)) / (1 - (froude(yx, q, n, b0, m)) ** 2)
    elif shape == 1:
        m = 0
        alpha, a, p, t, c, v = CircularChannelData(y1, b0, q)
        fr = froudeC(y1, q, n, b0, m)
        alpha, kappa = alpha_kappa(a, t, p, m, fr, sb)
        sxd = (sb - sf0(yx, q, n, b0, m)) / (1 - (froude(yx, q, n, b0, m)) ** 2)

    dtdy = 2 * m

    x = L - x1
    kappa = alpha_kappa(a, t, p, m, fr, sb)[1]
    alpha = (
        t
        / (a * fr * (1 - fr**2))
        * (
            (2 + (kappa - 1) * fr**2) * sb
            - (2 + (kappa - 1) * fr**2 - (a / t**2 * dtdy + kappa - 2) * fr**4) * sxx
        )
    )

    gamma = calculate_gamma(t, kappa, sb, fr, sxd, v, m)

    # Low frequencies
    td_ = x / (c + v)
    tu_ = x / (c - v)

    ad_ = backwater_area_downstream_direction(t, c, v, gamma, x)
    au_ = backwater_area_upstream_direction(t, c, v, gamma, x)

    # High frequencies
    p11_inf_ = calculate_p11_inf(t, c, v, fr, alpha, gamma, x)
    p12_inf_ = calculate_p12_inf(t, c, v, fr, alpha, gamma, x)
    p21_inf_ = calculate_p21_inf(t, c, v, fr, alpha, gamma, x)
    p22_inf_ = calculate_p22_inf(t, c, v, fr, alpha, gamma, x)

    # Part 4 -------------------Global model-----------------------------------
    if x1 == 0:  # If there is only backwater part
        td_hat = td_
        tu_hat = tu_

        ad_hat = ad_  # *(1+(ad/au_))
        au_hat = au_  # *(1+(au_/ad))

        p11_inf_hat = p11_inf_
        p12_inf_hat = p12_inf_
        p21_inf_hat = p21_inf_
        p22_inf_hat = p22_inf_

    else:
        td_hat = td + td_
        tu_hat = tu + tu_

        ad_hat = ad_ * (1 + (ad / au_))
        au_hat = au * (1 + (au_ / ad))

        p11_inf_hat = p11_inf + (p12_inf * p21_inf) / (p11_inf_ + p22_inf)
        p12_inf_hat = (p12_inf * p12_inf_) / (p11_inf_ + p22_inf)
        p21_inf_hat = (p21_inf * p21_inf_) / (p11_inf_ + p22_inf)
        p22_inf_hat = p22_inf_ + (p12_inf_ * p21_inf_) / (p11_inf_ + p22_inf)

    return (
        p11_inf_hat,
        p12_inf_hat,
        p21_inf_hat,
        p22_inf_hat,
        au_hat,
        ad_hat,
        tu_hat,
        td_hat,
        yn,
        x2,
    )


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
        self.y0 = np.linspace(self.y_nominal_down, self.y_nominal, self.n_level_nodes)
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


class GetIDZVariables:
    def __init__(
        self,
        length=0,
        h_b_up=0,
        h_b_down=0,
        q_nominal=0,
        width=0,
        y_nominal=0,
        side_slope=0,
        friction_coefficient=0,
    ):
        self.length = length
        self.h_b_up = h_b_up
        self.h_b_down = h_b_down
        self.q_nominal = q_nominal
        self.width = width
        self.y_nominal = y_nominal
        self.side_slope = side_slope
        self.friction_coefficient = friction_coefficient

    def getVariables(self):
        # s_b = (self.h_b_up - self.h_b_down) / self.length

        q = self.q_nominal
        n = self.friction_coefficient
        B = self.width
        m = self.side_slope
        Sb = (self.h_b_up - self.h_b_down) / self.length
        Y0 = self.y_nominal - self.h_b_down
        L = self.length
        (
            p11_inf_hat,
            p12_inf_hat,
            p21_inf_hat,
            p22_inf_hat,
            au_hat,
            ad_hat,
            tu_hat,
            td_hat,
            yn,
            x2,
        ) = IdzFun(q, n, B, m, Sb, Y0, L, 0)

        self.Delay_in_hour = (tu_hat + td_hat) / 2
        self.p11 = p11_inf_hat
        self.p12 = p12_inf_hat
        self.Ad = ad_hat
        self.Au = au_hat
        self.p21 = p21_inf_hat
        self.p22 = p22_inf_hat


# Test function
# res = IdzFun(2.7, 0.02, 7.0, 1.5, 0.0001, 2.1, 7000, 0)
