import numpy as np

def CircularChannelData(y, r, q):
    # Circular channel calculations
    diameter = 2 * r
    alpha = np.arccos(1 - y / r)
    area = (diameter**2) / 4 * (alpha - np.sin(2 * alpha) / 2)
    perimeter = alpha * diameter
    hydraulic_radius = (diameter / 4) * (1 - np.sin(2 * alpha) / (2 * alpha))

    # Celerity (Kinematic wave speed)
    c = np.sqrt(9.81 * hydraulic_radius)

    # Velocity
    v = q / area

    return alpha, area, perimeter, hydraulic_radius, c, v

def froudeC(y, q, n, r, m):
    g = 9.81  # Gravitational acceleration (m/s^2)

    # Circular channel calculations
    diameter = 2 * r
    alpha = np.arccos(1 - y / r)
    area = (diameter**2) / 4 * (alpha - np.sin(2 * alpha) / 2)
    perimeter = alpha * diameter
    hydraulic_radius = diameter / 4 * (1 - np.sin(2 * alpha) / (2 * alpha))

    # Celerity (Kinematic wave speed)
    c = np.sqrt(g * hydraulic_radius)

    # Velocity
    v = q / area

    # Froude number
    fr = v / c


    return fr

def froude(y, q, n, b0, m):
    g = 9.81  # Gravitational acceleration (m/s^2)

    # Trapezoidal channel calculations
    a = b0 * y + m * y**2  # Cross-sectional area
    t = b0 + 2 * y * m     # Top width
    c = np.sqrt(g * a / t) # Celerity (Kinematic wave speed)

    # Velocity
    v = q / a

    # Froude number
    fr = v / c

    return fr

def sf0C(y, q, n, r, m):
    # Circular channel calculations
    diameter = 2 * r
    alpha = np.arccos(1 - y / r)
    area = (diameter**2) / 4 * (alpha - np.sin(2 * alpha) / 2)
    perimeter = alpha * diameter
    hydraulic_radius = diameter / 4 * (1 - np.sin(2 * alpha) / (2 * alpha))

    # Friction slope
    sf0_ = (q**2 * n**2) / (area**2 * hydraulic_radius**(4 / 3))

    return sf0_


def sf0(y, q, n, b0, m):
    # Area of the trapezoidal cross-section
    area = b0 * y + m * y ** 2

    # Wetted perimeter of the trapezoidal cross-section
    perimeter = b0 + 2 * y * np.sqrt(1 + m ** 2)

    # Hydraulic radius
    hydraulic_radius = area / perimeter

    # Friction slope
    sf0_ = (q ** 2 * n ** 2) / (area ** 2 * hydraulic_radius ** (4 / 3))

    return sf0_




def IdzFun(q, n, B, m, Sb, Y0, L, shape):
    kh = 0  # initial height = 0
    sb = Sb  # bottom slope
    Man = n  # Manning coefficient??
    b0 = B  # Bottom width
    h = Y0  # initial height??
    yx = h  # downstream height

    g = 9.81

    # Normal depth
    if shape == 0:
        y1 = 0  # initial values??
        y2 = yx * 5
        for k in range(25):
            y = (y1 + y2) / 2  # average height
            a = b0 * y + m * y ** 2  # using average height to calculate area
            p = b0 + 2 * y * (1 + m ** 2) ** 0.5  # wet perimeter calculated trapezoidal cross-section
            r = a / p  # hydraulic radius
            dif = (q ** 2 * n ** 2) / (a ** 2 * r ** (4 / 3)) - sb  # partial derivative of the friction slope?? sb = bottom slope
            if dif < 0:
                y2 = (y1 + y2) / 2
            else:
                y1 = (y1 + y2) / 2
    elif shape == 1:
        y1 = 0  # initial values??
        y2 = yx * 5
        for k in range(25):
            y = (y1 + y2) / 2  # average height
            alpha, a, p, r, c, v = CircularChannelData(y, b0, q)
            dif = (q ** 2 * n ** 2) / (a ** 2 * r ** (4 / 3)) - sb  # partial derivative of the friction slope?? sb = bottom slope
            if dif < 0:
                y2 = (y1 + y2) / 2
            else:
                y1 = (y1 + y2) / 2

    yn = y  # Normal depth

    if shape == 0:
        sxx = (sb - sf0(yx, q, n, b0, m)) / (1 - (froude(yx, q, n, b0, m)) ** 2)
    elif shape == 1:
        sxx = (sb - sf0C(yx, q, n, b0, m)) / (1 - (froudeC(yx, q, n, b0, m)) ** 2)
    yu = yx - L * sxx  # good

    if sxx == 0:
        x1_calc = L
    else:
        x1_calc = max(0, L - (yx - yn) / sxx)
    x1 = x1_calc

    if x1 == 0:
        y1 = yx - sxx * L
    else:
        y1 = yn

    x2 = (L + x1) / 2

    if x1 == 0:
        y2 = yx - (L - x2) * sxx
    else:
        y2 = yn + (x2 - x1) * sxx

    sxu = 0

    # Part 2 ----- Upstream part - calculation of basic variables--------------
    if shape == 0:
        a = b0 * y1 + m * y1 ** 2
        p = b0 + 2 * y1 * (1 + m ** 2) ** 0.5
        t = b0 + 2 * y1 * m
        c = (g * a / t) ** 0.5
        v = q / a
    elif shape == 1:
        alpha, a, p, hydraulic_radius, c, v = CircularChannelData(y1, b0, q)

    cu = c
    vu = v

    sf = sf0(y1, q, n, b0, m)
    fr = froude(y1, q, n, b0, m)

    if shape == 0:
        kappa = 7 / 3 - 4 * a / (3 * t * p) * 2 * (1 + m ** 2) ** 0.5
        alpha = (t * (2 + (kappa - 1) * fr ** 2) * sb) / (a * fr * (1 - fr ** 2))
        beta = -2 * g / v * (sb - sxu)
        gamma = g * t * ((1 + kappa) * sb - (1 + kappa - fr ** 2 * (kappa - 2)) * sxu)
    elif shape == 1:
        kappa = 7 / 3 - 4 * a / (3 * hydraulic_radius * p) * 2 * (1 + m ** 2) ** 0.5
        alpha = (hydraulic_radius * (2 + (kappa - 1) * fr ** 2) * sb) / (a * fr * (1 - fr ** 2))
        beta = -2 * g / v * (sb - sxu)
        gamma = g * hydraulic_radius * ((1 + kappa) * sb - (1 + kappa - fr ** 2 * (kappa - 2)) * sxu)

    # Low frequencies
    x = x1
    td = x / (c + v)
    tu = x / (c - v)
    ad = (t ** 2 * (c ** 2 - v ** 2)) / gamma * (1 - np.exp(-gamma / (t * (c ** 2 - v ** 2)) * x))
    au = (t ** 2 * (c ** 2 - v ** 2)) / gamma * (np.exp(gamma / (t * (c ** 2 - v ** 2)) * x) - 1)

    # High frequencies
    p11_inf = 1 / (t * c * (1 - fr)) * ((1 + ((1 - fr) / (1 + fr)) ** 2 * np.exp(alpha * x)) / (1 + np.exp(alpha * x))) ** 0.5
    p12_inf = 2 / (t * c * (1 - fr ** 2)) * (np.exp(-(gamma) / (2 * t * (c ** 2 - v ** 2)) * x) / (1 + np.exp(alpha * x)) ** 0.5)
    p21_inf = 2 / (t * c * (1 - fr ** 2)) * (np.exp((gamma) / (2 * t * (c ** 2 - v ** 2)) * x) / (1 + np.exp(alpha * x)) ** 0.5)
    p22_inf = 1 / (t * c * (1 + fr)) * ((1 + ((1 + fr) / (1 - fr)) ** 2 * np.exp(alpha * x)) / (1 + np.exp(alpha * x))) ** 0.5

    # Part 3 ----Downstream part------------------------------------------------
    if shape == 0:
        a = b0 * y2 + m * y2 ** 2
        p = b0 + 2 * y2 * (1 + m ** 2) ** 0.5
        r = a / p
        t = b0 + 2 * y2 * m

        c = (g * a / t) ** 0.5
        v = q / a
        fr = v / c
        sf = (q ** 2 * n ** 2) / (a ** 2 * r ** (4 / 3))

        sxd = (sb - sf0(yx, q, n, b0, m)) / (1 - (froude(yx, q, n, b0, m)) ** 2)
    elif shape == 1:
        a = b0 * y2 + m * y2 ** 2
        p = b0 + 2 * y2 * (1 + m ** 2) ** 0.5
        r = a / p
        t = b0 + 2 * y2 * m

        c = (g * a / t) ** 0.5
        v = q / a
        fr = v / c
        sf = (q ** 2 * n ** 2) / (a ** 2 * r ** (4 / 3))

        sxd = (sb - sf0(yx, q, n, b0, m)) / (1 - (froude(yx, q, n, b0, m)) ** 2)

    dtdy = 2 * m
    dtdx = dtdy * sxd

    x = L - x1

    kappa = 7 / 3 - 4 * a / (3 * t * p) * 2 * (1 + m ** 2) ** 0.5
    alpha = t / (a * fr * (1 - fr ** 2)) * ((2 + (kappa - 1) * fr ** 2) * sb - (
                2 + (kappa - 1) * fr ** 2 - (
                    a / t ** 2 * dtdy + kappa - 2) * fr ** 4) * sxx)
    gamma1 = g * t * ((1 + kappa) * sb - (
                1 + kappa - fr ** 2 * (kappa - 2)) * sxd)  # paper
    gamma = v ** 2 * 2 * m * sxd + g * t * ((1 + kappa) * sb - (
                1 + kappa - fr ** 2 * (kappa - 2)) * sxd)  # paper tsst

    # Low frequencies
    td_ = x / (c + v)
    tu_ = x / (c - v)

    ad_ = (t ** 2 * (c ** 2 - v ** 2)) / gamma1 * (
                1 - np.exp(-gamma1 / (t * (c ** 2 - v ** 2)) * x))
    au_ = (t ** 2 * (c ** 2 - v ** 2)) / gamma1 * (
                np.exp(gamma1 / (t * (c ** 2 - v ** 2)) * x) - 1)

    # High frequencies
    p11_inf_ = 1 / (t * c * (1 - fr)) * (
                (1 + ((1 - fr) / (1 + fr)) ** 2 * np.exp(alpha * x)) / (
                    1 + np.exp(alpha * x))) ** 0.5
    p12_inf_ = 2 / (t * c * (1 - fr ** 2)) * (
                np.exp(-(gamma) / (2 * t * (c ** 2 - v ** 2)) * x) / (
                    1 + np.exp(alpha * x)) ** 0.5)
    p21_inf_ = 2 / (t * c * (1 - fr ** 2)) * (
                np.exp((gamma) / (2 * t * (c ** 2 - v ** 2)) * x) / (
                    1 + np.exp(alpha * x)) ** 0.5)
    p22_inf_ = 1 / (t * c * (1 + fr)) * (
                (1 + ((1 + fr) / (1 - fr)) ** 2 * np.exp(alpha * x)) / (
                    1 + np.exp(alpha * x))) ** 0.5

    # Part 4 -------------------Global model-----------------------------------
    if x1 == 0:
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

        p11_inf_hat = p11_inf + (p12_inf * p21_inf) / (p11_inf_ + p22_inf_)
        p12_inf_hat = (p12_inf * p12_inf_) / (p11_inf_ + p22_inf_)
        p21_inf_hat = (p21_inf * p21_inf_) / (p11_inf_ + p22_inf_)
        p22_inf_hat = p22_inf_ + (p12_inf_ * p21_inf_) / (p11_inf_ + p22_inf_)

    return p11_inf_hat, p12_inf_hat, p21_inf_hat, p22_inf_hat, au_hat, ad_hat, tu_hat, td_hat, yn, x2

q = 100.0
n = 0.02
B = 30.0
m = 0.0
Sb = 0.0002
Y0 = 3
L = 10000
shape = "Trapezoidal"
p11_inf_hat, p12_inf_hat, p21_inf_hat, p22_inf_hat, au_hat, ad_hat, tu_hat, td_hat, yn, x2 = IdzFun(q, n, B, m, Sb, Y0, L, 0)
print(p11_inf_hat, p12_inf_hat, p21_inf_hat, p22_inf_hat, au_hat, ad_hat)