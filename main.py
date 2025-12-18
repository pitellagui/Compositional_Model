import numpy as np
from barril.units import Array, Scalar
import math

P = Array('pressure', np.array([30.0]), "bar")
T = Array('temperature', np.array([220.0]), "K")
zi = Array('dimensionless', np.array([0.50, 0.50]), "-")

R = Scalar('molar gas constant', 8.314, 'J/mol.K')


class Components:
    nc = 0

    def __init__(self, name, PM, Tc, Pc, Vc, Zc, omega):
        Components.nc += 1
        self.id = Components.nc
        self.name = name
        self.PM = PM
        self.Tc = Tc
        self.Pc = Pc
        self.Vc = Vc
        self.Zc = Zc
        self.omega = omega

    def dados(self):
        print(f"Component {self.id}: {self.name}")
        print(f"  PM = {self.PM} [g/mol]")
        print(f"  Tc = {self.Tc} [K]")
        print(f"  Pc = {self.Pc} [Bar]")
        print(f"  Vc = {self.Vc} [L/kmol]")
        print(f"  Zc = {self.Zc} [-]")
        print(f"  ω = {self.omega} [-]")
        print("-" * 35)


C1 = Components("C1", 16.043, 190.560, 45.99, 98.6, 0.286, 0.011)
CO2 = Components("CO2", 44.010, 304.120, 73.74, 94.1, 0.274, 0.225)

nc = Components.nc
valor = 0.0259
kij = np.diag([valor] * nc)

PM = Array('mass per mol', np.array([C1.PM, CO2.PM]), "g/mol")
Tc = Array('temperature', np.array([C1.Tc, CO2.Tc]), "K")
Pc = Array('pressure', np.array([C1.Pc, CO2.Pc]), "bar")
Vc = Array('molar volume', np.array([C1.Vc, CO2.Vc]), "L/kmol")
Zc = Array('dimensionless', np.array([C1.Zc, CO2.Zc]), "-")
omega = Array('dimensionless', np.array([C1.omega, CO2.omega]), "-")

Ki = (Pc.GetValues('bar') / P.GetValues('bar')) * np.exp(5.37 *
                                                         (1 + omega.GetValues()) * (1 - (Tc.GetValues('K') / T.GetValues('K'))))


def Calculate_RR(zi, Ki):
    it, V_old, Erro = 0, 0.5, 100

    while Erro > 1e-6:
        fv = 0
        dfv = 0

        for i in range(nc):
            fv += zi.GetValues()[i] * (Ki[i] - 1) / (1 + V_old * (Ki[i] - 1))
            dfv += -zi.GetValues()[i] * (Ki[i] - 1)**2 / \
                (1 + V_old * (Ki[i] - 1))**2

        V_new = V_old - fv / dfv
        Erro = abs(V_new - V_old)
        V_old = V_new

    L = 1 - V_new
    V = V_new

    xi = zi.GetValues() / (1 + V * (Ki - 1))
    yi = (Ki * zi.GetValues()) / (1 + V * (Ki - 1))

    return V, L, xi, yi


def Calculate_EOS_Coef(T, P, Tc, Pc, omega, ui):
    Omega_a = 0.42748
    Omega_b = 0.08664

    Ai = np.zeros(nc)
    Bi = np.zeros(nc)

    Tri = T.GetValues('K') / Tc.GetValues('K')
    Pri = P.GetValues('bar') / Pc.GetValues('bar')

    for i in range(nc):

        m = 0.480 + 1.574 * omega[i] - 0.176 * omega[i]**2

        Ai[i] = Omega_a * (Pri[i] / (Tri[i]**2)) * \
            (1 + m * (1 - np.sqrt(Tri[i])))**2
        Bi[i] = Omega_b * (Pri[i] / Tri[i])

    Aij = np.zeros((nc, nc))

    for i in range(nc):
        for j in range(nc):
            Aij[i, j] = np.sqrt(Ai[i] * Ai[j]) * (1 - kij[i, j])

        A = np.sum(ui[:, None] * ui[None, :] * Aij)

        B = np.sum(Bi * ui)
    return A, B, Ai, Bi, Aij


def Calculate_Cubic_Equation(A, B):
    delta1 = 0
    delta2 = 1
    Delta = delta1 - delta2
    C1 = 1
    C2 = ((delta1 + delta2 - 1)*B - 1)
    C3 = (A + delta1*delta2*B**2 - (delta1 + delta2)*B*(B + 1))
    C4 = -(A*B + delta1*delta2*B**2*(B + 1))

    Z_old = 0.5
    Erro = 100
    while Erro > 1e-6:
        fz = C1*Z_old**3 + C2*Z_old**2 + C3*Z_old + C4
        dfz = 3*C1*Z_old**2 + 2*C2*Z_old + C3
        Z_new = Z_old - fz / dfz
        Erro = abs(Z_new - Z_old)
        Z_old = Z_new

    Z1 = Z_new

    a = 1
    b = (Z1 + C2)
    c = (Z1*(Z1 + C2) + C3)
    Disc = b**2 - 4*a*c

    if Disc < 0:
        Z2 = Z1
        Z3 = Z1
    elif Disc == 0:
        Z2 = Z3 = -(Z1 + C2)/2
    else:
        Z2 = (-(Z1 + C2) + math.sqrt(Disc))/2
        Z3 = (-(Z1 + C2) - math.sqrt(Disc))/2

    roots = np.array([Z1, Z2, Z3])

    Zmin = np.min(roots)
    Zmax = np.max(roots)

    ge_l = (Zmin - 1) - np.log(Zmin - B) - (A/(Delta*B)) * \
        np.log((Zmin + delta1*B)/(Zmin + delta2*B))
    ge_v = (Zmax - 1) - np.log(Zmax - B) - (A/(Delta*B)) * \
        np.log((Zmax + delta1*B)/(Zmax + delta2*B))

    Zfase = Zmax if ge_l < ge_v else Zmin

    return Zfase


def Calculate_Fugacidade(Aij, A, B, Bi, ui, P, Z):
    delta1 = 0
    delta2 = 1
    Delta = delta1 - delta2
    psi_i = np.zeros(nc)

    for i in range(nc):
        psi_i[i] = np.sum(Aij[i, :] * ui)

    phi = np.exp((Z - 1) * (Bi / B)
                 - np.log(Z - B)
                 - (A / (Delta * B)) * (2 * psi_i / A - Bi / B)
                 * np.log((Z + delta1 * B) / (Z + delta2 * B))
                 )

    fi = ui * phi * P.GetValues('bar')

    return fi


sumf = 100

while sumf > 1e-8:
    V, L, xi, yi = Calculate_RR(zi, Ki)
    A_liq, B_liq, Ai, Bi, Aij = Calculate_EOS_Coef(T, P, Tc, Pc, omega, xi)
    A_vap, B_vap, Ai, Bi, Aij = Calculate_EOS_Coef(T, P, Tc, Pc, omega, yi)
    Z_liq = Calculate_Cubic_Equation(A_liq, B_liq)
    Z_vap = Calculate_Cubic_Equation(A_vap, B_vap)
    fi_liq = Calculate_Fugacidade(Aij, A_liq, B_liq, Bi, xi, P, Z_liq)
    fi_vap = Calculate_Fugacidade(Aij, A_vap, B_vap, Bi, yi, P, Z_vap)

    Ki_new = Ki * (fi_liq / fi_vap)
    Ki = Ki_new

    sumf = np.sum((fi_liq / fi_vap - 1)**2)

Ki = yi/xi


def MixtureProperties01(ui, PM, Z, P, T, R):

    PM_medio = np.sum(ui * PM.GetValues('kg/mol'))

    rho = (PM_medio * P.GetValues('Pa')) / \
        (Z * R.GetValue() * T.GetValues('K'))

    V = PM_medio / rho

    return rho, V


def MixtureProperties02(rho_o, rho_g, Vg, Vo):
    rho_agua = 1000  # kg/m3
    rho_ar = 1.225  # kg/m3

    do = rho_o / rho_agua
    dg = rho_g / rho_ar

    API = 141.5 / do - 131.5

    RGO = Vg / Vo
    return do, dg, API, RGO


rho_g, Vg = MixtureProperties01(yi, PM, Z_vap, P, T, R)
rho_o, Vo = MixtureProperties01(xi, PM, Z_liq, P, T, R)
do, dg, API, RGO = MixtureProperties02(rho_o, rho_g, Vg, Vo)

print("Input Data:")
print(f"  T = {T}")
print(f"  P = {P}")
print(f"  zi = {zi}")
print("-" * 35)
C1.dados()
CO2.dados()
print("Output Data:")
print(f"  V = {V}")
print(f"  L = {L}")
print(f"  xi = {xi}")
print(f"  yi = {yi}")
print(f"  Z Líquido = {Z_liq}")
print(f"  Z Vapor = {Z_vap}")
print(f"  Ki = {Ki}")
print("-" * 35)
print("Thermodynamic Properties of Oil:")
print(f"  ρo = {rho_o} [kg/m3]")
print(f"  Vo = {Vo} [m3/mol]")
print(f"  do = {do} [-]")
print(f"  API = {API} [°]")
print("-" * 35)
print("Thermodynamic Properties of Gas:")
print(f"  Vg = {Vg} [m3/mol]")
print(f"  ρg = {rho_g} [kg/m3]")
print(f"  dg = {dg} [-]")
print(f"  RGO = {RGO} [m3/m3 std]")
