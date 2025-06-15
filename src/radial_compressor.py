from math import pi, radians, degrees
from numpy import sqrt, sin, cos, tan, array


cot = lambda x: 1 / tan(x)


class RadialCompressor:
    """Центробежный компрессор"""

    def __init__(self, D_, alpha1, betta2l, d1_):
        self.D_ = D_
        self.alpha1 = alpha1
        self.betta2l = betta2l
        self.d1_ = d1_

    def fit(self, rotation_velocity):
        while True:
            c2m_ = c2m0_ * sin(betta2l)  # = c2m/u2
            z_rotor = z_rotor_0 * sin(betta2l)
            r1_av_ = sqrt(0.5 * (1 + d1_**2))
            sigma_t = 1 - sqrt(sin(betta2l)) * z_rotor ** (-0.7)
            mu_betta = (sigma_t - c2m0_ * cos(betta2l)) / (1 - c2m0_ * cos(betta2l))
            h_k_ = (
                mu_betta
                + alpha_fr
                - c2m_ * (mu_betta / tan(betta2l) + D_ * r1_av_ / tan(alpha1) / km)
            )

            Cp = 1006
            q = 0
            effeff = (pipi ** ((k - 1) / k) - 1) / (
                pipi ** ((k - 1) / k / eff_n) - 1 + q / Cp / TT_inlet
            )
            Hks = Cp * TT_inlet * (pipi ** ((k - 1) / k) - 1)
            Hks_ = h_k_ * effeff
            u2 = sqrt(Hks / Hks_)
            c2m = c2m_ * u2
            c1m = c2m / km
            c1 = c1m / sin(alpha1)
            c1a = c1m * cos(gamma1)
            c1u = c1m / tan(alpha1)
            a_cr1 = sqrt(2 * k / (k + 1) * R * TT_inlet)
            lamda1m = c1m / a_cr1
            lambda1 = c1 / a_cr1

            sigma_inlet = 1 / (
                1 + dzetta_inlet * (k / (k + 1) * GDF("ro", lambda1) * lambda1**2)
            )

            PP_1 = PP_inlet * sigma_inlet
            P_1 = PP_1 * GDF("P", lambda1)
            TT_1 = TT_inlet
            T_1 = TT_1 * GDF("T", lambda1)
            roro_1 = PP_1 / (R * TT_1)
            ro_1 = roro_1 * GDF("ro", lambda1)

            F1a = G / (ro_1 * c1a * kG1)
            D1 = [0] * 3
            D1[2] = sqrt((4 * F1a) / (pi * (1 - d1_**2)))
            D1[1] = D1[2] * r1_av_
            D1[0] = D1[2] * d1_
            h1 = 0.5 * (D1[2] - D1[0])
            D2 = u2 / pi / rotation_velocity

            break

    def plot(self):
        pass


def GDF(parameter: str, l):
    if parameter == "T":
        return 1 - ((k - 1) / (k + 1)) * l**2
    if parameter == "P":
        return GDF("T", l) ** (k / (k - 1))
    else:
        return GDF("T", l) ** (1 / (k - 1))


def main():

    # yellow parameter
    c2m0_ = 0.3  # 0.2..0.45
    z_rotor_0 = 20  # 16..38
    alpha_fr = 0.04  # 0.03..0.06
    eff_n = 0.86  # 0.83..0.87
    km = 1  # 0.9..1.1
    gamma1 = radians(15)  # 0..30
    dzetta_i = 0
    dzetta_na = 0.01 + 0.5 * ((90 - degrees(alpha1)) / 100) ** 2
    dzetta_inlet = dzetta_i * sin(alpha1) ** 2 + dzetta_na
    kG1 = 0.98  # 0.96..0.99 коэф. загроможденности въодного устройства

    D_ = 0.05  # 0.45..0.65
    alpha1 = radians(70)  # 60..90
    betta2l = radians(70)  # 50..90
    d1_ = 0.4  # 0.35..0.6
    c = RadialCompressor(D_, alpha1, betta2l, d1_)

    TT_inlet = 300
    PP_inlet = 10**5
    G = 3
    k = 1.4
    R = 287
    pipi = 3
    rotation_velocity = 12_000
    c.fit(rotation_velocity)


if __name__ == "__main__":
    main()
