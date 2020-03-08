import math


def index_of(x, y, z, n):
    return ((z * n) + y) * n + x


class FluidCube:
    def __init__(self, size, diffusion, viscosity, dt):
        self.size = size
        self.dt = dt
        self.diff = diffusion
        self.visc = viscosity

        self.s = [0 for _ in range(size * size * size)]
        self.density = [0 for _ in range(size * size * size)]

        self.vx = [0 for _ in range(size * size * size)]
        self.vy = [0 for _ in range(size * size * size)]
        self.vz = [0 for _ in range(size * size * size)]

        self.vx0 = [0 for _ in range(size * size * size)]
        self.vy0 = [0 for _ in range(size * size * size)]
        self.vz0 = [0 for _ in range(size * size * size)]

    def fluid_step(self):
        n = self.size
        visc = self.visc
        diff = self.diff
        dt = self.dt
        Vx = self.vx
        Vy = self.vy
        Vz = self.vz
        Vx0 = self.vx0
        Vy0 = self.vy0
        Vz0 = self.vz0
        s = self.s
        density = self.density

        _diffuse(1, Vx0, Vx, visc, dt, 4, n)
        _diffuse(2, Vy0, Vy, visc, dt, 4, n)
        _diffuse(3, Vz0, Vz, visc, dt, 4, n)

        _project(Vx0, Vy0, Vz0, Vx, Vy, 4, n)

        _advect(1, Vx, Vx0, Vx0, Vy0, Vz0, dt, n)
        _advect(2, Vy, Vy0, Vx0, Vy0, Vz0, dt, n)
        _advect(3, Vz, Vz0, Vx0, Vy0, Vz0, dt, n)

        _project(Vx, Vy, Vz, Vx0, Vy0, 4, n)

        _diffuse(0, s, density, diff, dt, 4, n)
        _advect(0, density, s, Vx, Vy, Vz, dt, n)

    def add_density(self, x, y, z, amount):
        self.density[index_of(x, y, z, self.size)] += amount

    def add_velocity(self, x, y, z, amount_x, amount_y, amount_z):
        index = index_of(x, y, z, self.size)

        self.vx[index] += amount_x
        self.vy[index] += amount_y
        self.vz[index] += amount_z


def _set_bounds(b, x, n):
    """
    Keep the fluid from leaking out of the box by creating a mirror force
    """
    for j in range(1, n - 1):
        for i in range(1, n - 1):
            x[index_of(i, j, 0, n)] = -x[index_of(i, j, 1, n)] if b == 3 else x[index_of(i, j, 1, n)]
            x[index_of(i, j, 0, n - 1)] = -x[index_of(i, j, 1, n - 2)] if b == 3 else x[index_of(i, j, 1, n - 2)]
    for k in range(1, n - 1):
        for i in range(1, n - 1):
            x[index_of(i, 0, k, n)] = -x[index_of(i, 1, k, n)] if b == 2 else x[index_of(i, 1, k, n)]
            x[index_of(i, n - 1, 0, n - 1)] = -x[index_of(i, n - 2, k, n - 2)] if b == 2 else x[
                index_of(i, n - 2, k, n - 2)]
    for k in range(1, n - 1):
        for j in range(1, n - 1):
            x[index_of(0, j, k, n)] = -x[index_of(1, j, k, n)] if b == 1 else x[index_of(1, j, k, n)]
            x[index_of(n - 1, j, k, n - 1)] = -x[index_of(n - 2, j, k, n)] if b == 1 else x[
                index_of(n - 2, j, k, n)]

    x[index_of(0, 0, 0, n)] = 1 / 3 * (x[index_of(1, 0, 0, n)] + x[index_of(0, 1, 0, n)] + x[index_of(0, 0, 1, n)])
    x[index_of(0, n - 1, 0, n)] = 1 / 3 * (
            x[index_of(1, n - 1, 0, n)] + x[index_of(0, n - 2, 0, n)] + x[index_of(0, n - 1, 1, n)])
    x[index_of(0, 0, n - 1, n)] = 1 / 3 * (
            x[index_of(1, 0, n - 1, n)] + x[index_of(0, 1, n - 1, n)] + x[index_of(0, 0, n - 2, n)])
    x[index_of(0, n - 1, n - 1, n)] = 1 / 3 * (
            x[index_of(1, n - 1, n - 1, n)] + x[index_of(0, n - 2, n - 1, n)] + x[index_of(0, n - 1, n - 2, n)])
    x[index_of(n - 1, 0, 0, n)] = 1 / 3 * (
            x[index_of(n - 2, 0, 0, n)] + x[index_of(n - 1, 1, 0, n)] + x[index_of(n - 1, 0, 1, n)])
    x[index_of(n - 1, n - 1, 0, n)] = 1 / 3 * (
            x[index_of(n - 2, n - 1, 0, n)] + x[index_of(n - 1, n - 2, 0, n)] + x[index_of(n - 1, n - 1, 1, n)])
    x[index_of(n - 1, 0, n - 1, n)] = 1 / 3 * (
            x[index_of(n - 2, 0, n - 1, n)] + x[index_of(n - 1, 1, n - 1, n)] + x[index_of(n - 1, 0, n - 2, n)])
    x[index_of(n - 1, n - 1, n - 1, n)] = 1 / 3 * (
            x[index_of(n - 2, n - 1, n - 1, n)] + x[index_of(n - 1, n - 2, n - 1, n)] + x[
        index_of(n - 1, n - 1, n - 2, n)])


def _lin_solve(b, x, x0, a, c, iterations, n):
    """
    Solve linear differential equation
    """
    c_recip = 1 / c
    for k in range(0, iterations):
        for m in range(1, n - 1):
            for j in range(1, n - 1):
                for i in range(1, n - 1):
                    x[index_of(i, j, m, n)] = (x0[index_of(i, j, m, n)] + a * (x[index_of(i + 1, j, m, n)]
                                                                               + x[index_of(i - 1, j, m, n)]
                                                                               + x[index_of(i, j + 1, m, n)]
                                                                               + x[index_of(i, j - 1, m, n)]
                                                                               + x[index_of(i, j, m + 1, n)]
                                                                               + x[index_of(i, j, m - 1, n)]
                                                                               )) * c_recip
        _set_bounds(b, x, n)


def _diffuse(b, x, x0, diff, dt, iterations, n):
    """
    Diffusion
    """
    a = dt * diff * (n - 2) * (n - 2)
    _lin_solve(b, x, x0, a, 1 + 6 * a, iterations, n)


def _advect(b, d, d0, vx, vy, vz, dt, n):
    """
    Apply velocity
    """
    dtx = dt * (n - 2)
    dty = dt * (n - 2)
    dtz = dt * (n - 2)

    for k in range(1, n - 1):
        for j in range(1, n - 1):
            for i in range(1, n - 1):
                tmp1 = dtx * vx[index_of(i, j, k, n)]
                tmp2 = dty * vy[index_of(i, j, k, n)]
                tmp3 = dtz * vz[index_of(i, j, k, n)]
                x = i - tmp1
                y = j - tmp2
                z = k - tmp3

                if x < 0.5:
                    x = 0.5
                elif x > n + 0.5:
                    x = n + 0.5
                i0 = math.floor(x)
                i1 = i0 + 1

                if y < 0.5:
                    y = 0.5
                elif y > n + 0.5:
                    y = n + 0.5
                j0 = math.floor(y)
                j1 = j0 + 1.0

                if z < 0.5:
                    z = 0.5
                elif z > n + 0.5:
                    z = n + 0.5
                k0 = math.floor(z)
                k1 = k0 + 1

                s1 = x - i0
                s0 = 1 - s1
                t1 = y - j0
                t0 = 1 - t1
                u1 = z - k0
                u0 = 1 - u1

                i0i = int(i0)
                i1i = int(i1)
                j0i = int(j0)
                j1i = int(j1)
                k0i = int(k0)
                k1i = int(k1)

                d[index_of(i, j, k, n)] = s0 * (
                        t0 * (u0 * d0[index_of(i0i, j0i, k0i, n)] + u1 * d0[index_of(i0i, j0i, k1i, n)]) + (t1 * (
                        u0 * d0[index_of(i0i, j1i, k0i, n)] + u1 * d0[index_of(i0i, j1i, k1i, n)]))) + s1 * (
                                                  t0 * (u0 * d0[index_of(i1i, j0i, k0i, n)] + u1 * d0[
                                              index_of(i1i, j0i, k1i, n)]) + (t1 * (
                                                  u0 * d0[index_of(i1i, j1i, j0i, n)] + u1 * d0[
                                              index_of(i1i, j1i, k1i, n)])))
    _set_bounds(b, d, n)


def _project(vx, vy, vz, p, div, iterations, n):
    """
    Normalize fluid quantity
    """
    for k in range(1, n - 1):
        for j in range(1, n - 1):
            for i in range(1, n - 1):
                div[index_of(i, j, k, n)] = -0.5 * (
                        vx[index_of(i + 1, j, k, n)]
                        - vx[index_of(i - 1, j, k, n)]
                        + vy[index_of(i, j + 1, k, n)]
                        - vy[index_of(i, j - 1, k, n)]
                        + vz[index_of(i, j, k + 1, n)]
                        - vz[index_of(i, j, k - 1, n)]
                ) / n
                p[index_of(i, j, k, n)] = 0
    _set_bounds(0, div, n)
    _set_bounds(0, p, n)
    _lin_solve(0, p, div, 1, 6, iterations, n)

    for k in range(1, n - 1):
        for j in range(1, n - 1):
            for i in range(1, n - 1):
                vx[index_of(i, j, k, n)] -= 0.5 * (p[index_of(i + 1, j, k, n)] - p[index_of(i - 1, j, k, n)]) * n
                vy[index_of(i, j, k, n)] -= 0.5 * (p[index_of(i, j + 1, k, n)] - p[index_of(i, j - 1, k, n)]) * n
                vz[index_of(i, j, k, n)] -= 0.5 * (p[index_of(i, j, k + 1, n)] - p[index_of(i, j, k - 1, n)]) * n
    _set_bounds(1, vx, n)
    _set_bounds(2, vy, n)
    _set_bounds(3, vz, n)
