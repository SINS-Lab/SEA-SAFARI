import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import sys

def dlb_exp(r, a, b, c, d):
    return a * np.exp(-b * r) + c * np.exp(-d * r)

for i in range(1, len(sys.argv)):
    dir = sys.argv[i]
    far = np.loadtxt(f'./{dir}/far.log').transpose()
    near = np.loadtxt(f'./{dir}/near.log').transpose()

    ht_to_ev = 27.2114

    x_f = far[0]
    y_f = far[1]
    e_0 = y_f[-1]
    y_f -= e_0
    y_f *= ht_to_ev
    x_n = near[0]
    y_n = near[1] - e_0
    y_n *= ht_to_ev

    x_n = np.flip(x_n)
    y_n = np.flip(y_n)

    extras_x = x_f[x_f>x_n[-1]]
    extras_y = y_f[x_f>x_n[-1]]

    x_int = np.hstack((x_n, extras_x))
    y_int = np.hstack((y_n, extras_y))

    r_max = 10
    r_step = 0.002
    r_min = r_step
    x_arr = np.arange(r_min, r_max, r_step)
    from scipy.interpolate import CubicSpline as spline
    sp = spline(x_int, y_int)
    y_arr = sp(x_arr)

    y_r_max = sp(r_max)
    print(y_r_max)

    y_arr -= y_r_max
    y_n -= y_r_max
    y_f -= y_r_max

    plt.plot(x_arr, y_arr, label=f"{dir} Interpolated")
    plt.scatter(x_n, y_n, label=f"{dir} Way in", marker='o')
    plt.scatter(x_f, y_f, label=f"{dir} Way out", marker='x')

    V_r = y_arr
    dV_dr = -np.gradient(y_arr)
    args = dir.split("_")
    A = args[0]
    B = args[1]
    pots_file = open(dir+'.pots', 'w')
    header = '# Generated {} Potentials for {}-{}.\n#\n'+\
            '# Table Range in file (min, step, max): {}-{}-{}\n#\n'+\
            '# Atom1\tAtom2\tV_r\t-dV_dr\n'
    header = header.format("PySCF", A, B, r_min, r_step, r_max)
    pots_file.write(header)
    fmt = '{}\t{}\t{:.5e}\t{:.5e}\n'
    for i in range(len(V_r)):
        pots_file.write(fmt.format(A, B, V_r[i], dV_dr[i]))
    pots_file.close()

    popt, _ = curve_fit(dlb_exp, x_arr, y_arr)
    print(popt)

    

plt.xlabel("Distance (Å)")
plt.ylabel("Energy (eV)")
plt.legend()
plt.show()