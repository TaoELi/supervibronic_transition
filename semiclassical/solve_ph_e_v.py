# A simple code to reveal the supervibronic transition mechanism
# System: Ne electronic two-level systems + Nv vibraitonal harmonic oscillators + a single IR cavity mode


import numpy as np
import scipy
import columnplots as clp

au2fs = 0.02418884254

# Important parameters
# Number of electronic two-level systems
Ne = 1e10
# Frequency of the electronic two-level systems
omega_e = 0.1
# Transition dipole moment for each electronic two-level system
de = 0.5
# permanent ground-state dipole moment for each electronic two-level system
mug = 0.0
# permanent excited-state dipole moment for each electronic two-level system
mue = 1.0
# spontaneous emission rate
gamma_se = 1e-5

# Number of vibraitonal molecules
Nv = 1e10
# Frequency of the vibrational mode
omega_v = 0.01
Ndark = 500
omega_v_dark = np.linspace(0.007, 0.013, Ndark)
# Transition dipole moment for each vibrational molecule
dv = 0.01
# Damping coefficient
gamma_v = 2e-6 / Ndark**0.5

# Frequency of the photon mode
omega_c = 0.01
# Coupling strength for the photon mode
lambda_c = 2e-6
# Damping coefficient
gamma_c = 2e-5

# external pulse amplitude
E0_ext = 1e-2
sigma_pulse = 100.0
sigma_pulse2 = sigma_pulse**2
t0 = 500.0


# Define the electronic Hamiltonian
Hme = np.zeros((2, 2))
Mue_t = np.zeros((2, 2))
Mue_p = np.zeros((2, 2))
for i in range(1, 2):
    Hme[i, i] = omega_e
    Mue_t[0, i] = Mue_t[i, 0] = de
    Mue_p[i, i] = mue
Mue_p[0, 0] = mug
Mue = Mue_t + Mue_p
Mue2 = np.dot(Mue, Mue)

SigmaMinus = np.array([[0, 1], [0, 0]])
SigmaPlus = np.array([[0, 0], [1, 0]])
print("Number operator is", np.dot(SigmaPlus, SigmaMinus))


print("Hme = ")
print(Hme)
print("Mue_t = ")
print(Mue_t)
print("Mue_p = ")
print(Mue_p)
print("Mue = ")
print(Mue)
print("Mue2 = ")
print(Mue2)


def wrap_initial_state(Rho_e, xc, pc, xb, pb, xd, pd):
    initial_state = np.zeros(2*2 + 4 + Ndark*2, dtype=np.complex128)
    initial_state[0:4] = Rho_e.flatten()
    initial_state[4] = xc
    initial_state[5] = pc
    initial_state[6] = xb
    initial_state[7] = pb
    initial_state[8:8+Ndark] = xd
    initial_state[8+Ndark:] = pd
    return initial_state


def decode_initial_state(initial_state):
    Rho_e = initial_state[0:4].reshape((2, 2))
    xc = initial_state[4]
    pc = initial_state[5]
    xb = initial_state[6]
    pb = initial_state[7]
    xd = initial_state[8:8+Ndark]
    pd = initial_state[8+Ndark:]
    return Rho_e, xc, pc, xb, pb, xd, pd


# Define the initial electronic state as the ground state
Rho_e = np.zeros((2, 2))
Rho_e[0, 0] = 1.0
xc, pc = 0.0, 0.0
xb, pb = 0.0, 0.0
xd, pd = np.zeros(Ndark), np.zeros(Ndark)
initial_state = wrap_initial_state(Rho_e, xc, pc, xb, pb, xd, pd)


def gaussian_pulse(t):
    return E0_ext * np.sin(omega_e * t) * np.exp(-(t - t0)**2 / sigma_pulse2)


# Define the equations of motion for the coupled photon-electronic system
def func_dynamics(t, state):
    # Decode the state to xc, pc, Rho_e
    Rho_e, xc, pc, xb, pb, xd, pd = decode_initial_state(state)
    # equations of motion
    dxcdt = pc
    current_mue_single = np.trace(np.dot(Rho_e, Mue))
    current_muv = Nv ** 0.5 * xb * dv
    dpcdt = -omega_c**2 * xc - lambda_c * omega_c * (Ne * current_mue_single + current_muv) - gamma_c * pc
    Hsc = Hme + omega_c * lambda_c * xc * Mue + 0.5 * lambda_c**2 * (Mue2 + 2.0 * Mue * (Ne-1) * current_mue_single + 2.0 * Mue * current_muv)
    # add the external pulse to excite the system
    Hsc += gaussian_pulse(t) * Mue_t
    dRhoedt = -1j * np.dot(Hsc, Rho_e) + 1j * np.dot(Rho_e, Hsc)
    # Add a Lindblad term to describe the spontaneous emission and decoherence
    Lindblad = gamma_se * (np.dot(SigmaMinus, np.dot(Rho_e, SigmaPlus)) - 0.5 * np.dot(SigmaPlus, np.dot(SigmaMinus, Rho_e)) - 0.5 * np.dot(Rho_e, np.dot(SigmaPlus, SigmaMinus)))
    dRhoedt += Lindblad

    dmudr = Nv ** 0.5 * dv
    dxbdt = pb
    dpbdt = -omega_v ** 2 * xb - lambda_c * omega_c * xc * dmudr - gamma_v * np.sum(xd) - lambda_c ** 2 * (current_muv + current_mue_single * Ne) * dmudr
    dxddt = pd
    dpddt = -omega_v_dark ** 2 * xd - gamma_v * xb

    # Encode the equations of motion
    dstate = wrap_initial_state(dRhoedt, dxcdt, dpcdt, dxbdt, dpbdt, dxddt, dpddt)
    return dstate


# Solve the ODE for the coupled photon-electronic system
solution = scipy.integrate.RK45(fun=func_dynamics, t0=0.0, y0=initial_state,
                                t_bound=1000000.0, first_step=0.04,
                                rtol=1e-6, atol=1e-9)
# collect data
t_values = []
y_values = []
for i in range(30000):
    # get solution step state
    solution.step()
    t_values.append(solution.t)
    y_values.append(solution.y)
    # break loop after modeling is finished
    if solution.status == 'finished':
        break
    # print(solution.t, np.real(solution.y))

# analyze the data
xc_lst = []
Pg_lst = []
e_ph_lst = []
e_tot_lst = []
xb_lst = []
xd_lst = []
e_vib_lst = []
e_dark_lst = []
e_He_lst = []
for t, y in zip(t_values, y_values):
    Rho_e, xc, pc, xb, pb, xd, pd = decode_initial_state(y)
    xc_lst.append(xc)
    xb_lst.append(xb)
    xd_lst.append(np.mean(xd**2))
    Pg_lst.append(np.real(Rho_e[0, 0]))
    mu_current = Ne * np.real(np.trace(np.dot(Rho_e, Mue))) + Nv ** 0.5 * xb * dv
    e_ph = 0.5 * omega_c**2 * (xc + lambda_c / omega_c * mu_current)**2 + 0.5 * pc**2
    e_ph_lst.append(e_ph)
    e_Hme = np.real(np.trace(np.dot(Rho_e, Hme))) * Ne
    e_vib = 0.5 * omega_v ** 2 * xb ** 2 + 0.5 * pb ** 2
    e_dark = np.sum(0.5 * omega_v_dark ** 2 * xd ** 2 + 0.5 * pd ** 2)
    e_dark_lst.append(e_dark)
    e_vib_lst.append(e_vib)
    e_tot_lst.append(e_ph + e_Hme + e_vib + e_dark)
    e_He_lst.append(e_Hme)


t_values = np.array(t_values)
xc_lst = np.array(xc_lst)
xb_lst = np.array(xb_lst)
xd_lst = np.array(xd_lst)
Pe_lst = 1.0 - np.array(Pg_lst)
e_ph_lst = np.array(e_ph_lst) * 219474.63
e_vib_lst = np.array(e_vib_lst) * 219474.63
e_dark_lst = np.array(e_dark_lst) * 219474.63
e_tot_lst = np.array(e_tot_lst) * 219474.63
e_He_lst = np.array(e_He_lst) * 219474.63

axes = clp.initialize(7, 1, width=3.4, height=3.4*0.618*4, sharex=False)
clp.plotone([t_values * au2fs], [xc_lst], axes[0], ylabel="xc")
clp.plotone([t_values * au2fs], [xb_lst], axes[1], ylabel="xb")
clp.plotone([t_values * au2fs], [xd_lst], axes[2], ylabel="xd")
clp.plotone([t_values * au2fs], [Pe_lst], axes[3], ylabel="Pe")
clp.plotone([t_values * au2fs], [gaussian_pulse(t_values)], axes[4], ylabel="E(t)")
#clp.plotone([t_values * au2fs]*3, [e_ph_lst, e_He_lst*1e-1, e_tot_lst*1e-1], axes[5],
#            lw=1.0, colors=["r", "y", "g"], labels=["ph", "eHe*1e-1", "eTot*1e-1"],
#            ylabel="energy [cm-1]")
clp.plotone([t_values * au2fs]*2, [e_He_lst, e_ph_lst], axes[5],
            lw=1.0, colors=["r", "y", "g"], labels=["eTot", "ePh"],
            ylabel="energy [cm-1]")

clp.plotone([t_values * au2fs]*2, [e_vib_lst, e_dark_lst, e_vib_lst + e_dark_lst], axes[6],
            lw=1.0, colors=["b", "r", "k"], labels=["bright", "dark", "tot vib"],
            xlabel="Time (fs)", ylabel="energy [cm-1]")


clp.adjust()
