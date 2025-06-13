import numpy as np
import scipy
import columnplots as clp

au2fs = 0.02418884254

# Important parameters
# Number of electronic two-level systems
Ne = 1e6
# Frequency of the electronic two-level systems
omega_e = 0.2
# Transition dipole moment for each electronic two-level system
de = 0.5 * 0.0
# permanent ground-state dipole moment for each electronic two-level system
mug = 0.0
# permanent excited-state dipole moment for each electronic two-level system
mue = 0.2 * 0.0

# Number of vibraitonal molecules
Nv = 1e10
# Frequency of the vibrational mode
omega_v = 0.01
Ndark = 300
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
gamma_c = 2e-5 * 0.0

# external pulse amplitude
E0_ext = 1e-3
sigma_pulse = 200.0
sigma_pulse2 = sigma_pulse**2
t0 = 5e4


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

print("Hme = ")
print(Hme)
print("Mue_t = ")
print(Mue_t)
print("Mue_p = ")
print(Mue_p)
print("Mue2 = ")
print(Mue2)


def wrap_initial_state(Rho_e, xc, pc):
    initial_state = np.zeros(2*2 + 2, dtype=np.complex128)
    initial_state[0] = xc
    initial_state[1] = pc
    initial_state[2:] = Rho_e.flatten()
    return initial_state


def decode_initial_state(initial_state):
    xc = initial_state[0]
    pc = initial_state[1]
    Rho_e = initial_state[2:].reshape((2, 2))
    return xc, pc, Rho_e


def wrap_initial_state_phv(xc, pc, xb, pb, xd, pd):
    initial_state = np.array([xc, pc, xb, pb] + xd.tolist() + pd.tolist())
    return initial_state


def decode_initial_state_phv(initial_state):
    xc = initial_state[0]
    pc = initial_state[1]
    xb = initial_state[2]
    pb = initial_state[3]
    xd = initial_state[4:4+Ndark]
    pd = initial_state[4+Ndark:]
    return xc, pc, xb, pb, xd, pd


# Define the initial state of the system
xc, pc = 0.0, 100.0
xb, pb = 0.0, 0.0
xd, pd = 0.0, 0.0
xd, pd = np.zeros(Ndark), np.zeros(Ndark)
initial_state = wrap_initial_state_phv(xc, pc, xb, pb, xd, pd)


def gaussian_pulse(t):
    return E0_ext * np.sin(omega_e * t) * np.exp(-(t - t0)**2 / sigma_pulse2)


# Define the equations of motion for the coupled photon-electronic system
def func_dynamics(t, state):
    # Decode the state to xc, pc, Rho_e
    xc, pc, xb, pb, xd, pd = decode_initial_state_phv(state)
    # equations of motion
    dxcdt = pc
    current_mu = Nv**0.5 * xb * dv
    dmudr = Nv**0.5 * dv
    dpcdt = -omega_c**2 * xc - lambda_c * omega_c * current_mu - gamma_c * pc
    dxbdt = pb
    dpbdt = -omega_v**2 * xb - lambda_c * omega_c * xc * dmudr - lambda_c**2 * current_mu * dmudr - gamma_v * np.sum(xd)
    dxddt = pd
    dpddt = -omega_v_dark**2 * xd - gamma_v * xb
    # we need to assign energy transfer to let bright->dark state energy flow
    # Encode the equations of motion
    dstate = wrap_initial_state_phv(dxcdt, dpcdt, dxbdt, dpbdt, dxddt, dpddt)
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
xb_lst = []
xd_lst = []
e_ph_lst = []
e_vib_lst = []
e_dark_lst = []
for t, y in zip(t_values, y_values):
    xc, pc, xb, pb, xd, pd = decode_initial_state_phv(y)
    xc_lst.append(xc)
    xb_lst.append(xb)
    xd_lst.append(np.sum(xd))
    e_ph = 0.5 * omega_c**2 * (xc + lambda_c / omega_c * Nv**0.5 * dv * xb)**2 + 0.5 * pc**2
    e_vib = 0.5 * omega_v**2 * xb**2 + 0.5 * pb**2
    e_dark = np.sum(0.5 * omega_v**2 * xd**2 + 0.5 * pd**2)
    e_ph_lst.append(e_ph)
    e_dark_lst.append(e_dark)
    e_vib_lst.append(e_vib)

t_values = np.array(t_values)
xc_lst = np.array(xc_lst)
xb_lst = np.array(xb_lst)
xd_lst = np.array(xd_lst)
e_ph_lst = np.array(e_ph_lst) * 219474.63
e_vib_lst = np.array(e_vib_lst) * 219474.63
e_dark_lst = np.array(e_dark_lst) * 219474.63
e_tot_lst = e_ph_lst + e_vib_lst + e_dark_lst

axes = clp.initialize(7, 1, width=3.4, height=3.4*0.618*4, sharex=True)
clp.plotone([t_values * au2fs], [xc_lst], axes[0], ylabel="xc", xlim=[0, 5000])
clp.plotone([t_values * au2fs], [xb_lst], axes[1], ylabel="xb")
clp.plotone([t_values * au2fs], [xd_lst], axes[2], ylabel="xd")
#clp.plotone([t_values * au2fs], [Pe_lst], axes[3], ylabel="Pe")
clp.plotone([t_values * au2fs], [gaussian_pulse(t_values)], axes[4], ylabel="E(t)")
clp.plotone([t_values * au2fs]*2, [e_ph_lst, e_tot_lst*1e-2], axes[5],
            lw=1.0, colors=["r", "g"], labels=["ph", "eTot*1e-2"],
            ylabel="ph energy [cm-1]")
clp.plotone([t_values * au2fs]*2, [e_vib_lst, e_dark_lst], axes[6],
            lw=1.0, colors=["b", "k"], labels=["vib", "dark"],
            xlabel="Time (fs)", ylabel="ph energy [cm-1]")
clp.adjust()