import numpy as np
import scipy
import columnplots as clp

au2fs = 0.02418884254

# Important parameters
# Number of electronic two-level systems
Ne = 1e9
# Frequency of the electronic two-level systems
omega_e = 0.1
# Transition dipole moment for each electronic two-level system
de = 0.5
# permanent ground-state dipole moment for each electronic two-level system
mug = 0.0
# permanent excited-state dipole moment for each electronic two-level system
mue = 0.2
# spontaneous emission rate
gamma_se = 1e-5

# Frequency of the photon mode
omega_c = 0.01
# Coupling strength for the photon mode
lambda_c = 2e-6
# Damping coefficient
gamma_c = 2e-5

# external pulse amplitude
E0_ext = 1e-3
sigma_pulse = 1000.0
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

SigmaMinus = np.array([[0, 1], [0, 0]])
SigmaPlus = np.array([[0, 0], [1, 0]])
print("Number operator is", np.dot(SigmaPlus, SigmaMinus))


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


# Define the initial electronic state as the ground state
Rho_e = np.zeros((2, 2))
Rho_e[0, 0] = 1.0
xc, pc = 0.0, 0.0
initial_state = wrap_initial_state(Rho_e, xc, pc)


def gaussian_pulse(t):
    return E0_ext * np.sin(omega_e * t) * np.exp(-(t - t0)**2 / sigma_pulse2)


# Define the equations of motion for the coupled photon-electronic system
def func_dynamics(t, state):
    # Decode the state to xc, pc, Rho_e
    xc, pc, Rho_e = decode_initial_state(state)
    # equations of motion
    dxcdt = pc
    current_mu = np.trace(np.dot(Rho_e, Mue))
    dpcdt = -omega_c**2 * xc - lambda_c * omega_c * Ne * current_mu - gamma_c * pc
    Hsc = Hme + omega_c * lambda_c * xc * Mue + 0.5 * lambda_c**2 * (Mue2 + 2.0 * Mue * (Ne-1) * current_mu)
    # add electronic relaxation term using the nonlinear Hamiltonian
    V_dephasing = (1j * gamma_se / 2.0) * (np.trace(np.dot(SigmaPlus, Rho_e)) * SigmaMinus - np.trace(np.dot(SigmaMinus, Rho_e)) * SigmaPlus)
    Hsc += V_dephasing
    # add the external pulse to excite the system
    Hsc += gaussian_pulse(t) * Mue_t
    dRhoedt = -1j * np.dot(Hsc, Rho_e) + 1j * np.dot(Rho_e, Hsc)
    # Add a Lindblad term to describe the spontaneous emission and decoherence
    #Lindblad = gamma_se * (np.dot(SigmaMinus, np.dot(Rho_e, SigmaPlus)) - 0.5 * np.dot(SigmaPlus, np.dot(SigmaMinus, Rho_e)) - 0.5 * np.dot(Rho_e, np.dot(SigmaPlus, SigmaMinus)))
    #dRhoedt += Lindblad
    # Encode the equations of motion
    dstate = wrap_initial_state(dRhoedt, dxcdt, dpcdt)
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
for t, y in zip(t_values, y_values):
    xc, pc, Rho_e = decode_initial_state(y)
    xc_lst.append(xc)
    Pg_lst.append(np.real(Rho_e[0, 0]))
    e_ph = 0.5 * omega_c**2 * (xc + lambda_c / omega_c * Ne * np.real(np.trace(np.dot(Rho_e, Mue))))**2 + 0.5 * pc**2
    e_ph_lst.append(e_ph)
    e_Hme = np.real(np.trace(np.dot(Rho_e, Hme))) * Ne
    e_tot_lst.append(e_ph + e_Hme)


t_values = np.array(t_values)
xc_lst = np.array(xc_lst)
Pe_lst = 1.0 - np.array(Pg_lst)
e_ph_lst = np.array(e_ph_lst)
e_tot_lst = np.array(e_tot_lst)

axes = clp.initialize(4, 1, width=4, height=4*0.618*3, sharex=True)
clp.plotone([t_values * au2fs], [xc_lst], axes[0], ylabel="xc")
clp.plotone([t_values * au2fs], [Pe_lst], axes[1], ylabel="Pe")
clp.plotone([t_values * au2fs], [gaussian_pulse(t_values)], axes[2], ylabel="E(t)")
clp.plotone([t_values * au2fs]*2, [e_ph_lst * 219474.63, e_tot_lst * 219474.63], axes[3],
            colors=["r", 'b'], lw=0.5,
            xlabel="Time (fs)", ylabel="ph energy [cm-1]")
clp.adjust()