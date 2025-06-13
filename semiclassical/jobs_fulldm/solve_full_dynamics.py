# A simple code to reveal the supervibronic transition mechanism
# System: Ne electronic two-level systems + Nv vibraitonal harmonic oscillators + a single IR cavity mode

import numpy as np
import scipy
import columnplots as clp
import json
import argparse
from tqdm import tqdm

# input parameters from command line
parser = argparse.ArgumentParser(description="Simulate the full dynamics of the coupled e-v-ph system under VSC",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-i", "--input", help="input filename of the json parameters", default="none")
parser.add_argument("-p", "--plot",  help="boolean to control if plot the simulated results", default=False)
parser.add_argument("-s", "--save",  help="output filename of the simulation data", default="none")
args = parser.parse_args()
config = vars(args)
print(config)

if config["input"] == "none":
    print("Please specify the input file name such as params_template.json")
    exit()
else:
    input_file = config["input"]
    print("Input file to control the parameters is", input_file)

au2fs = 0.02418884254

# Load the parameters from the json file
with open(input_file, "r") as handle:
    params = json.load(handle)

# Important parameters
# Number of electronic two-level systems
Ne = params["Ne"]
# Frequency of the electronic two-level systems
omega_e = params["omega_e"]
# Transition dipole moment for each electronic two-level system
de = params["de"]
# permanent ground-state dipole moment for each electronic two-level system
mug = params["mug"]
# permanent excited-state dipole moment for each electronic two-level system
mue = params["mue"]
# spontaneous emission rate
gamma_se = params["gamma_se"]

# Number of vibraitonal molecules
Nv = params["Nv"]
# Frequency of the vibrational mode
omega_v = params["omega_v"]
Ndark = params["Ndark"]
omega_dstart = params["omega_dstart"]
omega_dend = params["omega_dend"]
omega_v_dark = np.linspace(omega_dstart, omega_dend, Ndark)
# Transition dipole moment for each vibrational molecule
dv = params["dv"]
# Damping coefficient
gamma_v = params["gamma_v"] / Ndark**0.5

# Frequency of the photon mode
omega_c = params["omega_c"]
# Coupling strength for the photon mode
lambda_c = params["lambda_c"]
# Damping coefficient
gamma_c = params["gamma_c"]

# external pulse amplitude
E0_ext = params["E0_ext"]
sigma_pulse = params["sigma_pulse"]
sigma_pulse2 = sigma_pulse**2
t0 = params["t0"]


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

# print("Number operator is", np.dot(SigmaPlus, SigmaMinus))
#print("Hme = ")
#print(Hme)
#print("Mue_t = ")
#print(Mue_t)
#print("Mue_p = ")
#print(Mue_p)
#print("Mue2 = ")
#print(Mue2)


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
Rho_e[0, 0] = 1.0 - params["Pe_start"]
Rho_e[1, 1] = params["Pe_start"]
Rho_e[0, 1] = np.sqrt(Rho_e[0, 0] * Rho_e[1, 1])
Rho_e[1, 0] = Rho_e[0, 1]

xc, pc = params["xc_start"], 0.0
xb, pb = 0.0, 0.0
xd, pd = np.zeros(Ndark), np.zeros(Ndark)
initial_state = wrap_initial_state(Rho_e, xc, pc, xb, pb, xd, pd)

print("Initial state is")
print("Rho_e = ", Rho_e)
print("xc = ", xc)
print("pc = ", pc)
print("xb = ", xb)
print("pb = ", pb)
print("xd = zeros")
print("pd = zeros")


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
                                first_step=0.01,
                                t_bound=10000000.0,
                                rtol=1e-8, atol=1e-10)
                                #rtol=1e-4, atol=1e-7)
# collect data
t_values = []
y_values = []
nstep_max = params['nstep_max']
#for i in range(nstep_max):
# showing the progress bar to make the progress more fancy
for i in tqdm(range(nstep_max)):
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
xc_lst = np.real(np.array(xc_lst))
xb_lst = np.real(np.array(xb_lst))
xd_lst = np.real(np.array(xd_lst))
Pe_lst = 1.0 - np.array(Pg_lst)
e_ph_lst = np.real(np.array(e_ph_lst)) * 219474.63
e_vib_lst = np.real(np.array(e_vib_lst)) * 219474.63
e_dark_lst = np.real(np.array(e_dark_lst)) * 219474.63
e_tot_lst = np.real(np.array(e_tot_lst)) * 219474.63
e_He_lst = np.real(np.array(e_He_lst)) * 219474.63


# plot the results
if config['plot'] == 'True' or config['plot'] == 'true' or config['plot'] == '1':
    axes = clp.initialize(7, 1, width=3.4, height=3.4*0.618*4, sharex=True)
    clp.plotone([t_values * au2fs], [xc_lst], axes[0], ylabel="xc", xlim=[0, 5000], showlegend=False)
    clp.plotone([t_values * au2fs], [xb_lst], axes[1], ylabel="xb", showlegend=False)
    clp.plotone([t_values * au2fs], [xd_lst], axes[2], ylabel="xd", showlegend=False)
    clp.plotone([t_values * au2fs], [Pe_lst], axes[3], ylabel="Pe", showlegend=False)
    clp.plotone([t_values * au2fs], [gaussian_pulse(t_values)], axes[4], ylabel="E(t)", showlegend=False)
    clp.plotone([t_values * au2fs]*3, [e_ph_lst, e_He_lst*1e-1, e_tot_lst*1e-1], axes[5],
            lw=1.0, colors=["r", "y", "g"], labels=["ph", "eHe*1e-1", "eTot*1e-1"],
            ylabel="energy [cm-1]")
    clp.plotone([t_values * au2fs]*2, [e_vib_lst, e_dark_lst], axes[6],
            lw=1.0, colors=["b", "r"], labels=["vib", "dark"],
            xlabel="Time (fs)", ylabel="energy [cm-1]")
    clp.adjust()

# save the simulation data
if config['save'] != "none":
    data_final = np.zeros((len(t_values), 10))
    data_final[:, 0] = t_values * au2fs
    data_final[:, 1] = np.real(Pe_lst)
    data_final[:, 2] = np.real(xc_lst)
    data_final[:, 3] = np.real(xb_lst)
    data_final[:, 4] = np.real(xd_lst)
    data_final[:, 5] = np.real(e_He_lst)
    data_final[:, 6] = np.real(e_ph_lst)
    data_final[:, 7] = np.real(e_vib_lst)
    data_final[:, 8] = np.real(e_dark_lst)
    data_final[:, 9] = np.real(e_tot_lst)
    header = "# t[fs], Pe [a.u.], xc [a.u.], xb [a.u.], mean(xd**2) [a.u.], eHe [cm-1], ePh [cm-1], eVib [cm-1], eDark [cm-1], eTot [cm-1]"
    np.savetxt(config['save'], data_final, fmt='%10.5E', header=header)
    print("Data saved to %s" % config['save'])


