import numpy as np
import columnplots as clp

params = dict()

# omegae dependence
params["angle"] = np.linspace(0, np.pi, 10)
params["angle_filenames"] = ["../semiclassical/jobs_fulldm/params_angle_%.3f.out" % angle for angle in params["angle"]]


def import_data(filename):
    data = np.loadtxt(filename)
    t, Pe, xc, xb, xd, eHe, ePh, eVib, eDark, eTot = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6], data[:,7], data[:,8], data[:,9]
    data_dict = {"t": t, "Pe": Pe, "xc": xc, "xb": xb, "xd": xd, "eHe": eHe, "ePh": ePh, "eVib": eVib, "eDark": eDark, "eTot": eTot}
    return data_dict


def get_long_time_energy(filename, t_fs=5000):
    data = import_data(filename)
    t_values = data["t"]
    Pe_lst = data["Pe"]
    xc_lst = data["xc"]
    xb_lst = data["xb"]
    xd_lst = data["xd"]
    e_He_lst = data["eHe"]
    e_ph_lst = data["ePh"]
    e_vib_lst = data["eVib"]
    e_dark_lst = data["eDark"]
    e_tot_lst = data["eTot"]
    idx = np.argmin(np.abs(t_values - t_fs))
    print(filename, "t = ", t_values[idx], "fs")
    return e_ph_lst[idx], e_He_lst[idx], e_vib_lst[idx], e_dark_lst[idx], e_tot_lst[idx], np.max(e_tot_lst)


def get_params_dependence(task="N", t_fs=5000):
    filenames = params[task + "_filenames"]
    e_ph_lst, e_He_lst, e_vib_lst, e_dark_lst, e_tot_lst, e_max_lst = [], [], [], [], [], []
    for filename in filenames:
        e_ph, e_He, e_vib, e_dark, e_tot, e_max = get_long_time_energy(filename, t_fs=t_fs)
        e_ph_lst.append(e_ph)
        e_He_lst.append(e_He)
        e_vib_lst.append(e_vib)
        e_dark_lst.append(e_dark)
        e_tot_lst.append(e_tot)
        e_max_lst.append(e_max)
    e_ph = np.array(e_ph_lst)
    e_He = np.array(e_He_lst)
    e_vib = np.array(e_vib_lst)
    e_dark = np.array(e_dark_lst)
    e_tot = np.array(e_tot_lst)
    e_max = np.array(e_max_lst)
    x = np.array(params[task])
    return x, e_ph, e_He, e_vib, e_dark, e_tot, e_max


def get_avg_angle_func(angle=np.pi, Ne=100):
    angle_lst = np.linspace(-angle, angle, Ne)
    r = np.mean(np.cos(angle_lst))**2
    return r


def subplot(axes, idx=0, task="N", t_fs=5000, Ne=100):
    x, e_ph, e_He, e_vib, e_dark, e_tot, e_max = get_params_dependence(task=task, t_fs=t_fs)
    # we plot the dark-mode energy as a function of different parameters
    xs = [x / np.pi]
    e_dark /= Ne
    e_max /= Ne
    ys = [e_dark]
    xlog = False #if task == "Ne" or task == "Nv" or task == "mue" or task == "E0" or task == "lambdac" else False
    ylog = False #if task == "Ne" or task == "Nv" or task == "E0" or task == "lambdac" else False
    xlabel = r"dipole angle threshold for e-TLSs $\theta_{\rm m}$ [$\pi$]" if task == "angle" else None

    # plot the analytical trend
    ys_analytical = np.array([get_avg_angle_func(angle=x_ele) for x_ele in x]) * np.max(e_dark)
    clp.plotone([x / np.pi], [ys_analytical], axes[idx], colors=[clp.navy_blue], linestyles=['--'],
                showlegend=False, xlabel=r"$N_{\rm e}$", lw=0.8)

    axes[idx].text(0.1, 0.7, r"$\propto \left\langle\cos\theta\right\rangle^2$", transform=axes[idx].transAxes,
                   fontsize=12, verticalalignment='top', color=clp.navy_blue)

    clp.plotone(xs, ys, axes[idx], colors=["ko"], xlabel=xlabel,
                ylabel=r"$E_{\rm D}$($t = $ 5 ps) / $N_{\rm e}$ [cm$^{-1}$]",
                xlog=xlog, ylog=ylog, showlegend=False)




fig, ax = clp.initialize(col=1, row=1, width=4.3, height=4.3*0.618, sharex=False, LaTeX=True,
                         return_fig_args=True, fontsize=12)

axes = [ax]

t_fs = 5000

subplot(axes, idx=0, task="angle", t_fs=t_fs)

clp.adjust(tight_layout=True, savefile="params_dependence_eangle.pdf")

