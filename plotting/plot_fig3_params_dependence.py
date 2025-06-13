import numpy as np
import columnplots as clp

params = dict()
# Ne = Nv dependence
params["Ne"] = [1e8, 2e8, 3e8, 4e8, 5e8, 6e8, 8e8, 1e9, 3e9, 1e10, 2e10, 3e10, 4e10, 6e10]
params["Ne_filenames"] = ["../semiclassical/jobs/params_Ne_%.0e.out" % N for N in params["Ne"]]
params["Nv"] = [1e8, 2e8, 3e8, 4e8, 5e8, 6e8, 8e8, 1e9, 3e9, 1e10, 2e10, 3e10, 4e10, 6e10]
params["Nv_filenames"] = ["../semiclassical/jobs/params_Nv_%.0e.out" % N for N in params["Nv"]]
# omegac dependence
params["omegac"] = [0.005, 0.006, 0.007, 0.008, 0.009, 0.010, 0.011, 0.012, 0.013, 0.014, 0.015]
params["omegac_filenames"] = ["../semiclassical/jobs/params_omegac_%.3f.out" % omegac for omegac in params["omegac"]]
# coupling dependence
params["lambdac"] = [3e-8, 6e-8, 1e-7, 2e-7, 4e-7, 6e-7, 8e-7, 1e-6, 2e-6, 3e-6, 4e-6]
params["lambdac_filenames"] = ["../semiclassical/jobs/params_lambdac_%.0e.out" % lambdac for lambdac in params["lambdac"]]
# E field dependence
params["E0"] = [1e-3, 2e-3, 4e-3, 6e-3, 8e-3, 1e-2, 2e-2, 4e-2]
params["E0_filenames"] = ["../semiclassical/jobs/params_E0_%.0e.out" % E0 for E0 in params["E0"]]
# mue dependence the excited state dipole moment
params["mue"] = [0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0]
params["mue_filenames"] = ["../semiclassical/jobs/params_mue_%.1f.out" % mue for mue in params["mue"]]


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


def subplot(axes, idx=0, task="N", t_fs=5000):
    x, e_ph, e_He, e_vib, e_dark, e_tot, e_max = get_params_dependence(task=task, t_fs=t_fs)
    # we plot the dark-mode energy as a function of different parameters
    xs = [x]
    e_dark *= 1e-10
    e_max *= 1e-10
    ys = [e_dark]
    xlog = True if task == "Ne" or task == "Nv" or task == "mue" or task == "E0" or task == "lambdac" else False
    ylog = True if task == "Ne" or task == "Nv" or task == "E0" or task == "lambdac" else False
    xlabel = None
    # we plot the analytic form
    if task == "Ne":
        clp.plotone([xs[0][:-4]], [e_dark[0] * (x[:-4] / x[0]) ** 2], axes[idx], colors=[clp.navy_blue], linestyles=['--'],
                    showlegend=False, xlabel=r"$N_{\rm e}$", lw=0.8)
        axes[idx].text(0.1, 0.78, r"$\propto N_{\rm e}^2$", transform=axes[idx].transAxes,
                       fontsize=12, verticalalignment='top', color=clp.navy_blue)
        axes[idx].text(0.5, 0.2, r"inversion", transform=axes[idx].transAxes,
                       fontsize=12, verticalalignment='top', color=clp.red)
    elif task == "Nv":
        axes[idx].set_xlabel(r"$N_{\rm v}$")
        #clp.plotone(xs, [e_dark[0] * (np.ones(x[0].size) * x[0][-1])], axes[idx], colors=[clp.navy_blue], linestyles=['--'],
        #            showlegend=False, xlabel=r"$N_{\rm v}$", lw=0.8)
        #axes[idx].text(0.2, 0.7, r"$\propto N_{\rm v}$", transform=axes[idx].transAxes,
        #               fontsize=12, verticalalignment='top', color=clp.navy_blue)
    elif task == "mue":
        clp.plotone(xs, [e_dark[1] * (x / x[1])**2], axes[idx], colors=[clp.navy_blue], linestyles=['--'],
                    showlegend=False, xlabel=r"$\Delta d$ [a.u.]", lw=0.8)
        axes[idx].text(0.14, 0.7, r"$\propto \Delta d^2$", transform=axes[idx].transAxes,
                       fontsize=12, verticalalignment='top', color=clp.navy_blue)
    elif task == "omegac":
        au2cminv = 219474.63
        xs = [(x - 0.01) * au2cminv]
        xlabel = r"$\omega_{\rm c} - \omega_{\rm v}$ [cm$^{-1}$]"
    elif task == "E0":
        clp.plotone(xs, [e_dark[0] * (x / x[0])**4], axes[idx], colors=[clp.navy_blue], linestyles=['--'],
                    showlegend=False, xlabel=r"$E_0$ [a.u.]", lw=0.8)
        axes[idx].text(0.14, 0.7, r"$\propto E_0^4$", transform=axes[idx].transAxes,
                       fontsize=12, verticalalignment='top', color=clp.navy_blue)
    elif task == "lambdac":
        xs = [x * 1e7]
        clp.plotone(xs, [e_dark[1] * (x / x[1])**2], axes[idx], colors=[clp.navy_blue], linestyles=['--'],
                    showlegend=False, xlabel=r"$\lambda_{\rm c}$ [$\times 10^{-7}$ a.u.]", lw=0.8)
        axes[idx].text(0.14, 0.7, r"$\propto \lambda_{\rm c}^2$", transform=axes[idx].transAxes,
                       fontsize=12, verticalalignment='top', color=clp.navy_blue)
        axes[idx].text(0.5, 0.2, r"inversion", transform=axes[idx].transAxes,
                       fontsize=12, verticalalignment='top', color=clp.red)
    clp.plotone(xs, ys, axes[idx], colors=["ko"], xlabel=xlabel,
                # ylabel=r"$E_{\rm D}$($t = $ 5 ps) [$\times 10^{10}$ cm$^{-1}$]" if idx % 3 == 0 else None,
                xlog=xlog, ylog=ylog, showlegend=False)
    if task == "Ne":
        clp.plotone([xs[0][-4:]], [ys[0][-4:]], axes[idx], colors=["ro"], showlegend=False)
    elif task == "lambdac":
        clp.plotone([xs[0][-3:]], [ys[0][-3:]], axes[idx], colors=["ro"], showlegend=False)
    # add efficiency subplots
    '''
    ax2 = axes[idx].twinx()  # instantiate a second axes that shares the same x-axis

    lines = clp.plotone(xs, [e_dark/e_max*100], ax2, colors=["k*"], xlabel=xlabel,
                xlog=xlog, ylog=ylog, showlegend=False, yscientific=True)
    for line in lines:
        line.set_color(clp.navy_blue)
    ax2.tick_params(axis='y', labelcolor=clp.navy_blue)
    if idx == 5:
        ax2.set_ylabel(r"energy conversion percentage [%]", color=clp.navy_blue, fontsize=12)
    '''


fig, axes = clp.initialize(col=2, row=3, width=6.0, height=4.3*0.618*2, sharex=False, LaTeX=True,
                           labelthem=True, labelthemPosition=[0.25, 0.95], labelsize=14,
                           return_fig_args=True, fontsize=12, commonY=[-0.01, 0.5, r"$E_{\rm D}$($t = $ 5 ps) [$\times 10^{10}$ cm$^{-1}$]"])
t_fs = 5000

axes = axes.reshape(6)

subplot(axes, idx=0, task="lambdac", t_fs=t_fs)
subplot(axes, idx=1, task="Ne", t_fs=t_fs)
subplot(axes, idx=2, task="E0", t_fs=t_fs)
subplot(axes, idx=3, task="mue", t_fs=t_fs)
subplot(axes, idx=4, task="omegac", t_fs=t_fs)
subplot(axes, idx=5, task="Nv", t_fs=t_fs)


clp.adjust(tight_layout=True, savefile="params_dependence.pdf")

