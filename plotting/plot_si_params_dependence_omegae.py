import numpy as np
import columnplots as clp

params = dict()

# omegae dependence
params["omegae"] = [0.015, 0.020, 0.025, 0.035, 0.040, 0.050, 0.060, 0.070, 0.080, 0.100, 0.12, 0.14, 0.16, 0.20]
params["omegae_filenames"] = ["../semiclassical/jobs/params_omegae_%.3f.out" % omegae for omegae in params["omegae"]]

params["omegae2"] = [0.015, 0.020, 0.025, 0.035, 0.040, 0.050, 0.060, 0.070, 0.080, 0.100, 0.12, 0.14, 0.16, 0.20]
params["omegae2_filenames"] = ["../semiclassical/jobs/params_omegae_%.3f_mue0.out" % omegae2 for omegae2 in params["omegae2"]]

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
    if task == "omegae" or task == "omegae2":
        xs = [x/ 0.01]
        xlabel = r"$\omega_{\rm e}/\omega_{\rm v}$"
        axes[idx].axvline(x=1, color=clp.red, linestyle="--", lw=0.8)
        axes[idx].set_xlim(0, 20)
    clp.plotone(xs, ys, axes[idx], colors=["ko"], xlabel=xlabel,
                # ylabel=r"$E_{\rm D}$($t = $ 5 ps) [$\times 10^{10}$ cm$^{-1}$]" if idx % 3 == 0 else None,
                xlog=xlog, ylog=ylog, showlegend=False)

fig, axes = clp.initialize(col=1, row=2, width=8.6, height=4.3*0.618, sharex=False, LaTeX=True,
                           labelthem=True, labelthemPosition=[0.94, 0.95], labelsize=14,
                           return_fig_args=True, fontsize=12, commonY=[-0.01, 0.5, r"$E_{\rm D}$($t = $ 5 ps) [$\times 10^{10}$ cm$^{-1}$]"])

t_fs = 5000

subplot(axes, idx=0, task="omegae", t_fs=t_fs)

subplot(axes, idx=1, task="omegae2", t_fs=t_fs)

axes[0].text(0.2, 0.8, r"$\Delta d = 1.0$ a.u.", transform=axes[0].transAxes,
             fontsize=12, verticalalignment='top', color=clp.red)
axes[1].text(0.2, 0.8, r"$\Delta d = 0$", transform=axes[1].transAxes,
             fontsize=12, verticalalignment='top', color=clp.red)

clp.adjust(tight_layout=True, savefile="params_dependence_omegae.pdf")

