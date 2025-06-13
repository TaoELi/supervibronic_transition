import numpy as np
import columnplots as clp

au2ev = 27.211399
ev2ps = 4.13567 * 1e-3

params = dict()
# electronic depopulation dependence
params["gammae"] = [1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2]
params["gammae_filenames"] = ["../semiclassical/jobs/params_gammae_%.0e.out" % gammae for gammae in params["gammae"]]
params["gammae_label"] = r"$1/\gamma_{\rm e}$ [ps]"
# cavity loss dependence
params["gammac"] = [1e-6, 3e-6, 1e-5, 2e-5, 5e-5, 1e-4, 3e-4, 1e-3]
params["gammac_filenames"] = ["../semiclassical/jobs/params_gammac_%.0e.out" % gammac for gammac in params["gammac"]]
params["gammac_label"] = r"$1/\gamma_{\rm c}$ [ps]"
# dark mode coupling dependence
params["gammav"] = [1e-7, 3e-7, 1e-6, 2e-6, 5e-6, 1e-5, 3e-5]
params["gammav_filenames"] = ["../semiclassical/jobs/params_gammav_%.0e.out" % gammav for gammav in params["gammav"]]
params["gammav_label"] = r"$\gamma_{\rm v} \sqrt{N_{\text{dark}}}$ [a.u.]"


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
    xs = [0.02418884254e-3 / x] if task != "gammav" else [x]
    e_dark *= 1e-10
    e_max *= 1e-10
    ys = [e_dark]
    xlog = True
    ylog = True
    xlabel = params[task + "_label"]
    clp.plotone(xs, ys, axes[idx], colors=["ko"], xlabel=xlabel, xlim=None if task != "gammae" else [1e-3, 5e1],
                ylabel=r"$E_{\rm D}$($t = $ 5 ps) [$\times 10^{10}$ cm$^{-1}$]" if idx == 0 else None,
                xlog=xlog, ylog=ylog, showlegend=False)
    if task == "gammae":
        x, y = 0.02418884254e-3 / 1e-5, 7.472
    if task == "gammac":
        x, y = 0.02418884254e-3 / 2e-5, 7.472
    elif task == "gammav":
        x, y = 2e-6, 7.472
    axes[idx].plot(x, y, "ro", mfc='none')


fig, axes = clp.initialize(col=1, row=3, width=5.0, height=5.0*0.618, sharex=False, LaTeX=True,
                           labelthem=True, labelthemPosition=[0.3, 0.94], labelsize=14,
                           return_fig_args=True, fontsize=12, sharey=True)
t_fs = 5000

subplot(axes, idx=0, task="gammae", t_fs=t_fs)
subplot(axes, idx=1, task="gammac", t_fs=t_fs)
subplot(axes, idx=2, task="gammav", t_fs=t_fs)

clp.adjust(tight_layout=True, savefile="params_dephasing.pdf")

