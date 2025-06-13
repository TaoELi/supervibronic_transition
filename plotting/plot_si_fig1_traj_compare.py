import numpy as np
import columnplots as clp
import matplotlib.pyplot as plt

params = dict()
params["lambdac"] = [0.0, 1e-6, 2e-6, 3e-6, 4e-6]
params["lambdac_labels"] = ["decoupled", r"$\lambda_{\rm c} = 1\times 10^{-6}$ a.u.",
                            r"$\lambda_{\rm c} = 2\times 10^{-6}$ a.u.", r"$\lambda_{\rm c} = 3\times 10^{-6}$ a.u.",
                            r"$\lambda_{\rm c} = 4\times 10^{-6}$ a.u."]
params["lambdac_filenames"] = ["../semiclassical/jobs/params_lambdac_%.0e.out" % lambdac for lambdac in params["lambdac"]]

def import_data(filename):
    data = np.loadtxt(filename)
    t, Pe, xc, xb, xd, eHe, ePh, eVib, eDark, eTot = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6], data[:,7], data[:,8], data[:,9]
    data_dict = {"t": t, "Pe": Pe, "xc": xc, "xb": xb, "xd": xd, "eHe": eHe, "ePh": ePh, "eVib": eVib, "eDark": eDark, "eTot": eTot}
    return data_dict


# prepare data
xs, y1s, y2s = [], [], []
for filename in params["lambdac_filenames"]:
    data = import_data(filename)
    nskip = 1
    t_values = data["t"][:: nskip]
    Pe_lst = data["Pe"][:: nskip]
    xc_lst = data["xc"][:: nskip]
    xb_lst = data["xb"][:: nskip]
    xd_lst = data["xd"][:: nskip]
    e_He_lst = data["eHe"][:: nskip] * 1e-13
    e_ph_lst = data["ePh"][:: nskip] * 1e-13
    e_vib_lst = data["eVib"][:: nskip] * 1e-13
    e_dark_lst = data["eDark"][:: nskip] * 1e-13
    e_tot_lst = data["eTot"][:: nskip] * 1e-13
    xs.append(t_values)
    y1s.append(e_He_lst)
    y2s.append(e_ph_lst)

fig, ax = clp.initialize(col=1, row=1, width=3.4, height=3.4*0.618, return_fig_args=True, LaTeX=True, fontsize=12)

clp.plotone(xs, y1s, ax, ylabel=r"$E_{\rm e}$ [$\times 10^{13}$ cm$^{-1}$]", labels=params["lambdac_labels"],
            xlim=[0, 80], ylim=[0, np.max(y1s[0]) * 1.08], lw=1.0, xlabel=r"time [fs]", showlegend=True)

# change the color of the lines
leg = ax.get_legend()
colormap = plt.cm.hot
colors = [colormap(i) for i in np.linspace(0, 0.6, len(xs))]
for i, j in enumerate(ax.lines):
    j.set_color(colors[i])
    leg.legendHandles[i].set_color(colors[i])
# Put a legend to the right of the current axis
leg = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=True, fancybox=True,
                facecolor='inherit', edgecolor='inherit',
                fontsize=8)

clp.adjust(tight_layout=True, savefile="traj_compare.pdf", includelegend=leg)

