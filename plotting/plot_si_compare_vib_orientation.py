import numpy as np
import columnplots as clp
import matplotlib.pyplot as plt

def import_data(filename):
    data = np.loadtxt(filename)
    t, Pe, xc, xb, xd, eHe, ePh, eVib, eDark, eTot = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6], data[:,7], data[:,8], data[:,9]
    data_dict = {"t": t, "Pe": Pe, "xc": xc, "xb": xb, "xd": xd, "eHe": eHe, "ePh": ePh, "eVib": eVib, "eDark": eDark, "eTot": eTot}
    return data_dict


def subplot(axes, filename, idx=0, xlim=[0, 5000], Ne=500):
    data = import_data(filename)
    # get dt from the data
    nskip = 1
    t_values = data["t"][:: nskip]
    Pe_lst = data["Pe"][:: nskip]
    xc_lst = data["xc"][:: nskip]
    xb_lst = data["xb"][:: nskip]
    xd_lst = data["xd"][:: nskip]
    e_He_lst = data["eHe"][:: nskip] / Ne
    e_ph_lst = data["ePh"][:: nskip] / Ne
    e_vib_lst = data["eVib"][:: nskip] / Ne
    e_dark_lst = data["eDark"][:: nskip] / Ne
    e_vib_tot_lst = e_vib_lst + e_dark_lst
    e_tot_lst = data["eTot"][:: nskip] / Ne
    clp.plotone([t_values] * 2, [e_ph_lst, e_He_lst], axes[idx][0], ylabel=r"energy / $N_{\rm e}$ [cm$^{-1}$]",
                xlim=xlim, showlegend=True if idx==0 else False, labels=["ph", "e"], legendloc="center right",
                colors=[clp.red, clp.navy_blue], lw=1,
                ylim=[0, np.max(e_He_lst)*1.08])
    clp.plotone([t_values], [e_vib_tot_lst], axes[idx][1], legendloc="center right",
                ylabel=r"energy / $N_{\rm e}$ [cm$^{-1}$]", showlegend=True if idx==0 else False, labels=["total vib"],
                colors=[clp.brown], xlim=xlim, lw=1,
                ylim=[0, np.max(e_vib_tot_lst)*1.08], yscientific=True)



fig, axes = clp.initialize(col=3, row=2, width=6, height=6*0.618*1.5, sharex=True, LaTeX=True,
                           labelthem=True, labelthemPosition=[0.95, 0.95], labelsize=14,
                           return_fig_args=True, fontsize=12)

subplot(axes, "../semiclassical/jobs_vib_orientations/params_angle_0.000.out", idx=0, xlim=[0, 5000])
subplot(axes, "../semiclassical/jobs_vib_orientations/params_angle_0.000_vangle.out", idx=1, xlim=[0, 5000])
subplot(axes, "../semiclassical/jobs_vib_orientations/params_angle_3.142_vangle.out", idx=2, xlim=[0, 5000])

axes[2, 0].set_xlabel("time [fs]")
axes[2, 1].set_xlabel("time [fs]")

# add annotation


axes[0, 0].text(0.13, 0.7, "uniform vib dipoles\nvib-B/D modes", horizontalalignment='left', verticalalignment='bottom',
                fontsize=12, color="0.65",
                transform=axes[0, 0].transAxes)
axes[1, 0].text(0.13, 0.7, "uniform vib dipoles\nvib-oscillators", horizontalalignment='left', verticalalignment='bottom',
                fontsize=12, color="0.65",
                transform=axes[1, 0].transAxes)
axes[2, 0].text(0.13, 0.7, "random vib dipoles\nvib-oscillators", horizontalalignment='left', verticalalignment='bottom',
                fontsize=12, color="0.65",
                transform=axes[2, 0].transAxes)

clp.adjust(tight_layout=True, savefile="compare_vib_orientation.pdf")

