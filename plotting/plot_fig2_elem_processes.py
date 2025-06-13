import numpy as np
import columnplots as clp


def import_data(filename):
    data = np.loadtxt(filename)
    t, Pe, xc, xb, xd, eHe, ePh, eVib, eDark, eTot = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4], data[:,5], data[:,6], data[:,7], data[:,8], data[:,9]
    data_dict = {"t": t, "Pe": Pe, "xc": xc, "xb": xb, "xd": xd, "eHe": eHe, "ePh": ePh, "eVib": eVib, "eDark": eDark, "eTot": eTot}
    return data_dict


def subplot(axes, filename, idx=0, xlim=[0, 5000], add_insert=False):
    data = import_data(filename)
    # get dt from the data
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
    clp.plotone([t_values], [Pe_lst], axes[idx][0], ylabel=r"excited-state population $P_{\rm e}$" if idx == 1 else None,
                xlim=xlim, showlegend=False,
                colors=[clp.navy_blue], lw=1,
                ylim=[0, np.max(Pe_lst)*1.08])
    clp.plotone([t_values] * 4, [e_ph_lst, e_He_lst, e_vib_lst*10, e_dark_lst*100], axes[idx][1],
                ylabel=r"energy [$\times 10^{13}$ cm$^{-1}$]" if idx == 1 else None, showlegend=False,
                colors=[clp.red, clp.navy_blue, clp.dark_green, clp.black], xlim=xlim, lw=1,
                ylim=[0, np.max(e_He_lst)*1.08], yscientific=True)
    if add_insert and idx == 1:
        # add inserted figure in the first plot, zooming in the lines
        from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
        from mpl_toolkits.axes_grid1.inset_locator import mark_inset
        import matplotlib.pyplot as plt
        ax = axes[idx][0]
        axins = zoomed_inset_axes(ax, 180, loc=7)  # zoom = 6
        clp.plotone([t_values], [Pe_lst], axins, showlegend=False, colors=[clp.navy_blue],
                    xlim=[2000,2010], ylim=[0.05, 0.1], lw=0.8)

        # sub region of the original image
        x1, x2, y1, y2 = 100, 105, 0.085, 0.115
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        axins.set_box_aspect(1.0)

        plt.xticks(visible=True)
        plt.yticks(visible=True)
        plt.xticks(fontsize=5)
        plt.yticks(fontsize=5)

        # draw a bbox of the region of the inset axes in the parent axes and
        # connecting lines between the bbox and the inset axes area
        mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

    if add_insert and idx == 1:
        # add inserted figure in the first plot, zooming in the lines
        from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
        from mpl_toolkits.axes_grid1.inset_locator import mark_inset
        import matplotlib.pyplot as plt
        ax = axes[idx][1]
        axins = zoomed_inset_axes(ax, 180, loc=7)  # zoom = 6
        clp.plotone([t_values], [e_ph_lst], axins, showlegend=False, colors=[clp.red],
                    xlim=[100,105], ylim=[0., 0.6], lw=0.8)

        # sub region of the original image
        x1, x2, y1, y2 = 100, 105, 0.0, 0.6
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        axins.set_box_aspect(1.0)

        plt.xticks(visible=True)
        plt.yticks(visible=True)
        plt.xticks(fontsize=5)
        plt.yticks(fontsize=5)

        # draw a bbox of the region of the inset axes in the parent axes and
        # connecting lines between the bbox and the inset axes area
        mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")



# We plot the following four different dynamics
# Fig. a, an isolated electronic subsystem under the pulse excitation: Pe dynamics + energy dynamics (ph, e, vib, dark)
# Fig. b, an isolated coupled VSC subsystem with xc initially excited: xc dynamics + energy dynamics (ph, e, vib, dark)
# Fig. c, a coupled ph+e subsystem under the pulse excitation: Pe dynamics + energy dynamics (ph, e, vib, dark)
# Fig. d, a coupled ph+e+v system under the pulse excitation: Pe dynamics + energy dynamics (ph, e, vib, dark)


fig, axes = clp.initialize(col=3, row=2, width=6, height=6*0.618*1.2, sharex=True, LaTeX=True,
                           labelthem=True, labelthemPosition=[0.95, 0.95], labelsize=14,
                           return_fig_args=True, fontsize=12)
subplot(axes, "../semiclassical/jobs/params_decoupled_e_withfield.out", idx=0, xlim=[0, 5000])
subplot(axes, "../semiclassical/jobs/params_decoupled_eph.out", idx=1, xlim=[0, 5000], add_insert=True)
subplot(axes, "../semiclassical/jobs/params_template.out", idx=2, xlim=[0, 5000], add_insert=True)

axes[2, 0].set_xlabel("time [fs]")
axes[2, 1].set_xlabel("time [fs]")

# add annotation

axes[1, 1].text(0.1, 0.21, "ph", color=clp.red, transform=axes[1, 1].transAxes, fontsize=12)
axes[2, 1].text(0.63, 0.33, r"vib-D$\times 10^{2}$", color=clp.black, transform=axes[2, 1].transAxes, fontsize=12)
axes[2, 1].text(0.03, 0.27, r"vib-B$\times 10$", color=clp.dark_green, transform=axes[2, 1].transAxes, fontsize=12)
axes[2, 1].annotate("", xy=(200, 0.75), xytext=(200, 0.2), arrowprops=dict(arrowstyle="->", color=clp.dark_green))

axes[0, 0].text(0.03, 0.03, 'isolated e-TLSs', horizontalalignment='left', verticalalignment='bottom',
                fontsize=12, color="0.65",
                transform=axes[0, 0].transAxes)
axes[1, 0].text(0.03, 0.03, 'e-TLSs + ph', horizontalalignment='left', verticalalignment='bottom',
                fontsize=12, color="0.65",
                transform=axes[1, 0].transAxes)
axes[2, 0].text(0.03, 0.03, 'e-TLSs + ph + vib', horizontalalignment='left', verticalalignment='bottom',
                fontsize=12, color="0.65",
                transform=axes[2, 0].transAxes)

clp.adjust(tight_layout=True, savefile="elem_processes.pdf")

