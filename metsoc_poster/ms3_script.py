import duckdb, matplotlib.pyplot as plt, pandas as pd

con = duckdb.connect("ms3_data/ms3_data.duckdb")
ms1 = con.execute("SELECT * FROM MS1 WHERE mz BETWEEN 282.1186 AND 282.1208").fetchdf()
ms2 = con.execute("SELECT * FROM MS2 WHERE premz BETWEEN 282.1186 AND 282.1208 ORDER BY int DESC").fetchdf()
ms3 = con.execute("""SELECT * FROM MS3 WHERE prepremz BETWEEN 282.1186 AND 282.1208 AND 
                  premz BETWEEN 166.0718 AND 166.0731""").fetchdf()
con.close()

# fig, axs = plt.subplots(3, 1, figsize=(9, 6))
# axs[0].plot(ms1['rt'], ms1['int'])
# axs[0].vlines(ms2['rt'], ymin=0, ymax=ms1['int'].max())
# axs[1].stem(ms2['fragmz'], ms2['int'])
# axs[2].stem(ms3['fragmz'], ms3['int'])
# plt.show()

ms1["int"] /= 1e6
ms2["int"] /= 1e6
ms3["int"] /= 1e6

fig, axs = plt.subplots(3, 1, figsize=(9, 6), sharex=False)
fig.subplots_adjust(hspace=0.4)

axs[0].plot(ms1['rt'], ms1['int'], color="black", lw=1)
axs[0].vlines(ms2['rt'], ymin=0, ymax=ms1['int'].max(), color="#a41118", lw=1.5)
axs[0].vlines(sorted(ms2['rt'].unique())[2], ymin=0, ymax=ms1['int'].max(), color="#e0670b", lw=1.5)
ms1_lab = """MS$^1$ Extracted
Ion Chromatogram"""
axs[0].text(0.01, 0.95, ms1_lab, va="top", ha="left", transform=axs[0].transAxes)
axs[0].set_xlabel("Retention time")
axs[0].set_ylabel(r"Intensity ($10^6$)")

axs[1].stem(ms2['fragmz'], ms2['int'], linefmt="#e0670b", basefmt='k')
axs[1].stem(ms2['fragmz'][0], ms2['int'][0], linefmt="#f4bb23", basefmt='k')
axs[1].text(0.01, 0.95, r"MS$^2$ Spectrum", va="top", ha="left", transform=axs[1].transAxes)
axs[1].set_xlabel("")
axs[1].set_ylabel(r"Intensity ($10^6$)")

axs[2].stem(ms3['fragmz'], ms3['int'], linefmt='#f4bb23', basefmt='k')
axs[2].text(0.01, 0.95, r"MS$^3$ Spectrum", va="top", ha="left", transform=axs[2].transAxes)
axs[2].set_xlabel(r"Fragment m/z")
axs[2].set_ylabel(r"Intensity ($10^6$)")

axs[2].sharex(axs[1])

for ax in axs:
    ax.ticklabel_format(style='plain', axis='x')
    ax.ticklabel_format(style='plain', axis='y')
    xticks = ax.get_xticks()
    if len(xticks) < 2:
        xmin, xmax = ax.get_xlim()
        ax.set_xticks([xmin, xmax])
    yticks = ax.get_yticks()
    if len(yticks) < 2:
        ymin, ymax = ax.get_ylim()
        ax.set_yticks([ymin, ymax])


plt.savefig("ms3_duckdb_python_fig.png", dpi=600, bbox_inches="tight")
