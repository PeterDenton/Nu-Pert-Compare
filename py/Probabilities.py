import matplotlib.pyplot as plt
import numpy as np

# get the expressions
fname = "data/Probabilities.txt"
dataf = open(fname, "r")
expressions = dataf.readline().split()
dataf.close()

# get the data
dts = ["E"]
dts += expressions
dts += [expression + "bar"  for expression in expressions]
dt = [(d, "f") for d in dts]
data = np.loadtxt(fname, dtype = dt, skiprows = 1)

# don't plot every function
exclude = [r"AM^2", "AKS", "Diag", "AJLOS(31)", "MP", "MF", "FL", "Zeroth", "Vacuum", r"DMP^1", "AKT"]

# two plots
f, (ax1, ax2) = plt.subplots(1, 2, sharey = True)
f.subplots_adjust(wspace = 0)

for i in xrange(1, len(dts)):
	d = dts[i]

	# make the exact expression identifiable
	if "ZS" in d:
		color = "k"
		zorder = 10
		ls = "--"
	else:
		color = None
		zorder = 0
		ls = "-"

	if d in exclude: continue
	if "bar" in d: continue

	l, = ax1.plot(data["E"], data[d], color = color, zorder = zorder, ls = ls)
	ax2.plot(data["E"], data[d + "bar"], color = l.get_color(), label = r"${\rm %s}$" % d, zorder = zorder, ls = ls)

# log energy scale
ax1.set_xscale("log")
ax2.set_xscale("log")

# make the axes range
v = list(ax1.axis())
v[0] = 0.5
v[1] = round(data["E"][-1], 1)
v[2] = max(v[2], 0)
v[3] = 0.10

# reasonable looking xticks
xticks = [1., 10.]
xtick_labels = [r"$%.1f$" % xtick for xtick in xticks]
for i in xrange(len(xticks)):
	if xticks[i] >= 10.:
		xtick_labels[i] = r"$%i$" % xticks[i]
ax1.set_xticks(xticks)
ax1.set_xticklabels(xtick_labels)
ax2.set_xticks(xticks)
ax2.set_xticklabels(xtick_labels)

ax1.get_xaxis().majorTicks[1].label1.set_horizontalalignment("right")

# label axes
ax1.set_xlabel(r"$E{\rm\ [GeV]}$")
ax2.set_xlabel(r"$E{\rm\ [GeV]}$")
ax1.set_ylabel(r"$P_{\mu e}$")

# set the axes dimensions
ax1.axis(v)
ax2.axis(v)

# describe which figure is neutrinos and which is antineutrinos
ax1.text(3, 0.005, r"$\nu$")
ax2.text(3, 0.005, r"$\bar\nu$")

# add a legend to the nubar figure
ax2.legend()

# save the figure
plt.savefig("fig/Probabilities.pdf")

