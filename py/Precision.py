import matplotlib.pyplot as plt
import numpy as np

# get the expressions
fname = "data/Precision.txt"
dataf = open(fname, "r")
expressions = dataf.readline().split()
dataf.close()

# get the data
dts = ["E"]
dts += expressions
dt = [(d, "f") for d in dts]
data = np.loadtxt(fname, dtype = dt, skiprows = 1)

# form several groups for several plots
madrid = ["Madrid", "FL", "AJLOS(31)"]
main = ["Madrid", r"DMP^0", r"DMP^1", r"AM^{5/2}", "AJLOS(48)", "AKT"]
exacts = ["Diag"]
ignores = ["ZS", "Zeroth", "Vacuum"]
others = list(set(expressions) - set(madrid) - set(main) - set(exacts) - set(ignores))

def p(which):
	if which == "madrid":	include = madrid
	if which == "main":		include = main
	if which == "others":	include = others
	if which == "exacts":	include = exacts

	# set zorder and plot
	for i in xrange(1, len(dts)):
		d = dts[i]
		if d not in include: continue
		zo = 100 - i
		if "DMP" in d: zo += 1e5
		plt.plot(data["E"], data[d], label = r"${\rm %s}$" % d, zorder = zo)

	# log-log scale
	plt.xscale("log")
	plt.yscale("log")

	# set the axes
	v = list(plt.axis())
	v[0] = data["E"][0]
	v[1] = data["E"][-1]

	v[2] = 1e-3
	v[3] = 2e0
	if which == "main": v[2] = 1e-7
#	if which == "madrid": v[3] = 1e0
	if which == "exacts":
		v[2] = 1e-16
		v[3] = 1e-12
	plt.axis(v)

	# reasonable looking xticks
	xticks = [0.3, 1., 10.]
	xtick_labels = [r"$%.1f$" % xtick for xtick in xticks]
	for i in xrange(len(xticks)):
		if xticks[i] >= 10.:
			xtick_labels[i] = r"$%i$" % xticks[i]
	plt.xticks(xticks, xtick_labels)

	# lots of yticks
	yticks = [10 ** i for i in xrange(int(np.ceil(np.log10(v[2]))), int(np.floor(np.log10(v[3]))) + 1)]
	plt.yticks(yticks)

	# label axes
	plt.xlabel(r"$E{\rm\ [GeV]}$")
	plt.ylabel(r"$|\Delta P|/P$")

	# include a legend
	legend_loc = (1.02, 1.00)
	plt.legend(ncol = 2, loc = 1, bbox_to_anchor = legend_loc, fontsize = 13, columnspacing = 0.5)

	# save the figure
	if which == "madrid":	plt.savefig("fig/Precision_Madrid.pdf")
	if which == "main":		plt.savefig("fig/Precision.pdf")
	if which == "others":	plt.savefig("fig/Precision_Others.pdf")
	if which == "exacts":	plt.savefig("fig/Precision_Exacts.pdf")
	plt.clf()

p("main")
p("madrid")
p("others")
p("exacts")

