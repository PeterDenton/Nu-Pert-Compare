import matplotlib.pyplot as plt
import numpy as np

# read in the data
fname = "data/Peak_Precision.txt"
dts = ["name", "dEE1", "dPP1", "dEE2", "dPP2"]
dt = [(d, "f") for d in dts]
dt[0] = (dts[0], "U16")
data = np.loadtxt(fname, dtype = dt)

excludes = ["ZS", "Diag", "Vacuum", "Zeroth"]
mask = [name not in excludes for name in data["name"]]

n = len(data["name"]) - len(excludes)

colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

# plot the data slightly offset horizontally for visual purposes
shift = 0.10
zo = 10
plt.plot([x - shift for x in xrange(n)], data["dPP1"][mask], "*", color = colors[0], zorder = zo)
plt.plot([x - shift for x in xrange(n)], data["dEE1"][mask], "x", color = colors[0], zorder = zo)
plt.plot([x + shift for x in xrange(n)], data["dPP2"][mask], "*", color = colors[1], zorder = zo)
plt.plot([x + shift for x in xrange(n)], data["dEE2"][mask], "x", color = colors[1], zorder = zo)

# logscale in the precision
plt.yscale("log")

# label the yaxis
plt.ylabel(r"$|\Delta x|/x$")

v = list(plt.axis())
v[0] = -0.5
v[1] = n - 0.5
v[2] = 10 ** np.floor(np.log10(v[2]) + 1)
v[3] = 1e-1
plt.axis(v)

# plot outside the figure for the legend
plt.plot([v[1] + 5], [v[3] + 5], "k*", label = r"$P$")
plt.plot([v[1] + 5], [v[3] + 5], "kx", label = r"$E$")
plt.plot([v[1] + 5], [v[3] + 5], ".", color = colors[0], ms = 10, label = r"${\rm First}$")
plt.plot([v[1] + 5], [v[3] + 5], ".", color = colors[1], ms = 10, label = r"${\rm Second}$")

# draw a line at 1%
plt.plot(v[:2], [1e-2, 1e-2], "k:")

# names as xticks
plt.xticks(xrange(n), [r"${\rm %s}$" % name for name in data["name"][mask]], rotation = "vertical")

# lots of yticks
yticks = [10 ** i for i in xrange(int(np.ceil(np.log10(v[2]))), int(np.floor(np.log10(v[3]))) + 1)]
plt.yticks(yticks)

# include a legend
plt.legend(frameon = True, loc = 0, fontsize = 14)

# save the figure
plt.savefig("fig/Peak_Precision.pdf")

