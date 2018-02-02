from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

def openF(s):
    p = []
    with open(s) as input_file:
        for line in input_file:
            line = line.strip()
            for number in line.split():
                p.append(float(number))
    return p
p = openF("test.txt")
a = [x for i, x in enumerate(p) if i % 3 == 0]
b = [x for i, x in enumerate(p) if i % 3 == 1]
c = [x for i, x in enumerate(p) if i % 3 == 2]
#a, b = np.meshgrid(a, b)
fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.scatter3D(a, b, c, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
plt.show()


