from matplotlib import pyplot as plt
import numpy as np
import os
import glob

path = os.getcwd()
files = glob.glob(path + "/delta_ms_*")

for file in files:
    name = file.split("/")[-1]
    f = open(file, "r")
    lines = f.readlines()
    f.close()
    x_axis = list(map(lambda x: x[x.find("=")+1:x.find(".")], lines))
    y_axis = list(map(lambda y: y.split()[-1], lines))
    print(x_axis, y_axis)
    plt.plot(x_axis, y_axis)
    plt.xlabel("N")
    plt.ylabel("delta_ms")
    plt.title(name)
    plt.show()
