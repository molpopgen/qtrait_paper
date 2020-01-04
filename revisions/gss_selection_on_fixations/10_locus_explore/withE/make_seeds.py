import numpy as np
import sys

nseeds = int(sys.argv[1])
seeds = []

for i in range(nseeds):
    x = np.random.randint(0, int(4e7), 1)[0]
    while x in seeds:
        x = np.random.randint(0, int(4e7), 1)[0]
    seeds.append(x)

for i in seeds:
    print(i)
