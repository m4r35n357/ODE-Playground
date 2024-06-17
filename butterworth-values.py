from sys import argv
from math import pi, sin

order = int(argv[1])
for i in range(order):
    print(f'{2.0 * sin((2.0 * (i + 1.0) - 1.0) * pi / (2 * order)): .6f} ', end='')
print()
