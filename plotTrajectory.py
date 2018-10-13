#!/usr/bin/env python

from json import loads
from math import sqrt, sin, cos, fabs, pi, atan2, acos, log10
from sys import argv, stderr, stdin
from visual import display, sphere, curve, rate, ellipsoid, ring, color, label, points


def error_colour (error):
    if error < -72:
        return (0.35, 0.25, 0.35)
    elif error < -36:
        return (0.5, 0.0, 0.5)
    elif error < -18:
        return color.blue
    elif error < -12:
        return (0.0, 0.5, 0.0)
    elif error < -9:
        return (0.0, 0.5, 0.5)
    elif error < -6:
        return (0.5, 0.5, 0.0)
    elif error < -3:
        return color.orange
    else:
        return color.red


def main():
    print "Trajectory Plotter: {}".format(argv)
    t = int(argv[1])
    x = int(argv[2])
    y = int(argv[3])
    z = int(argv[4]) if len(argv) == 5 else 0
    #  set up the scene
    window_name = 'orbit'
    my_scene = display(title=window_name)
    my_scene.center = centre = (0.0, 0.0, 0.0)
    my_scene.width = my_scene.height = 1024
    my_scene.range = (50.0, 50.0, 50.0)
    # radial line
    radial = curve(pos=[(0.0, 0.0, 0.0), (0.0, 0.0, 0.0)], color=color.gray(0.2))
    hud = label(pos=(0.0, 0.0, 0.0), xoffset=340, yoffset=330, line=False, border=10, font='Monospace', height=16,
                color=(0.5, 0.5, 0.0), linecolor=color.gray(0.1))
    ball = sphere(radius=0.03, color=color.white)  # Origin
    ball = sphere(radius=0.05, color=color.red)  # Particle
    count = 0
    data_line = stdin.readline()
    while data_line:  # build raw data arrays
        rate(60)
        if count % 1000 == 0:
            ball.visible = False
            ball = sphere(radius=0.05, color=color.red)  # Particle
            ball.trail = curve(size=1)  #  trail
        data = data_line.split(' ')
        count += 1
        x2 = float(data[x])**2
        y2 = float(data[y])**2
        z2 = float(data[z])**2
        ball.pos = (float(data[x]), float(data[y]), float(data[z]))
        radial.pos = ((0.0, 0.0, 0.0), ball.pos)
        error = -16
        ball.trail.append(pos=ball.pos, color=(error_colour(error)))
        hud.text = u"t  %.1f\nr  %.1f" % (float(data[t]), sqrt(x2 + y2 + z2))
        #popen('import -window ' + windowName + ' -compress None VPythonOutput/' + str(counter).zfill(4) + '.png')
        data_line = stdin.readline()


if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"
