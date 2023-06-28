#!/usr/bin/env python3
#
#  (c) 2018-2023 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

from sys import stderr, argv, stdin
from math import sin, cos, radians
from pi3d import Sphere, Display, Camera, Shader, Keyboard, screenshot, Lines, Font, String, Mouse

class Body(Sphere):
    def __init__(self, shader, colour, radius, position=(0.0, 0.0, 0.0), track_shader=None, track_max=1000):
        super(Body, self).__init__(radius=radius, x=position[0], y=position[1], z=position[2])
        super(Body, self).set_draw_details(shader, [])
        self.pos = position
        self.track_shader = track_shader
        self.vertices = []
        self.track_length = 0
        self.track_max = track_max
        self.track = None
        self.set_material(colour)
        self.pointer = Lines(vertices=[(0, 0, 0,), (0, 0, 0,)])

    def position_and_draw(self, trace_material):
        # body
        self.position(self.pos[0], self.pos[1], self.pos[2])
        self.draw()
        self.pointer = Lines(vertices=[(0, 0, 0,), self.pos])
        self.pointer.draw()
        # track
        if self.track_shader:
            self.track_length += 1
            self.vertices.append(self.pos)
            self.track = Lines(vertices=self.vertices, material=trace_material)
            self.track.set_shader(self.track_shader)
            if self.track_length > self.track_max:
                self.vertices = self.vertices[-self.track_max:]
                self.track_length = self.track_max
            if self.track:
                self.track.draw()

def main():
    print(f'Multi 3D ODE Plotter: {argv}', file=stderr)
    argc = len(argv) - 1
    if argc == 0 or argc == 1:  # single particle plot from stdin, optional arg is track length
        files = [stdin]
    elif argc == 7:  # called by IC script, 1 + 6 particle file names
        files = []
        for arg in argv[1:]:
            files.append(open(arg))
    else:
        raise Exception('>>> ERROR! Please supply 1, 2 or 7 arguments! <<<')
    # Setup display and initialise pi3d
    display = Display.create(x=0, y=0, frames_per_second=60)
    display.set_background(0, 0, 0, 1)  # r,g,b,alpha
    # Camera
    camera = Camera()
    rot, tilt = 135.0, 90.0 - 54.73561
    cam_rad = 50.0
    font = Font('/usr/share/fonts/truetype/liberation2/LiberationMono-Regular.ttf', color='green', codepoints='-0123456789. txyz:=+', font_size=18)
    font.blend = True
    hud = String(camera=Camera(is_3d=False), font=font, is_3d=False, string=f' t {0.0:-5.1f}  x {0.0:-5.1f}  y {0.0:-5.1f}  z {0.0:-5.1f}')
    hud.set_shader(Shader('uv_flat'))
    (lt, bm, ft, rt, tp, bk) = hud.get_bounds()
    hud.position((-display.width + rt - lt) / 2.0, (0.9 * display.height - tp + bm) / 2.0, 1.0)
    hud.draw()  # NB has to be drawn before quick_change() is called as buffer needs to exist
    particles = []
    if argc == 0:
        particles.append(Body(Shader('mat_light'), (0.0, 1.0, 1.0), 0.1, track_shader=Shader('mat_flat')))
    elif argc == 1:
        particles.append(Body(Shader('mat_light'), (0.0, 1.0, 1.0), 0.05, track_shader=Shader('mat_flat'), track_max=int(argv[1])))
    elif argc == 7:
        particles.append(Body(Shader('mat_light'), (0.0, 1.0, 1.0), 0.1, track_shader=Shader('mat_flat')))
        particles.append(Body(Shader('mat_light'), (0.0, 1.0, 1.0), 0.1, track_shader=Shader('mat_flat')))
        particles.append(Body(Shader('mat_light'), (1.0, 1.0, 0.0), 0.1, track_shader=Shader('mat_flat')))
        particles.append(Body(Shader('mat_light'), (1.0, 1.0, 0.0), 0.1, track_shader=Shader('mat_flat')))
        particles.append(Body(Shader('mat_light'), (1.0, 0.0, 1.0), 0.1, track_shader=Shader('mat_flat')))
        particles.append(Body(Shader('mat_light'), (1.0, 0.0, 1.0), 0.1, track_shader=Shader('mat_flat')))
        particles.append(Body(Shader('mat_light'), (1.0, 1.0, 1.0), 0.1, track_shader=Shader('mat_flat')))
    # Enable key presses and mouse
    keys = Keyboard()
    mouse = Mouse(restrict=False)
    mouse.start()
    omx, omy = mouse.position()
    # Display scene
    while display.loop_running():
        # prepare for next iteration
        lines = []
        for file in files:
            lines.append(file.readline())
        if not lines[0]:
            display.stop()
        data = [[float(item) for item in lines[0].split()[:4]]]
        if argc == 7:
            data = []
            for line in lines:
                data.append([float(item) for item in line.split()[:4]])
        hud.quick_change(f' t{data[0][3]:-5.1f}  x{data[0][0]:-5.1f}  y{data[0][1]:-5.1f}  z{data[0][2]:-5.1f}')
        hud.draw()
        # camera control
        camera.reset()
        camera.rotate(-tilt, rot, 0)
        r_rot, r_tilt = radians(rot), radians(tilt)
        camera.position((cam_rad * sin(r_rot) * cos(r_tilt), cam_rad * sin(r_tilt), - cam_rad * cos(r_rot) * cos(r_tilt)))
        # plot the particles
        particles[0].pos = [data[0][0], data[0][1], data[0][2]]
        particles[0].position_and_draw(trace_material=(0.0, 0.25, 0.0))
        if argc == 7:
            particles[1].pos = [data[1][0], data[1][1], data[1][2]]
            particles[1].position_and_draw(trace_material=(0.0, 0.25, 0.0))
            particles[2].pos = [data[2][0], data[2][1], data[2][2]]
            particles[2].position_and_draw(trace_material=(0.4, 0.0, 0.0))
            particles[3].pos = [data[3][0], data[3][1], data[3][2]]
            particles[3].position_and_draw(trace_material=(0.4, 0.0, 0.0))
            particles[4].pos = [data[4][0], data[4][1], data[4][2]]
            particles[4].position_and_draw(trace_material=(0.0, 0.0, 0.5))
            particles[5].pos = [data[5][0], data[5][1], data[5][2]]
            particles[5].position_and_draw(trace_material=(0.0, 0.0, 0.5))
            particles[6].pos = [data[6][0], data[6][1], data[6][2]]
            particles[6].position_and_draw(trace_material=(0.25, 0.25, 0.25))
        # process mouse & keyboard input
        mx, my = mouse.position()
        if mouse.button_status() == mouse.LEFT_BUTTON:
            rot -= (mx - omx) * 0.2
            tilt -= (my - omy) * 0.2
        elif mouse.button_status() == mouse.RIGHT_BUTTON:
            cam_rad += (my - omy) * 0.1
        omx, omy = mx, my
        key = keys.read()
        if key > -1:
            if key == 112:  # 'p'
                screenshot('trajectory.jpg')
            elif key == 27:  # 'ESC'
                keys.close()
                mouse.stop()
                display.stop()
                break

main()
