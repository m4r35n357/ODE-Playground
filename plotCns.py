#!/usr/bin/env python3
#
#  (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#
from math import sin, cos, pi, radians
from sys import stdin, stderr, argv
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
    print("CNS Plotter: {}".format(argv))
    if len(argv) < 2:
        raise Exception('>>> ERROR! Please supply two quantities to plot <<<')
    file1 = open(argv[1])
    file2 = open(argv[2])

    # Setup display and initialise pi3d
    display = Display.create(x=0, y=0, frames_per_second=60)
    display.set_background(0, 0, 0, 1)  # r,g,b,alpha
    # Camera
    camera = Camera()
    rot = tilt = 0
    cam_rad = 50.0
    hud_font = Font('/usr/share/fonts/truetype/liberation2/LiberationMono-Regular.ttf', color='green',
                    codepoints='-0123456789. txyz:=+', font_size=18)
    hud_font.blend = True
    hud_data = f' t {0.0:-5.1f}  x {0.0:-5.1f}  y {0.0:-5.1f}  z {0.0:-5.1f} '
    hud_string = String(camera=Camera(is_3d=False), font=hud_font, is_3d=False, string=hud_data)
    hud_string.set_shader(Shader('uv_flat'))
    (lt, bm, ft, rt, tp, bk) = hud_string.get_bounds()
    hud_string.position((-display.width + rt - lt) / 2.0, (0.9 * display.height - tp + bm) / 2.0, 1.0)
    hud_string.draw()  # NB has to be drawn before quick_change() is called as buffer needs to exist

    if len(argv) > 3:
        particle1 = Body(Shader('mat_light'), (0.0, 1.0, 1.0), 0.2, track_shader=Shader('mat_flat'),
                        track_max=int(argv[3]))
        particle2 = Body(Shader('mat_light'), (1.0, 1.0, 0.0), 0.2, track_shader=Shader('mat_flat'),
                         track_max=int(argv[3]))
    else:
        particle1 = Body(Shader('mat_light'), (0.0, 1.0, 1.0), 0.2, track_shader=Shader('mat_flat'))
        particle2 = Body(Shader('mat_light'), (1.0, 1.0, 0.0), 0.2, track_shader=Shader('mat_flat'))
    # Enable key presses and mouse
    keys = Keyboard()
    mouse = Mouse(restrict=False)
    mouse.start()
    omx, omy = mouse.position()
    # Display scene
    line1 = file1.readline()
    line2 = file2.readline()
    while display.loop_running():
        data1 = line1.split(' ')
        data2 = line2.split(' ')
        hud_string.quick_change(
            f' t {float(data1[3]):-5.1f}  x {float(data1[0]):-5.1f}  y {float(data1[1]):-5.1f}  z {float(data1[2]):-5.1f} ')
        hud_string.draw()
        # camera control
        camera.reset()
        camera.rotate(-tilt, rot, 0)
        r_rot, r_tilt = radians(rot), radians(tilt)
        camera.position((cam_rad * sin(r_rot) * cos(r_tilt), cam_rad * sin(r_tilt), - cam_rad * cos(r_rot) * cos(r_tilt)))
        # plot the entities
        particle1.pos = [float(data1[0]), float(data1[1]), float(data1[2])]
        particle1.position_and_draw(trace_material=(0.0, 0.5, 0.0))
        particle2.pos = [float(data2[0]), float(data2[1]), float(data2[2])]
        particle2.position_and_draw(trace_material=(0.5, 0.0, 0.0))
        # process mouse and keyboard input
        mx, my = mouse.position()
        if mouse.button_status() == mouse.LEFT_BUTTON:
            rot += (mx - omx) * 0.2
            tilt += (my - omy) * 0.2
        omx, omy = mx, my
        key = keys.read()
        if key > -1:
            if key == 112:  # 'p'
                screenshot('trajectory.jpg')
            elif key == 119:  # 'w' rotate camera up
                tilt += 2.0
            elif key == 115:  # 's' down
                tilt -= 2.0
            elif key == 97:  # 'a' left
                rot -= 2.0
            elif key == 100:  # 'd' right
                rot += 2.0
            elif key == 61:  # '-' zoom in
                cam_rad -= 0.5
            elif key == 45:  # '=' zoom out
                cam_rad += 0.5
            elif key == 27:  # 'ESC'
                keys.close()
                mouse.stop()
                display.stop()
                break
        # prepare for next iteration
        line1 = file1.readline()
        line2 = file2.readline()
        if not line1:
            display.stop()

if __name__ == "__main__":
    main()
else:
    print(__name__ + " module loaded", file=stderr)
