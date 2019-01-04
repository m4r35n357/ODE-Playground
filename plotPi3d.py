#!/usr/bin/python3

from math import sin, cos, pi
from sys import stdin, stderr, argv
from pi3d import Sphere, Display, Camera, Shader, Keyboard, screenshot, Lines, Font, String


class Body(Sphere):
    def __init__(self, shader, colour, radius, position=(0.0, 0.0, 0.0), track_shader=None, track_max=5000):
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
        self.pointer = Lines(vertices=[(0, 0, 0,), self.pos], material=(1.0, 1.0, 0.0))
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


def radians(degrees):
    return pi * degrees / 180.0


def main():
    # Setup display and initialise pi3d
    display = Display.create(x=0, y=0, frames_per_second=60)
    display.set_background(0, 0, 0, 1)  # r,g,b,alpha
    # Camera
    camera = Camera()
    rot = tilt = 0
    rot_tilt = True
    cam_rad = 50.0
    hud_font = Font('/usr/share/fonts/truetype/liberation2/LiberationMono-Regular.ttf', color='green',
                    codepoints='-0123456789. txyz:=+', font_size=18)
    hud_font.blend = True
    hud_data = " t {:-5.1f}  x {:-5.1f}  y {:-5.1f}  z {:-5.1f} ".format(0.0, 0.0, 0.0, 0.0)
    hud_string = String(camera=Camera(is_3d=False), font=hud_font, is_3d=False, string=hud_data)
    hud_string.set_shader(Shader("uv_flat"))
    (lt, bm, ft, rt, tp, bk) = hud_string.get_bounds()
    hud_string.position((-display.width + rt - lt) / 2.0, (0.9 * display.height - tp + bm) / 2.0, 1.0)
    hud_string.draw()  # NB has to be drawn before quick_change() is called as buffer needs to exist

    if len(argv) > 1:
        particle = Body(Shader("mat_light"), (1.0, 0.0, 0.0), 0.05, track_shader=Shader("mat_flat"),
                        track_max=int(argv[1]))
    else:
        particle = Body(Shader("mat_light"), (1.0, 0.0, 0.0), 0.05, track_shader=Shader("mat_flat"))
    # Enable key presses
    keys = Keyboard()
    # Display scene
    line = stdin.readline()
    while display.loop_running():
        data = line.split(' ')
        hud_string.quick_change(
            " t {:-5.1f}  x {:-5.1f}  y {:-5.1f}  z {:-5.1f} ".format(float(data[3]), float(data[0]), float(data[1]),
                                                                      float(data[2])))
        hud_string.draw()
        # camera control
        if rot_tilt:
            camera.reset()
            camera.rotate(-tilt, rot, 0)
            camera.position((cam_rad * sin(radians(rot)) * cos(radians(tilt)), cam_rad * sin(radians(tilt)),
                             -cam_rad * cos(radians(rot)) * cos(radians(tilt))))
            rot_tilt = False
        # plot the entities
        particle.pos = [float(data[0]), float(data[1]), float(data[2])]
        particle.position_and_draw(trace_material=(0.0, 0.25, 0.25))
        # process keyboard input
        key = keys.read()
        if key > -1:
            rot_tilt = True
            if key == 112:
                screenshot("trajectory.jpg")
            elif key == 119:  # key W rotate camera up
                tilt += 2.0
            elif key == 115:  # kry S down
                tilt -= 2.0
            elif key == 97:  # key A left
                rot -= 2
            elif key == 100:  # key D right
                rot += 2
            elif key == 61:  # key += in
                cam_rad -= 0.5
            elif key == 45:  # key _- out
                cam_rad += 0.5
            elif key == 27:
                keys.close()
                display.stop()
                break
        # prepare for next iteration
        line = stdin.readline()
        if not line:
            display.stop()


if __name__ == "__main__":
    main()
else:
    print(__name__ + " module loaded", file=stderr)
