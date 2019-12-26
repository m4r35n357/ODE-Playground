#!/usr/bin/env python3

#
#  (c) 2018-2020 m4r35n357@gmail.com (Ian Smith), for licencing see the LICENCE file
#

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

    def position_and_draw(self, trace_material, origin=(0, 0, 0,)):
        # body
        self.position(self.pos[0], self.pos[1], self.pos[2])
        self.draw()
        self.pointer = Lines(vertices=[origin, self.pos], material=(1.0, 1.0, 0.0))
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
    # Setup display and initialise pi3d
    display = Display.create(x=0, y=0, frames_per_second=60)
    display.set_background(0, 0, 0, 1)  # r,g,b,alpha
    # Camera
    camera = Camera()
    cam_rad = 50.0
    hud_font = Font('/usr/share/fonts/truetype/liberation2/LiberationMono-Regular.ttf', color='green',
                    codepoints='-0123456789. ehtxyz:=+', font_size=18)
    hud_font.blend = True
    hud_headings = " t {:-5.1f} "
    hud_data = hud_headings.format(0.0)
    hud_string = String(camera=Camera(is_3d=False), font=hud_font, is_3d=False, string=hud_data)
    hud_string.set_shader(Shader("uv_flat"))
    (lt, bm, ft, rt, tp, bk) = hud_string.get_bounds()
    hud_string.position((-display.width + rt - lt) / 2.0, (0.9 * display.height - tp + bm) / 2.0, 1.0)
    hud_string.draw()  # NB has to be drawn before quick_change() is called as buffer needs to exist

    if len(argv) > 1:
        particle1 = Body(Shader("mat_light"), (1.0, 0.0, 0.0), 0.05, track_max=int(argv[1]))
        particle2 = Body(Shader("mat_light"), (1.0, 0.0, 0.0), 0.05, track_shader=Shader("mat_flat"), track_max=int(argv[1]))
    else:
        particle1 = Body(Shader("mat_light"), (1.0, 0.0, 0.0), 0.05)
        particle2 = Body(Shader("mat_light"), (1.0, 0.0, 0.0), 0.05, track_shader=Shader("mat_flat"))
    # Enable key presses
    keys = Keyboard()
    # Display scene
    camera.reset()
    camera.position((0.0, 0.0, -cam_rad))
    camera_moved = False
    line = stdin.readline()
    while display.loop_running():
        data = line.split(' ')
        hud_string.quick_change(hud_headings.format(float(data[4])))
        hud_string.draw()
        # camera control
        if camera_moved:
            camera.reset()
            camera.position((0.0, 0.0, -cam_rad))
            camera_moved = False
        # plot the entities
        particle1.pos = [float(data[0]), float(data[1]), 0.0]
        particle1.position_and_draw(trace_material=(0.0, 0.25, 0.25))
        particle2.pos = [float(data[2]), float(data[3]), 0.0]
        particle2.position_and_draw(trace_material=(0.0, 0.25, 0.25), origin=particle1.pos)
        # process keyboard input
        key = keys.read()
        if key > -1:
            if key == 112:  # 'p'
                screenshot("trajectory.jpg")
            elif key == 61:  # '-' zoom in
                cam_rad -= 0.5
                camera_moved = True
            elif key == 45:  # '=' zoom out
                cam_rad += 0.5
                camera_moved = True
            elif key == 27:  # 'ESC'
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
