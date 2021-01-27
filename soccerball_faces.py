import numpy as np
import xyzMath

def soccerball_face_centers():
   C0 = 3 * (np.sqrt(5) - 1) / 4
   C1 = 9 * (9 + np.sqrt(5)) / 76
   C2 = 9 * (7 + 5 * np.sqrt(5)) / 76
   C3 = 3 * (1 + np.sqrt(5)) / 4

   return [
      xyzMath.Vec(0.0, C0, C3),
      xyzMath.Vec(0.0, C0, -C3),
      xyzMath.Vec(0.0, -C0, C3),
      xyzMath.Vec(0.0, -C0, -C3),
      xyzMath.Vec(C3, 0.0, C0),
      xyzMath.Vec(C3, 0.0, -C0),
      xyzMath.Vec(-C3, 0.0, C0),
      xyzMath.Vec(-C3, 0.0, -C0),
      xyzMath.Vec(C0, C3, 0.0),
      xyzMath.Vec(C0, -C3, 0.0),
      xyzMath.Vec(-C0, C3, 0.0),
      xyzMath.Vec(-C0, -C3, 0.0),
      xyzMath.Vec(C1, 0.0, C2),
      xyzMath.Vec(C1, 0.0, -C2),
      xyzMath.Vec(-C1, 0.0, C2),
      xyzMath.Vec(-C1, 0.0, -C2),
      xyzMath.Vec(C2, C1, 0.0),
      xyzMath.Vec(C2, -C1, 0.0),
      xyzMath.Vec(-C2, C1, 0.0),
      xyzMath.Vec(-C2, -C1, 0.0),
      xyzMath.Vec(0.0, C2, C1),
      xyzMath.Vec(0.0, C2, -C1),
      xyzMath.Vec(0.0, -C2, C1),
      xyzMath.Vec(0.0, -C2, -C1),
      xyzMath.Vec(1.5, 1.5, 1.5),
      xyzMath.Vec(1.5, 1.5, -1.5),
      xyzMath.Vec(1.5, -1.5, 1.5),
      xyzMath.Vec(1.5, -1.5, -1.5),
      xyzMath.Vec(-1.5, 1.5, 1.5),
      xyzMath.Vec(-1.5, 1.5, -1.5),
      xyzMath.Vec(-1.5, -1.5, 1.5),
      xyzMath.Vec(-1.5, -1.5, -1.5)
   ]

def cgo_arrow(c1, c2, r, col=(0.7, 0.5, 0.5), col2=None):
   c1 = (c1.x, c1.y, c1.z)
   c2 = (c2.x, c2.y, c2.z)
   print(c1)
   print(c2)
   if col2 is None:
      col2 = col
   c = 1.0
   a = 0.5
   b = 0.0
   return [
      cgo.CYLINDER, c * c1[0] + (1 - c) * c2[0], c * c1[1] + (1 - c) * c2[1],
      c * c1[2] + (1 - c) * c2[2], a * c1[0] + (1 - a) * c2[0], a * c1[1] + (1 - a) * c2[1],
      a * c1[2] + (1 - a) * c2[2], r, col[0], col[1], col[2], col2[0], col2[1], col2[2], cgo.CONE,
      a * c1[0] + (1 - a) * c2[0], a * c1[1] + (1 - a) * c2[1], a * c1[2] + (1 - a) * c2[2],
      b * c1[0] + (1 - b) * c2[0], b * c1[1] + (1 - b) * c2[1], b * c1[2] + (1 - b) * c2[2],
      2 * r, 0, col[0], col[1], col[2], col2[0], col2[1], col2[2], 1, 0
   ]

def make_arrow_ab(u, v, r):
   view = cmd.get_view()
   a = com('chain A')
   b = com('chain B')
   cmd.load_cgo(cgo_arrow(u * a + (1 - u) * b, v * a + (1 - v) * b, r), 'arrow')
   cmd.set_view(view)

def draw_soccerball_arrows():
   coords = soccerball_face_centers()
   cgo = []
   for coord in coords:
      cgo += cgo_arrow(coord * 13, coord * 10, 1)
   cmd.load_cgo(cgo, 'soccerball_arrows')
