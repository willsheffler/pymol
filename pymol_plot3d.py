from pymol import cmd
from pymol.cgo import *

def plot3d(ary):
   assert ary.shape[1] == 3
   cmd.delete("plot")
   print("plotting", ary.shape[0], "points")
   obj = list()
   obj.extend([BEGIN, POINTS, COLOR, 1.0, 1.0, 1.0])
   for i in range(ary.shape[0]):
      obj.extend([VERTEX, ary[i, 0], ary[i, 1], ary[i, 2]])
   cmd.load_cgo(obj, "plot")
