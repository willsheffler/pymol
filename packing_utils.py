import pymol
from pymol import cmd
import sys
import os
import random
from math import floor

rainbow = ['0xCCFF00', '0xFF0000', '0xCC00FF', '0x00FF66', '0x80FF00', '0x0019FF', '0x00FF19', '0x33FF00', '0xFF0099', '0xFF004D',
           '0xFF00E6', '0x00FFFF', '0x0066FF', '0x8000FF', '0x00B3FF', '0xFFE500', '0x00FFB2', '0xFF4C00', '0x3300FF', '0xFF9900']


def useOccColors(sel="all"):
    d = {}
    for a in cmd.get_model().atom:
        d[a.q] = True
    colors = rainbow
   # random.shuffle(colors)
   # colors *= len(d)/len(colors)+1
    for ii in range(len(d.keys())):
        cmd.color(colors[ii], "%s and q=%i" % (sel, d.keys()[ii]))


def useTempColors(sel="all"):
    for a in cmd.get_model(sel).atom:
        q = a.b
        c = intcolors[int(floor(q)) % len(intcolors)]
        cmd.color(c, "%s and resi %s and name %s" % (sel, a.resi, a.name))


def useOccRadii(sel="all"):
    for a in cmd.get_model(sel).atom:
        q = a.q
        if q >= 3:
            print "shrik radius"
            q < - 0.1
        cmd.alter("%s and resi %s and name %s" %
                  (sel, a.resi, a.name), "vdw=%f" % (q))
    cmd.rebuild()


def useTempRadii(sel="all"):
    for ii in range(30):
        radius = "%0.1f" % (float(ii + 1) / 10)
        cmd.alter(sel + " and b=" + radius, "vdw=" + radius)
    cmd.rebuild()


def loadPackingPDB(file, name=None, native=None):
    """
    usage: loadPackingPDB <file> , [<name for object>]
    loads a foo_packing.pdb file and colors it all pretty-like
    creates two selections along with the loaded object called
    NAMEcavities and NAMEprotein which are the heteratoms representing
    holes and everything else, respectively. Names can get pretty long,
    by pymol lets you do good stuff like "select NA*cav*", which will
    match a selection called NAMEISREALLYLONGcavities.
    """

    if name is None:
        name = name = os.path.basename(file)
    if name.endswith('.gz'):
        name = name[:-3]
    if name.endswith('.pdb'):
        name = name[:-4]
    if name.endswith('.'):
        name = name[:-1]
    if name.endswith("_packing"):
        name = name[:-8]

    zload(file, name)
    cmd.hide('everything', name)

    if native is not None:
        cmd.align(name, native)
        cmd.zoom(native)

    useRosettaRadii()

    cavselname = name + "cavities"
    protselname = name + "protein"

    cmd.select(cavselname, "resn WSS and %s" % (name))
    cmd.select(protselname, "(not resn WSS) and %s" % (name))

    useTempRadii(cavselname)
    useOccColors(cavselname)
    cmd.color("white", protselname)

    cmd.show('spheres', cavselname)
    cmd.show("cartoon", protselname)
    cmd.show("lines", protselname)

    cmd.select("none")
    cmd.delete("sele*")
    cmd.move('z', -50)

    return name


cmd.extend("loadPackingPDB", loadPackingPDB)
cmd.extend("useOccRadii", useOccRadii)
cmd.extend("useOccColors", useOccColors)
cmd.extend("useTempRadii", useTempRadii)
cmd.extend("useTempColors", useTempColors)
