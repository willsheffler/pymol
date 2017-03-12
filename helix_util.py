# -*- mode:python;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
import sys
import os
import inspect
import functools
newpath = os.path.dirname(inspect.getfile(
    inspect.currentframe()))  # script directory
if not newpath in sys.path:
    sys.path.append(newpath)
import string
import re
import gzip
import itertools
import collections
from sym_util import *
import operator as op
from xyzMath import *
from itertools import product, ifilter
import cProfile


def xform_covers_all_coms(axis, cen, coms, xform, nfold, show=False):
    # TODO: add NFOLD!!!!!!!!!!!!!!
    syms = [RAD(axis, i * 360.0 / nfold, cen) for i in range(nfold)]
    seenit = [False] * len(coms)
    error2 = 0.0
    for i, in_com in enumerate(coms):
        hcom = coms[0]
        for j in range(222):
            for sym in syms:
                testcom = sym * hcom
                if show:
                    showsphere(testcom)
                if in_com.distance_squared(testcom) < 0.1:
                    # print in_com.distance(hcom)
                    seenit[i] = True
                    error2 += in_com.distance_squared(testcom)
                    break
            hcom = xform * hcom
    if all(seenit):
        return sqrt(error2 / len(coms))
    else:
        return None


def get_correction_angle(axis, cen, coms, unit_xform, error):
    correction_ang = None
    nturns = None
    correction_ang_axis_dis = 0.0
    for i in range(len(coms) - 1, 1, -1):
        in_com = coms[i]
        print "try to get correction_ang", i, in_com, len(coms)
        hcom = coms[0]
        axis_dis = in_com.dot(axis)
        if axis_dis < correction_ang_axis_dis:
            continue
        for j in range(222):
            # if i==len(coms)-1: showsphere( hcom )
            if in_com.distance(hcom) < 2.0:
                # print coms[-1].distance(hcom)
                v1_0 = in_com - cen
                v2_0 = hcom - cen
                v1 = (v1_0 - v1_0.dot(axis) * axis).normalized()
                v2 = (v2_0 - v2_0.dot(axis) * axis).normalized()
                # print v1.dot(axis)
                # print v2.dot(axis)
                correction_ang = acos(v1.dot(v2))
                correction_ang_axis_dis = axis_dis
                # TODO: check this logic!
                if v1.cross(v2).dot(axis) < 0:
                    correction_ang *= -1.0
                print "CORRECTION ANGLE:", i, correction_ang, axis_dis
                nturns = j
                break
            hcom = unit_xform * hcom
    return correction_ang, nturns, error


def determine_helix_geometry(sele='vis', show=0, verbose=True, nfolds=None, n_to_contacts=None):
    """
    # delete all; load /Users/sheffler/Downloads/N4_C1_DR53_001.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; determine_helix_geometry(show=40)
    # delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N2_C1_DR04_065.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; determine_helix_geometry(show=50)
    # delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N2_C3_DR04_003.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; determine_helix_geometry(show=50)
    # delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N3_C1_DR04_002.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; determine_helix_geometry(show=50)
    # delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N3_C2_DR04_001.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; determine_helix_geometry(show=50)
    # delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N3_C4_DR04_001.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; determine_helix_geometry(show=50)
    # delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N4_C3_DR04_001.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; determine_helix_geometry(show=50)
    # delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N5_C1_DR04_002.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; determine_helix_geometry(show=50)
    # delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N4_C1_DR04_002.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; determine_helix_geometry(show=50)

    # delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test2/N1_C2_DR18_001_0023.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; determine_helix_geometry(show=50)
    # delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test2/N1_C5_DR04_035_0769.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; determine_helix_geometry(show=50)
    # delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test2/N3_C1_DR54_151_0104.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; determine_helix_geometry(show=50)
    # delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test2/N5_C1_DR10_029_0419.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; determine_helix_geometry(show=50)

    # delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_chainAB/N1_C2_DR79_078_0275.pose.pdb.gz; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; determine_helix_geometry(show=50)

    # STILL A PROBLEM:
    # N2_C1_DR14_001
    # N2_C3_DR14_001
    # N2_C3_DR14_008
    # N2_C3_DR64_004
    # N3_C1_DR05_008
    # N3_C1_DR10_006
    # N3_C1_DR14_003
    # N3_C1_DR14_004
    # N3_C2_DR14_004
    # N3_C4_DR10_003
    # N3_C4_DR10_006
    # N3_C4_DR14_002
    # N3_C4_DR14_003
    # N3_C4_DR14_004
    # N3_C4_DR14_007
    # N3_C4_DR49_002
    # N3_C4_DR49_008
    # N3_C4_DR54_002
    # N3_C4_DR54_004
    # N3_C4_DR54_005
    # N3_C4_DR54_006
    # N3_C4_DR54_007
    # N3_C4_DR76_005
    # N3_C4_DR76_006
    # N4_C1_DR10_005
    # N4_C1_DR14_004
    # N4_C1_DR14_008
    # N5_C1_DR14_001 -- div by zero
    # N5_C1_DR14_003 --
    # 35
    # /Users/sheffler/tmp/waiting_list/N2_C3_DR08_001.pdb
    # delete all; load /Users/sheffler/tmp/waiting_list/N2_C1_DR05_006.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; determine_helix_geometry(show=50)

    """

    if nfolds is None:
        nfolds = range(1, 7)
    elif isinstance(nfolds, int):
        nfolds = (nfolds,)

    if n_to_contacts is None:
        n_to_contacts = range(1, 4)
    elif isinstance(n_to_contacts, int):
        n_to_contacts = range(1, n_to_contacts)

    error = None

    chains = cmd.get_chains(sele)
    print "chains:", chains
    # make sure all chains have same len
    nres = [len(getres("((%s) and name ca and chain %s)" % (sele, c)))
            for c in chains]
    print "nres foreach chain:", nres
    if min(nres) != max(nres):
        raise Exception("all chains must be same length")

    # get centers-of-mass of all chains
    coms = [com("((%s) and name ca and chain %s)" % (sele, c)) for c in chains]

    # make sure sorted along axis
    if True:
        axis_tmp, ang1_tmp, cen_tmp = getrelframe('((%s) and name ca and chain %s)' % (
            sele, chains[1]), '((%s) and name ca and chain %s)' % (sele, chains[0])).rotation_axis_center()
        com_all = com("((%s) and name ca)" % sele)
        if axis_tmp.dot(com_all - cen_tmp) < 0:
            axis_tmp = Vec(0, 0, 0) - axis_tmp
        coms_tosort = [(axis_tmp.dot(coms[i] - cen_tmp), coms[i], chains[i])
                       for i in range(len(coms))]
        coms = [xyz for t, xyz, c in sorted(coms_tosort)]
        chains = [c for t, xyz, c in sorted(coms_tosort)]

    axis, ang1, cen = None, None, None
    if True:  # this is just to limit scope
        # get initial gusss at helix axis, rotation angle, center
        axis, ang1, cen = getrelframe('((%s) and name ca and chain %s)' % (
            sele, chains[1]), '((%s) and name ca and chain %s)' % (sele, chains[0])).rotation_axis_center()
        # make sure axis points toward bulk of helix
        com_all = com("((%s) and name ca)" % sele)
        if axis.dot(com_all - cen) < 0:
            axis = Vec(0, 0, 0) - axis
            ang1 = 2.0 * math.pi - ang1

        # this causes some problem I can't figure out...
        # # try to get better estimate of axis and cen by checking all chains vs. A
        # # probably not necessary?
        # print "AXIS0:", axis
        # print "CEN0:", cen
        # print "ANG0:", ang1
        # axis0 = axis
        # cen0 = cen
        # ang1_0 = ang1
        # ncen = 1
        # for c in chains[2:]:
        #   axisB, ang1B, cenB = getrelframe( '((%s) and name ca and chain %s)'%(sele,c), '((%s) and name ca and chain %s)'%(sele,chains[0]) ).rotation_axis_center()
        #   if not cenB: continue
        #   if   axis.dot(axisB) >  0.9: axis += axisB
        #   elif axis.dot(axisB) < -0.9: axis -= axisB
        #   if cen.distance( cenB ) < 0.1:
        #       cen += cenB
        #       ncen += 1
        # axis.normalize()
        # cen = cen / float( ncen )
        # ang1 = dihedral( coms[0], cen, cen+axis, coms[1] )
        # if abs( ang1 - ang1_0 ) > 0.01:   ang1 = 2.0*math.pi + ang1
        # if axis.dot( axis0 ) < 0.999: raise Exception("axis and axis0 are too far apart for sanity")
        # if cen.distance( cen0 ) > 0.01: raise Exception("axis and axis0 are too far apart for sanity")
        # if abs( ang1 - ang1_0 ) > 0.01: raise Exception("ang1 and ang1_0 are too far apart for sanity")

    if verbose:
        print "AXIS:", axis
    if verbose:
        print "CEN:", cen
    if verbose:
        print "ANG:", ang1

    if verbose:
        print "axis", axis, "sub2 angle", ang1 * 180.0 / math.pi

    # get unitary translation along axix
    mintrans = 9e9
    for c1 in coms:
        for c2 in coms:
            if c1 is not c2:
                trans = abs(axis.dot(c2 - c1))
                if trans > 0.5:  # assume is < 0.5, is cyclic sym related subunit
                    mintrans = min(mintrans, trans)
    # numtrans = 0
    # tmptrans = 0
    # for c1 in coms:
    #   for c2 in coms:
    #       if c1 is not c2:
    #           trans = abs( axis.dot(c2-c1) )
    #           if trans - mintrans < 0.5:
    #               if verbose: print "trans",trans
    #               tmptrans += trans
    #               numtrans += 1
    # mintrans = tmptrans / numtrans

    # refine esimate min translation along axis by looking at all translations
    for i in range(4):
        tmptrans = 0
        numtrans = 0
        for xyz in coms[1:]:
            mult_actual = (xyz - coms[0]).dot(axis) / mintrans
            if abs(mult_actual - round(mult_actual)) < 0.1 and mult_actual > 0.5:
                # if verbose: print "ARST", mult_actual
                tmptrans += mintrans * mult_actual / round(mult_actual)
                numtrans += 1
            # elif abs(mult_actual-round(mult_actual))-0.5 < 0.1:
            #   # if verbose: print "ASDF", mult_actual
            #   tmptrans += mintrans * ( (mult_actual+0.5) / round( mult_actual+0.5 ) )
            #   numtrans += 1
                # if verbose: print mintrans * mult_actual / round( mult_actual
                # )
        if verbose:
            print "old mintrans", mintrans
        mintrans = tmptrans / numtrans
        if verbose:
            print "new mintrans", mintrans

    if verbose:
        print "mintrans along axis is", mintrans

    # check if min translation along axis is "unitary" translation, or a
    # multiple...
    unit_trans = 0
    for div in range(1, 4):
        allok = True
        if verbose:
            print "======", div, "======"
        for xyz in coms[1:]:
            mult_actual = (xyz - coms[0]).dot(axis) / mintrans
            if verbose:
                print mult_actual * div
            if abs(round(mult_actual * div) - mult_actual * div) > 0.05:
                allok = False
        if allok:
            break
    if not allok:
        if not allok:
            raise Exception("TROUBLE DETERMINING HELIX SPACING FROM COMs")
    mintrans /= div

    unit_trans = 0
    num_unit_trans = 0
    for xyz in coms[1:]:
        mult_actual = (xyz - coms[0]).dot(axis) / mintrans
        mult = round(mult_actual)
        if mult > 0:
            unit_trans += (xyz - coms[0]).dot(axis) / mult
            num_unit_trans += 1
            if verbose:
                print mult_actual, mult, (xyz - coms[0]).dot(axis), unit_trans, (xyz - coms[0]).dot(axis) / mult
    unit_trans /= num_unit_trans
    sub2_trans = (coms[1] - coms[0]).dot(axis)
    if verbose:
        print mintrans, unit_trans, sub2_trans, sub2_trans / unit_trans

    # check that each translation is even multiple of unit_trans
    for xyz in coms[1:]:
        test_mult = (xyz - coms[0]).dot(axis) / unit_trans
        if abs(test_mult - round(test_mult)) > 0.03:  # * test_mult:
            raise Exception(
                "(xyz-coms[0]).dot(axis) / unit_trans == %f" % test_mult)

    unit_ang1 = ang1 / (sub2_trans / unit_trans)
    unit_xform1 = rotation_around(axis, unit_ang1, cen)
    unit_xform1.t += unit_trans * axis

    # now recompute rotation angle based on last subunit
    # this corrects small errors from pdb coordinates or whatever...
    # or maybe unit_xform1.t += unit_trans*axis is somehow wrong and I'm
    # dumb...
    if verbose:
        print "========== TEST PRIMARY ROTATION", unit_ang1, "=========="
    correction_ang, nturns, error = get_correction_angle(
        axis, cen, coms, unit_xform1, error)
    if not correction_ang and not error:
        error = "correction_ang is None"
    if correction_ang:
        if verbose:
            print "correcting unit_ang1 and ang1"
        unit_ang1 = unit_ang1 - correction_ang / nturns
        ang1 = unit_ang1 * (sub2_trans / unit_trans)
    unit_xform1 = rotation_around(axis, unit_ang1, cen)
    unit_xform1.t += unit_trans * axis

    # angs        = [ang1]
    unit_angs = [unit_ang1]
    unit_xforms = [unit_xform1]
    for i in n_to_contacts:
        ang2 = 2.0 * i * math.pi + ang1
        # this is crappy duplicated code....
        unit_ang2 = ang2 / (sub2_trans / unit_trans)
        unit_xform2 = rotation_around(axis, unit_ang2, cen)
        unit_xform2.t += unit_trans * axis
        if verbose:
            print "========== GET GEOMETRY FOR PRIMARY ROTATION", unit_ang2, "=========="
        correction_ang2, nturns2, error2 = get_correction_angle(
            axis, cen, coms, unit_xform2, error)
        if not correction_ang2 and not error2:
            error = "correction_ang2 is None"
        if correction_ang2:
            if verbose:
                print "correcting unit_ang and ang"
            unit_ang2 = unit_ang2 - correction_ang2 / nturns2
            ang2 = unit_ang2 * (sub2_trans / unit_trans)
        unit_xform2 = rotation_around(axis, unit_ang2, cen)
        unit_xform2.t += unit_trans * axis
        # angs       .append(      ang2 )
        unit_angs  .append(unit_ang2)
        unit_xforms.append(unit_xform2)

    # now figure out Nfold based on coverage
    unit_xform = None
    unit_ang = None
    # ang = None
    if not error:
        error = "error determining helix symmetry (and check for ang errors)"

    # loop over nfolds, then over the unit_xform possibilities
    # accept first that covers all actual COMs from input structure
    # check up to nfold 6
    for nfold in nfolds:
        for i in range(len(unit_xforms)):
            if verbose:
                print "========== TEST COM COVERAGE FOR NFOLD", nfold, ", NANG", i
            unit_xform = unit_xforms[i]
            unit_ang = unit_angs[i]
            # ang = angs[i]
            if verbose:
                print "checking symmetry Nfold", nfold, "ang_option", i
            rms = xform_covers_all_coms(axis, cen, coms, unit_xform, nfold)
            if rms:
                error = None
                if verbose:
                    print "RMS to input coms is", rms
                break
        if not error:
            break
    if not error and verbose:
        print "nfold is", nfold

    # optional graphics
    if show > 0:
        cgo, cgo2 = [], []
        xyz = coms[0]
        for i in range(int(math.ceil(show / nfold))):
            syms = [RAD(axis, i * 360.0 / nfold, cen) for i in range(nfold)]
            for sym in syms:
                cgo.extend(cgo_sphere(sym * xyz))
            xyz = unit_xform * xyz
        for c in coms:
            cgo2.extend(cgo_sphere(c, col=(1, 0, 0)))
        cgo.extend(cgo_sphere(cen))
        cgo.extend(cgo_segment(cen, cen + axis * (xyz.dot(axis))))
        cmd.load_cgo(cgo, "determine_helix_geometry")
        cmd.load_cgo(cgo2, "existing_coms")

    # if verbose: print "unitary rotation angle", unit_ang*180.0/math.pi

    if error:
        raise Exception(error)

    return axis, cen, unit_ang, nfold, unit_trans, unit_xform


def make_chainAB(sele='vis', show=0, **args):
    axis, cen, unit_ang, nfold, unit_trans, unit_xform = determine_helix_geometry(
        sele, show=show, *args)
    cmd.remove('not chain A')
    cmd.create("chainB", "chain A")
    xform("chainB", unit_xform)
    cmd.alter('chainB', "chain='B'")


def make_chainAB_all_in_dir(pattern):
    cmd.delete("all")
    path = os.path.split(pattern)[0]
    for fn in glob.glob(pattern):
        newfn = fn
        while newfn.endswith(".gz"):
            newfn = newfn[:-3]
        while newfn.endswith(".pdb"):
            newfn = newfn[:-4]
        newfn += "_chainAB.pdb"
        if not os.path.exists(newfn):
            try:
                cmd.load(fn)
                make_chainAB(verbose=False)
                cmd.save(newfn)
            except Exception as e:
                print "========================= ERROR ============================"
                print fn
                print e
                print "================================++++========================"
            cmd.delete("all")


def make_helix_symdef(sele='vis', show=0):
    """
    # delete all; load /Users/sheffler/Downloads/N4_C1_DR53_001.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; make_helix_symdef(show=40)
    # delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N2_C1_DR04_065.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; make_helix_symdef(show=50)
    # delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N2_C3_DR04_003.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; make_helix_symdef(show=50)
    # delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N3_C1_DR04_002.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; make_helix_symdef(show=50)
    # delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N3_C2_DR04_001.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; make_helix_symdef(show=50)
    # delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N3_C4_DR04_001.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; make_helix_symdef(show=50)
    # delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N4_C3_DR04_001.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; make_helix_symdef(show=50)
    # delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N5_C1_DR04_002.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; make_helix_symdef(show=50)
    # delete all; load /Users/sheffler/Dropbox/project/hao_helix/hao_test/N4_C1_DR04_002.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; make_helix_symdef(show=50)


    # delete all; load /Users/sheffler/tmp/waiting_list/N2_C1_DR05_001.pdb; hide lin; show rib; util.cbc; run /Users/sheffler/pymol/helix_util.py; make_helix_symdef(show=50)

    """
    cgo = []

    axis, cen, unit_ang, nfold, unit_trans, unit_xform = determine_helix_geometry(
        sele, show=show)

    coms = [com("((%s) and name ca and chain %s)" % (sele, c))
            for c in cmd.get_chains(sele)]

    # now generate symfile with xyzs on axis in intervals linked to subunits of same height
    # will need to generate below and above

    maxtrans = max(c.dot(axis) for c in coms)

    hcom, hcom2 = coms[0], coms[0]
    for j in range(100):
        if j * unit_trans - 0.1 > maxtrans:
            break
        syms = [RAD(axis, i * 360.0 / nfold, cen) for i in range(nfold)]
        for sym in syms:
            a = cen + j * unit_trans * axis
            a2 = cen - j * unit_trans * axis
            x = axis
            cgo.extend(cgo_sphere(a))
            if j > 0:
                cgo.extend(cgo_sphere(a2))
            for i, in_com in enumerate(coms):
                if in_com.distance(sym * hcom) < 0.1:
                    y = sym * hcom - a
                    y2 = ~sym * hcom2 - a2
                    # print "xyz",
                    cgo.extend(cgo_sphere(a + y))
                    cgo.extend(cgo_segment(a, a + y))
                    if j > 0:
                        cgo.extend(cgo_sphere(a2 + y2))
                    if j > 0:
                        cgo.extend(cgo_segment(a2, a2 + y2))
                    break
        hcom = unit_xform * hcom
        hcom2 = unit_xform.inverse() * hcom2

    if show:
        cmd.load_cgo(cgo, "make_helix_symdef")
    print "make_helix_symdef done"


def makeh(sele='vis', n=30, objname='HELIX', nfolds=None, show=False):
    cmd.delete(objname)
    axis, cen, unit_ang, nfold, unit_trans, unit_xform = determine_helix_geometry(
        sele, nfolds=nfolds, show=show)
    print "================================== makeh ========================================"
    print "unit_nag:", unit_ang, "unit_trans:", unit_trans, "nfold:", nfold
    if n < 0:
        n = -n
        unit_xform = unit_xform.inverse()
        print unit_xform.pretty()

    n = n / nfold
    cmd.delete(objname)
    v = cmd.get_view()
    cmd.create('tmp', sele + ' and chain A and name n+ca+c')
    count = 0
    for nf in range(nfold):
        x = Xform()
        print nf, x.pretty()
        symrotang = 360.0 / nfold * nf
        for i in range(n):
            cmd.create('Htmp%i' % count, 'tmp')
            xform('Htmp%i' % count, x)
            trans('Htmp%i' % count, -cen)
            rot('Htmp%i' % count, axis, symrotang)
            trans('Htmp%i' % count, cen)
            cmd.alter('Htmp%i' % count, "chain='%s'" % ROSETTA_CHAINS[count])
            count += 1
            x = x * unit_xform
    cmd.create(objname, 'Htmp*')
    cmd.delete("Htmp*")
    cmd.delete('tmp')
    cmd.hide('ev', objname)
    cmd.show('lines', objname)
    util.cbc(objname)
    cmd.set_view(v)
    return axis, cen, unit_ang, nfold, unit_trans, unit_xform


cmd.extend('makeh', makeh)


def determine_pore_size(sele='vis', n=30, objname='HELIX', nfolds=None, show=False):
    axis, cen, unit_ang, nfold, unit_trans, unit_xform = makeh(
        sele, n, objname, nfolds, show)
    min_dis_to_axis = 9e9
    for a in cmd.get_model(objname + " and name CA").atom:
        xyz = Vec(a.coord)
        # project point to zero-plane perp to helix axis
        xyz = xyz - xyz.dot(axis) * axis
        dist = xyz.distance(cen)
        min_dis_to_axis = min(dist, min_dis_to_axis)
    print "min_dis_to_axis", min_dis_to_axis
    return min_dis_to_axis


def makeh_from_2chains(n=10, chain1='A', chain2='B'):
    x0 = getrelframe("chain " + chain1, "chain " + chain2)
    x = x0
    for i in range(n):
        cmd.create("tmp%i" % i, "not tmp* and chain " + chain2)
        xform("tmp%i" % i, x)
        x = x0 * x
        cmd.alter("tmp%i" % i, "chain='%s'" %
                  "CDEFGHIJKLMNOPQRSTUVWXYZ"[i % 24])


def do_cprofile(func):
    def profiled_func(*args, **kwargs):
        profile = cProfile.Profile()
        try:
            profile.enable()
            result = func(*args, **kwargs)
            profile.disable()
            return result
        finally:
            profile.print_stats()
    return profiled_func

#@do_cprofile


def test_all(loc, maxnum=99999999):
    failures = []
    with open("pymol_helix_test.log", 'w') as out:
        for fn in glob.glob(loc)[:maxnum]:
            try:
                cmd.delete("all")
                cmd.load(fn)
                cmd.hide("all")
                cmd.show("rib")
                util.cbc()
                make_helix_symdef(show=0)
                out.write(fn + "\n")
                out.flush()
            except Exception as e:
                out.write(fn + " / " + str(e) + "\n")
                out.flush()
                failures.append((fn, e))
        print "=================================== test_all DONE ==================================="
        for f, e in failures:
            print f, e


TEST_LIST = []
TEST_POS = -1


def helix_test(loc, advance=True):
    global TEST_LIST
    global TEST_POS
    if not TEST_LIST:
        assert loc
        for fn in glob.glob(loc):
            TEST_LIST.append(fn)
        TEST_POS = -1
    TEST_POS += 1
    current_fn = TEST_LIST[TEST_POS]

    print "============================================================================================================================="
    print "START"
    print TEST_POS
    print current_fn
    print "============================================================================================================================="

    cmd.delete("all")
    cmd.load(current_fn)
    cmd.hide("all")
    cmd.show("rib")
    util.cbc()
    make_helix_symdef(show=50)

    print "============================================================================================================================="
    print "FINISH"
    print TEST_POS
    print current_fn
    print "============================================================================================================================="


def load_tests(loader, tests, ignore):
    tests.addTests(doctest.DocTestSuite())
    return tests


if __name__ == '__main__':
    import doctest
    r = doctest.testmod()
    print r
