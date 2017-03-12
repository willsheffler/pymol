import os

import pyrosetta
import rosetta

from xyzMath import *


def main():
    pyrosetta.init()
    print "foo"
    dirname = os.path.dirname(__file__)
    testfile = dirname + '/1pgx.pdb'
    pose = rosetta.core.import_pose.pose_from_file(testfile)
    print pose
    x = Xform(Mat(), Vec(10, 0, 0))
    for ir in range(1, pose.size() + 1):
        for ia in range(1, pose.residue(ir).natoms() + 1):
            aid = rosetta.core.id.AtomID(ia, ir)
            oldxyz = Vec(pose.xyz(aid))
            # manipulate the atom position here!!!!
            newxyz = x * oldxyz
            pose.set_xyz(aid, newxyz.to_rosetta())
    print 'dumping' + dirname + '/test_out.pdb'
    pose.dump_pdb(dirname + '/test_out.pdb')

if __name__ == '__main__':
    main()
