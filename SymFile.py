import sys
from xyzMath import Vec


class VirtXYZ(object):
    """
    >>> v1 = VirtXYZ( "test1", Vec(1,0,0), Vec(0,1,0), Vec(0,0,0) )
    >>> v2 = VirtXYZ( "test2", Vec(-1,0,0), Vec(0,1,0), Vec(0,0,0) )
    >>> print v1.symline()
    xyz test1                +1.000000,+0.000000,+0.000000  +0.000000,+1.000000,+0.000000  +0.000000,+0.000000,+0.000000
    >>> print v2.symline()
    xyz test2                -1.000000,+0.000000,+0.000000  +0.000000,+1.000000,+0.000000  +0.000000,+0.000000,+0.000000
    """

    def __init__(self, name, x, y, cen):
        super(VirtXYZ, self).__init__()
        self.name = name
        self.x = x
        self.y = y
        self.cen = cen

    def symline(self):
        s = "xyz %s %+f,%+f,%+f  %+f,%+f,%+f  %+f,%+f,%+f" % (self.name.ljust(
            20), self.x.x, self.x.y, self.x.z, self.y.x, self.y.y, self.y.z, self.cen.x, self.cen.y, self.cen.z)
        return s.rstrip()


class ConnectVirtual(object):
    """
    >>> v1 = VirtXYZ( "test1", Vec(1,0,0), Vec(0,1,0), Vec(0,0,0) )
    >>> v2 = VirtXYZ( "test2", Vec(-1,0,0), Vec(0,1,0), Vec(0,0,0) )
    >>> print ConnectVirtual( "test_jump", v1, v2 ).symline()
    connect_virtual            test_jump test1                test2
    """

    def __init__(self, name, virt1, virt2=None, chain=None):
        super(ConnectVirtual, self).__init__()
        self.name = name
        self.virt1 = virt1
        self.virt2 = virt2
        self.chain = chain

    def symline(self):
        if self.chain:
            s = "connect_virtual %s %s %s" % (self.name.rjust(
                20), self.virt1.name.ljust(20), "SUBUNIT", chain)
        else:
            assert self.virt2
            s = "connect_virtual %s %s %s" % (self.name.rjust(
                20), self.virt1.name.ljust(20), self.virt2.name.ljust(20))
        return s.rstrip()


if __name__ == '__main__':
    import doctest
    r = doctest.testmod()
    print r


# some not-finished work making a symdef parsers

# from pyparsing import Dict,Group,Literal,Regex,Word,ZeroOrMore,OneOrMore,alphas,delimitedList,nums,LineStart,LineEnd,Optional,Keyword,NotAny,StringEnd

# def printf(*args):
# 	for i,a in enumerate(args): print i,a

# paren   = Literal("(")|Literal(")")
# tag     = Regex(r'[a-zA-Z0-9_]+')
# num   = Word(nums).setParseAction( lambda s,l,t: int(t[0]) )
# real    = Regex(r'-?\d+(\.\d*)?([eE]\d+)?')
# vec     = real+","+real+","+real
# sname   = Keyword("symmetry_name" )+tag
# anchor  = Keyword("anchor_residue")+tag
# efirst  = (real+"*"+tag)
# erest   = ("+"+real+"*"+paren+tag+":"+tag+paren)
# eline   = (Keyword("E")+"=" +efirst+OneOrMore(erest))
# vcoord  = (Keyword('xyz') + tag + vec + vec + vec)
# vcoords = Keyword("virtual_coordinates_start") + OneOrMore(vcoord) + Keyword("virtual_coordinates_stop")
# subunit = Keyword("SUBUNIT") + Optional(Regex("[a-zA-Z0-9_]")) + Optional(tag|num)
# connect = Keyword("connect_virtual")+tag+tag+(NotAny("SUBUNIT")+tag|subunit)

# doftype = Keyword("x")|Keyword("y")|Keyword("z")|Keyword("angle_x")|Keyword("angle_y")|Keyword("angle_z")
# dof     = (doftype+Optional(paren+real+paren))
# dofs    = Keyword("set_dof")+tag+OneOrMore(dof   )

# jgroup  = Keyword("set_jump_group")+tag+OneOrMore(tag)
# symline = (sname|eline|anchor|vcoords|connect|dofs  |jgroup)+LineEnd()
# symfile = ZeroOrMore(symline)+StringEnd()

# real   .setParseAction( lambda s,l,t: float(t[0]) )
# sname  .setParseAction(tuple)
# anchor .setParseAction(tuple)
# efirst .setParseAction(lambda s,l,t: (t[0],t[2]))
# erest  .setParseAction(lambda s,l,t: ( t[1],t[4],t[6]))
# eline  .setParseAction(lambda s,l,t: ('E',tuple(t[2:])))
# vcoord .setParseAction(lambda s,l,t: tuple(t[1:]))
# vcoords.setParseAction(lambda s,l,t: ('xyz',tuple(t[1:-1])))
# subunit.setParseAction(lambda s,l,t: ("subunit",t[1:]))
# connect.setParseAction(lambda s,l,t: ("connect",tuple([x for x in t[1:] if not x=='\n'])))
# dof    .setParseAction(lambda s,l,t: ('dof',(t[0],t[2]) if len(t)>2 else (t[0],)))
# dofs   .setParseAction(lambda s,l,t: ('dofs',tuple(t[1:])))

# jgroup .setParseAction(lambda s,l,t: ('jump_group',t[1:]))
# symline.setParseAction(lambda s,l,t: tuple([x for x in t if not x=='\n']))

# print dof   .parseString("x(20.4244497891662)")
# print dof   .parseString("angle_x")
# print dofs  .parseString("set_dof JUMP0_0_to_com x(20.4244497891662)")
# #print dofs  .parseString("set_dof JUMP0_0_to_subunit angle_x angle_y angle_z")
# #print dofs  .parseString("set_dof JUMP0_0 x(3.81281197834553) angle_x")
# # print connect.parseString("connect_virtual JUMP0_0_to_com VRT0_0 VRT0_0_base")
# # print connect.parseString("connect_virtual JUMP0_0_to_subunit TEST SUBUNIT")
# # print connect.parseString("connect_virtual JUMP1_0_to_subunit VRT1_0_base SUBUNIT A")
# # print connect.parseString("connect_virtual JUMP1_0_to_subunit VRT1_0_base SUBUNIT A 8")
# sys.exit()
# s = "\n".join(open("/data/pdb/D3/1WAAA_D3.symm").readlines())
# for p in symline.searchString(s):
# 	print p


# testData = """
# +-------+------+------+------+------+------+------+------+------+
# |       |  A1  |  B1  |  C1  |  D1  |  A2  |  B2  |  C2  |  D2  |
# +=======+======+======+======+======+======+======+======+======+
# | min   |   7  |  43  |   7  |  15  |  82  |  98  |   1  |  37  |
# | max   |  11  |  52  |  10  |  17  |  85  | 112  |   4  |  39  |
# | ave   |   9  |  47  |   8  |  16  |  84  | 106  |   3  |  38  |
# | sdev  |   1  |   3  |   1  |   1  |   1  |   3  |   1  |   1  |
# +-------+------+------+------+------+------+------+------+------+
# """
# underline = Word("-=")
# number = Word(nums).setParseAction( lambda s,l,t:(l,[int(t[0])]) )
# vert = Literal("|").suess()

# rowDelim = ("+" + ZeroOrMore( underline + "+" ) ).suess()
# columnHeader = Group(vert + vert + delimitedList(Word(alphas + nums), "|") + vert)

# heading = rowDelim + columnHeader.setResultsName("columns") + rowDelim
# rowData = Group( vert + Word(alphas) + vert + delimitedList(number,"|") + vert )
# trailing = rowDelim

# datatable = heading + Dict( ZeroOrMore(rowData) ) + trailing

# for i in datatable.parseString(testData): print i
