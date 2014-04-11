# run /Users/sheffler/pymol/una.py; make_d3oct("test*","o33*",depth=3)

def make_d3oct(d3,cage,cage_trimer_chain="A", depth=4, maxrad=9e9):
	print cmd.super("(("+cage+") and (chain "+cage_trimer_chain+"))" ,"(("+d3+") and (chain A))")
	zcagecen = com(cage+" and name ca").z
	print zcagecen
	#return
	x = alignvectors(Vec(1,1,1),Vec(1,-1,0),Vec(0,0,1),Vec(1,0,0))
	# print x * Vec(1,1,1), x*Vec(1,-1,0)
	# RAD(Ux,180), RAD(Uy,120),
	G=[ RAD(Ux,180), RAD(Uz,120), RAD(x*Vec(1,0,0),90,Vec(0,0,zcagecen)), RAD(x*Vec(1,1,0),180,Vec(0,0,zcagecen)), ]
	makesym(G,sele="(("+d3+") and ((chain A+B) and name CA))",depth=depth,maxrad=maxrad)
	cmd.show("sph","MAKESYM")
	# cmd.disable("all")
	cmd.enable("MAKESYM")

def make_d3tet(d3,cage,cage_trimer_chain="A", depth=4, maxrad=9e9):
	print cmd.super("(("+cage+") and (chain "+cage_trimer_chain+"))" ,"(("+d3+") and (chain A))")
	zcagecen = com(cage+" and name ca").z
	print zcagecen
	#return
	x = alignvectors(Vec(1,1,1),Vec(1,-1,0),Vec(0,0,1),Vec(1,0,0))
	# print x * Vec(1,1,1), x*Vec(1,-1,0)
	# RAD(Ux,180), RAD(Uy,120),
	G=[ RAD(Ux,180), RAD(Uz,120), RAD(x*Vec(1,0,0),180,Vec(0,0,zcagecen)), ]
	makesym(G,sele="(("+d3+") and ((chain A+B) and name CA))",depth=depth,maxrad=maxrad)
	cmd.show("sph","MAKESYM")
	# cmd.disable("all")
	cmd.enable("MAKESYM")
