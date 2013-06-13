def sign(x):
	if x > 0: return  1.0
	if x < 0: return -1.0
	return 0.0

def gp120_fusion():
	cen = com('chain A+C+E')
	axis = cen.normalized()
	cmd.load("/Users/sheffler/Downloads/gp120_trimer_model.pdb","gp")
	alignc3('gp',tgtaxis=axis)
	cen2 = com('gp')
	l1 = cen.length()
	l2 = sign(cen.dot(cen2))*cen2.length()
	rot('gp',axis,random.random()*120.0)
	trans('gp',axis*( l1 - l2 + 100))
	cmd.remove("gp and not chain A")
	cmd.remove("not gp and not chain A+B")
	cmd.alter('gp','chain="C"')
	cmd.hide('ev')
	cmd.set('ribbon_radius',"7")
	cmd.set('ribbon_sampling',"3")
	cmd.show('rib')
	util.cbc()
	cmd.color("cyan" ,"chain A")
	cmd.color("green","chain B")
	cmd.color("yellow","chain C")
	makeicos("all",n=60)
	return axis

def gp120_fusion_oct():
	cen = com('chain B+J+L')
	axis = cen.normalized()
	cmd.load("/Users/sheffler/Downloads/gp120_trimer_model.pdb","gp")
	alignc3('gp',tgtaxis=axis)
	cen2 = com('gp')
	l1 = cen.length()
	l2 = sign(cen.dot(cen2))*cen2.length()
	rot('gp',axis,random.random()*120.0)
	trans('gp',axis*( l1 - l2 + 85))
	cmd.remove("gp and not chain A")
	cmd.remove("not gp and not chain A+B")
	cmd.alter('gp','chain="C"')
	cmd.hide('ev')
	cmd.set('ribbon_radius',"7")
	cmd.set('ribbon_sampling',"3")
	cmd.show('rib')
	util.cbc()
	cmd.color("yellow","chain C")
	makeoct("all",n=24)
	return axis