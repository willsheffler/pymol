from pymol import cmd



def marshmallows(sel):
	cmd.delete("surf*")
	cmd.delete("map*")
	COLs = ("green","cyan","magenta","yellow","pink")
	COLs = COLs + COLs + COLs
	cmd.do("remove hydro")
	cmd.do("bg_color white")
	cmd.do("hide everything")
	cmd.do("set surface_quality, 1")
	cmd.do("alter all, b=50")
	cmd.do("alter all, q=1")
	cmd.do("set gaussian_resolution,10")

	for i,c in enumerate(cmd.get_chains(sel)):
		cmd.do("map_new map%s, gaussian, 1, (%s and chain %s), 10"%(c,sel,c))
		cmd.do("isosurface surf%s, map%s"%(c,c))
		cmd.do("color %s, surf%s"%(COLs[i%2],c))

	cmd.do("set antialias, 2")
	cmd.do("set ray_trace_gain, 0.4")
	cmd.do("set ray_shadows, 0")
	cmd.do("set specular, 0")
	cmd.do("show surface, surf*")
