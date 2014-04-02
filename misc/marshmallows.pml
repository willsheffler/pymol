remove hydro
bg_color white
hide everything
set surface_quality, 1
alter all, b=50
alter all, q=1
set gaussian_resolution,10
map_new mapA, gaussian, 1, chain A, 10
map_new mapB, gaussian, 1, chain B, 10
map_new mapC, gaussian, 1, chain C, 10
map_new mapE, gaussian, 1, chain E, 10
map_new mapG, gaussian, 1, chain G, 10
map_new mapJ, gaussian, 1, chain J, 10
map_new mapL, gaussian, 1, chain L, 10
isosurface surfA, mapA
isosurface surfB, mapB
isosurface surfC, mapC
isosurface surfE, mapE
isosurface surfG, mapG
isosurface surfJ, mapJ
isosurface surfL, mapL
color green, surfA
color palegreen, surfC
color palegreen, surfE
color palegreen, surfG
color cyan, surfB
color palecyan, surfJ
color palecyan, surfL
set antialias, 2
set ray_trace_gain, 0.4
set ray_shadows, 0
set specular, 0
show surface, surf*
