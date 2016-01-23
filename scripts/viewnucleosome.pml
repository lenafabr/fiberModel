hide all
show cartoon
select axtrace, name AX
select ends, name XC+XT+XB

show sticks, axtrace
show sticks, ends

run scripts/displayCylinders.py
makeCylinders AC, AC1, AC2, 52, 55, 0.3

color magenta, ends
show spheres, ends

# visualize disco charges if present
select pcharge, name PC
select mcharge, name MC
color slate, mcharge
color red, pcharge
show spheres, pcharge
show spheres, mcharge
set sphere_scale=2, mcharge

select tail, name FT
color violet, tail
show spheres, tail
set sphere_scale=1, tail
show sticks, tail