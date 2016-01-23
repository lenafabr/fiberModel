
hide all

select bb, name A1+A2+F1+F2
select cent, name A3+F3
select prot, name AC
select filler, name F1+F2+F3

show lines, bb
set stick_radius=10
show sticks, cent

# display nucleosomes as cylinders
#import displayCylinders
#displayCylinders.makeCylinders()

spectrum count, rainbow
color gray50, bb
#color wheat, prot

# display nucleosomes as cylinders
run scripts/displayCylinders.py
makeCylinders AC, AC1, AC2, 52, 55
hide bb

# showing flexible tails
select tails, name FT
show spheres, tails
color slate, tails
set sphere_scale=3, tails
show lines, tails

# show charge positions
select charges, name CH
show spheres, charges
color slate, charges
set sphere_scale=3, charges
