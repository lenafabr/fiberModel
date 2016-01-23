
hide all

select bb, name A1+A2+F1+F2
select cent, name A3+F3
select prot, name AC
select filler, name F1+F2+F3


show lines, bb
set stick_radius=2
show sticks, cent

color gray50, bb
color forest, name A3
color skyblue, name F3

# display nucleosomes as cylinders
run scripts/displayCylinders.py
makeCylinders AC, AC1, AC2, 42, 45, 0.5
