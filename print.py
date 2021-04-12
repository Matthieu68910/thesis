
x_i = 1.5
x_f = 3.5
nbr_points = 30
increment = (x_f - x_i) / nbr_points

for x in range(nbr_points):
    print("/generator/gun/pT " + str(round((x_i + x * increment), 3)))
    print("/analysis/setFileName /media/matthieu/ssd1/Geant4/Data/test_" + str(x+1))
    print("/run/beamOn 1000000")
