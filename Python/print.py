
# x_i = 1.6
# x_f = 2.2
# nbr_points = 20
# increment = (x_f - x_i) / nbr_points
numbers = [0.00000,39.19081,19.59615,13.06493,9.79957,7.84055,6.53470,5.60210,4.90277,4.35896,3.92401,3.56823,3.27184,3.02112,2.80628,2.62016,2.45737,2.31379,2.18622,2.07213,1.96950,1.87669,1.79237,1.71542,1.64493,1.58012,1.52033,1.46501,1.41368,1.36592,1.32139]

for x in range(len(numbers)):
    print("/generator/momentum/SetMomentum " + str(round(numbers[x], 5)))
    print("/analysis/setFileName /media/matthieu/ssd1/Geant4/Data/test_" + str(x))
    print("/run/beamOn 500000")

