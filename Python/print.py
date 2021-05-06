
# x_i = 1.6
# x_f = 2.2
# nbr_points = 20
# increment = (x_f - x_i) / nbr_points
numbers = [ 1.0000,
            1.1000,
            1.2000,
            1.3000,
            1.4000,
            1.5000,
            1.6000,
            1.7000,
            1.7500,
            1.8000,
            1.8500,
            1.9000,
            1.9500,
            2.0000,
            2.0500,
            2.1000,
            2.1500,
            2.2000,
            2.3000,
            2.4000,
            2.5000,
            2.6000,
            2.7000,
            2.8000,
            2.9000,
            3.0000,
            3.1000,
            3.2000,
            3.3000,
            3.4000,
            3.5000]

for x in range(len(numbers)):
    print("/generator/momentum/SetMomentum " + str(round(numbers[x], 5)))
    print("/analysis/setFileName /media/matthieu/ssd1/Geant4/Data/DataSet_4/test_" + str(x))
    print("/run/beamOn 100000")

