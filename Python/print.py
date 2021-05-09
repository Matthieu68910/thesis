
# x_i = 1.6
# x_f = 2.2
# nbr_points = 20
# increment = (x_f - x_i) / nbr_points
numbers = [ 15.074491,
            13.064930,
            12.248565,
            10.887989,
            9.799568,
            7.260152,
            6.324103,
            6.324103,
            5.160402,
            4.705516,
            4.264396,
            3.166686,
            2.975539,
            2.806284,
            2.519976,
            2.427232,
            2.139088,
            1.989192,
            1.931279,
            1.876692,
            1.792368,
            1.715421,
            1.580117,
            1.486635,
            1.321387]

#numbers.sort()

for x in range(len(numbers)):
    print("/generator/momentum/SetMomentum " + str(round(numbers[x], 5)))
    print("/analysis/setFileName /media/matthieu/ssd1/Geant4/Data/Data_figure17-19-newCC/data_" + str(x))
    print("/run/beamOn 10000")

