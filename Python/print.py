
# x_i = 1.6
# x_f = 2.2
# nbr_points = 20
# increment = (x_f - x_i) / nbr_points
numbers = [ 1.32089,
            1.48578,
            1.57853,
            1.71422,
            1.79152,
            1.87568,
            1.93064,
            1.98904,
            2.13847,
            2.42703,
            2.51978,
            2.80491,
            2.97495,
            3.16560,
            1.20000,
            1.30000,
            1.40000,
            1.50000,
            1.60000,
            1.70000,
            1.80000,
            1.90000,
            2.00000,
            2.10000,
            2.20000,
            2.30000,
            2.40000,
            2.50000,
            2.60000,
            2.70000,
            2.80000,
            2.90000,
            3.00000,
            3.10000,
            3.20000,
            1.75000,
            1.85000,
            1.95000,
            2.05000,
            2.15000,
            2.25000]

numbers.sort()

for x in range(len(numbers)):
    print("/generator/momentum/SetMomentum " + str(round(numbers[x], 5)))
    print("/analysis/setFileName /media/matthieu/ssd1/Geant4/Data/Data_figure20/data_" + str(x))
    print("/run/beamOn 10000")

