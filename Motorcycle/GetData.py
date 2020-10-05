#import numpy as np
import matplotlib.pyplot as plt
import codecs
import re
from os import path

class geda:
    # default_path =

    def __init__(self):
        self.name = geda

    #def SmoothingCurve(self,):

    def GetRaman(self,DataFile):
        test_parameter = []
        data_list = []
        wavelength = []
        intensity = []

        f = codecs.open(DataFile, 'rb', 'utf-8', 'ignore')  # Open file, using codecs to uniform coding type
        line = f.readline()                                 # Read the first line
        while line:
            if re.search("#",line):                         # Searching for pattern "#"
                test_parameter.append(line)
            else:
                data = list(map(float,line.split()))        # Reading data and split them apart, every pair of data will be a list with [wavelength,intensity] form
                data_list.append(data)                      # Saving the whole list for potential use
                wavelength.append(data[0])
                intensity.append(data[1])
            line = f.readline()                             # Read the next line before ending the loop, keeping the loop functioning until finishing reading the whole file
        f.close()                                           # Closing file

        return wavelength,intensity,test_parameter,data_list

    def MapRaman(self,DataFile):
        x = self.GetRaman(DataFile)[0]
        y = self.GetRaman(DataFile)[1]
        plt.plot(x,y)
        plt.show()
        return


gd = geda()

file_path = path.dirname(path.dirname(__file__))+'/TestFile/240 50%-1.txt'
#file_path = '/TestFile/210 10%-1.txt'
#print(gd.GetRaman(file_path)[2])
#gd.MapRaman(file_path)

