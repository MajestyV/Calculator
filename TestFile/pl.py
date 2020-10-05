import argparse
import codecs
import numpy as np
import matplotlib.pyplot as plt

def plot():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True)
    # parser.add_argument('-w', '--wavelength scale', default='700', required=False)
    args = parser.parse_args()
    
    with codecs.open(args.input, 'rb', 'utf-8', 'ignore') as f:
        dataRaw = f.readlines()
        coef = dataRaw[:42]
        data = dataRaw[42:]
    
    x=[]
    y=[]
    for xy in data:
        x.append(float(xy.split()[0]))
        y.append(float(xy.split()[1]))
    ymax = max(y)
    xmax = x[y.index(max(y))]
    
    plt.figure()
    plt.xlabel('Wavelength(nm)')
    plt.ylabel('Intensity(a.u.)')
    plt.plot(x,y,label=args.input.split('.')[0])
    plt.xticks(np.arange(min(x), max(x), 20.0))
    plt.yticks([]) #hide y-axis
    plt.plot(xmax, ymax, 'ro')
    plt.annotate(round(xmax), xytext=(xmax+2, ymax), xy=(xmax, ymax))
    plt.legend(loc='best')
    plt.show()

if __name__=="__main__":
    plot()