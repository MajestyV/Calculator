import argparse
import codecs
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

def raman():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True)
    parser.add_argument('-t', '--threshold', default='500', required=False)
    parser.add_argument('-d', '--distance', default='10', required=False)
    # parser.add_argument('-w', '--wavelength scale', default='700', required=False)
    args = parser.parse_args()
    
    with codecs.open(args.input, 'rb', 'utf-8', 'ignore') as f:
        dataRaw = f.readlines()
        coef = dataRaw[:42]
        data = dataRaw[42:]
    
    x = []
    y = []
    for xy in data:
        x.append(float(xy.split()[0]))
        y.append(float(xy.split()[1]))
    ymax = max(y)
    xmax = x[y.index(max(y))]
    
    plt.figure()
    plt.xlabel('Raman Shift(/cm)')
    plt.ylabel('Intensity(a.u.)')
    plt.plot(x, y, label=args.input.split('.')[0])
    plt.xticks(np.arange(round(min(x)), round(max(x)), 50.0))
    # plt.yticks([]) # hide y-axis
    peaks, _ = find_peaks(y, height=float(args.threshold), distance=int(args.distance))
    print(peaks)
    for peak in list(peaks):  # peaks: return original index of y
        # print(x[peak],y[peak])
        plt.plot(x[peak], y[peak], 'ro')
        plt.annotate(round(x[peak]), xytext=(x[peak] + 2, y[peak]), xy=(x[peak], y[peak]))
    plt.plot(xmax, ymax, 'ro')  # Guarantee at least one maximum
    plt.annotate(round(xmax), xytext=(xmax + 2, ymax), xy=(xmax, ymax))
    plt.legend(loc='upper left')
    plt.show()

if __name__=="__main__":
    raman()