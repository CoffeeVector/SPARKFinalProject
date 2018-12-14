#!/usr/bin/env python
# coding: utf-8
import numpy as np
# Import numerical python modules as "np"; this stops us from having to type "numpy" all the time

import matplotlib.pyplot as plt
# Import plotting modules
from scipy import optimize
import scipy.stats as sp

import os

def load_osc_csv(fn, offset=True):
    with open(fn) as f:
        f.readline() #The first line is headers, just skip
        parts = f.readline().split(',') #The second line includes the time increment data
        data = np.loadtxt(f, delimiter=',', usecols=np.arange(len(parts)-2), unpack=True)

    start = float(parts[-2]) #Second to last entry is the time start
    inc = float(parts[-1]) #Last entry is the time increment
    data[0] *= inc #Multiply the time axis by the increment

    if offset: data[0] += start

    return data

def sin_fit_func(x, A, f, p, y0):
    return y0 + A * np.sin(2*np.pi*(f * (x - p)))

def getOpt(t, V, p0):
    popt, pconv = optimize.curve_fit(sin_fit_func, t, V, p0)
    plt.plot(t, V, label='data')
    plt.plot(t, sin_fit_func(t, *p0), 'b', label='original guess')
    plt.plot(t, sin_fit_func(t, *popt), 'r', label='best fit')
    plt.legend()
    plt.show()
    return popt;

def getDeltaPhi(t, V, p0):
    popt = getOpt(t, V, p0)
    print(str(popt[0])+ 'sin(' + str(popt[1]) +'t' + '-' + str(popt[2])+')')
    deltaPhi = popt[2]
#     if popt[0] < 0:
#         deltaPhi = deltaPhi - 1/25000000/2
    return deltaPhi
#If the phase shift is ever greater than the period, we have to either increase the frequency or do something else to fix it.
def getTimeDelay(filename, filenamebackground):
    t, V1, V2 = load_osc_csv(filename, offset=True)
    _, _,background = load_osc_csv(filenamebackground,offset=True)
    plt.plot(t, background, label='background')
    plt.plot(t, V2, label='Data')
    plt.legend()
    plt.show()
    meanV1 = np.mean(V1)
    meanV2 = np.mean(V2)
    meanB = np.mean(background)
    ampV1 = 0
    ampV2 = 0
    ampB = 0
    if(meanV1 < 0):
        ampV1 = meanV1 - np.min(V1)
    else:
        ampV1 = np.max(V1) - meanV1
    if(meanV2 < 0):
        ampV2 = meanV2 - np.min(V2)
    else:
        ampV2 = np.max(V2) - meanV2
    if(meanB < 0):
        ampB = meanB - np.min(background)
    else:
        ampB = np.max(background) - meanB
    freq = 25000000
    p1= (ampV1,freq, 0, meanV1)
    p2= (ampV2,freq, 0, meanV2)
    p3= (ampB,freq, 0, meanB)
    phi1 = getDeltaPhi(t, V1, p1)
    phi2 = getDeltaPhi(t, V2, p2)
    bopt = getOpt(t, background, p3)
    L = []
    for vt in zip(t,V2):
        L.append(vt[1] - sin_fit_func(vt[0], *bopt))
    ampL = 0
    meanL = np.mean(L)
    if(meanL < 0):
        ampL = meanL - np.min(L)
    else:
        ampL = np.max(L) - meanL
    p4= (ampL,freq, 0, meanL)
    plt.plot(t, L, label='difference')
    try:
        phi3 = getDeltaPhi(t, L, p4)
    except:
        print("Fail!")
        print('deltaPhi: ' + str(phi2-phi1))
        print('\n\n\n')
        return phi2-phi1

    print('phi1: ' + str(phi1))
    print('phi2: ' + str(phi2))
    print('phi3: ' + str(phi3))
    if(phi3 < phi1):
        print('deltaPhi: ' + str((1/25000000 + phi3)-phi1))
        print('DID A PHASE SHIFT')
        print('\n\n\n')
        return (1/25000000 + phi3)-phi1;
    else:
        print('deltaPhi: ' + str(phi3-phi1))
        print('\n\n\n')
        return phi2-phi1

x = []
phi = []
dir = "cats4/"
for filename in os.listdir(dir + 'lm'):
    distance = int(filename.replace("mm.csv",''))/1000
    x.append(float(distance))
    print('filename: ' + filename)
    print("Distance: " + str(distance))
    phi.append(getTimeDelay(dir + 'lm/' + filename, dir + 'm/' + filename))

xf = np.linspace(0, 3.5e-8, 100)
plt.plot(phi, x, 'ob')
slope, intercept, r_value, p_value, std_err = sp.linregress(phi, x)
plt.plot(xf, (xf)*slope + intercept , 'b', label='best fit')
print('speed of light: ' + str(slope))
print('standard error: ' + str(std_err))
print('speed of light confidence interval: (' + str(slope - 2*std_err) + ', ' + str(slope + 2*std_err) + ')')
print("R Squared: " + str(r_value**2))
plt.show()
