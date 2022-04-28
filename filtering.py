import numpy as np
import pandas as pd
from scipy import signal
import matplotlib.pyplot as plt
from scipy.signal import welch

#https://stackoverflow.com/questions/33973717/computing-a-power-spectrum

def sine_generator(fs, sinefreq, duration):
    T = duration
    nsamples = T/fs
    w = 2. * np.pi * sinefreq
    t_sine = np.linspace(0, T, int(nsamples), endpoint=False)
    y_sine = np.sin(w * t_sine)
    return y_sine

def butter_highpass(cutoff, fs, order=5):
    # N the order of the filter
    # Wn is the critical frequency 
    # btype can be lowpass, highpass, bandpass, bandstop
    # output can be ba, zpk, or sos. Apparently sos should be used
    # for general-purpose filtering
    # fs is the sampling frequency of the digital system
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='high', analog=False)
    return b, a

def butter_highpass_filter(data, cutoff, fs, order=5):
    b, a = butter_highpass(cutoff, fs, order=order)
    y = signal.filtfilt(b, a, data)
    return y


# Yearly
sine_yr = 1/(12*30*24*3600) 
duration = 5*12*24*30*3600 #5 years
fps = 24*3600 # this is going to be monthly 
sine_yearly = sine_generator(fps,sine_yr,duration)

# Seasonaly
sine_sea = 1/(3*30*24*3600) 
duration = 5*12*24*30*3600 #5 years
fps = 24*3600 # this is going to be monthly 
sine_season = sine_generator(fps,sine_sea,duration)

# monthly
sine_mon = 1/(30*24*3600) #monthly frequency
duration = 5*12*24*30*3600 #5 years
fps = 24*3600 # this is going to daily
sine_monthly = sine_generator(fps,sine_mon,duration)

sine = sine_yearly + sine_monthly + sine_season

f, psd = welch(sine,
               fs=fps,  # sample rate
               window='hanning',   # apply a Hanning window before taking the DFT
               nperseg=360,        # compute periodograms of 256-long segments of x
               detrend='constant') # detrend x by subtracting the mean

plot(psd)

f, psd = welch(sine,
               fs=1/fps,  # sample rate
               window='hanning',   # apply a Hanning window before taking the DFT
               nperseg=1200,        # compute periodograms of 256-long segments of x
               noverlap=0,
               detrend=False,
                average='mean') # detrend x by subtracting the mean

figure()
loglog(f,psd)

#filtered_sine = butter_highpass_filter(sine,sine_fq,fps)
# I need to remember that yearly signal has a smaller frequency
# than monthly signal, meaning that, if I want to keep monthly 
# frequency I need to do btype='high'

# I think I got it for yearly
b,a = signal.butter(5,1/(6*30*24*3600),btype='low',fs=1/(24*3600))
y = signal.filtfilt(b, a, sine)
plot(y)

# for monthly
b,a = signal.butter(5,1/(7*30*24*3600),btype='high',fs=1/(24*3600))
y = signal.filtfilt(b, a, sine)
#plot(y)

b,a = signal.butter(5,1/(1.5*30*24*3600),btype='high',fs=1/(24*3600))
yy = signal.filtfilt(b, a, y)
plot(yy)

# for seasonally
b,a = signal.butter(5,1/(7*30*24*3600),btype='high',fs=1/(24*3600))
y = signal.filtfilt(b, a, sine)
#plot(y)

b,a = signal.butter(5,1/(1.5*30*24*3600),btype='low',fs=1/(24*3600))
yy = signal.filtfilt(b, a, y)

plot(yy)

