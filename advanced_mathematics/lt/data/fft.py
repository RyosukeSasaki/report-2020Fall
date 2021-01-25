#########################################################
#                                                       #
#   Author: Ryosuke Sasaki (rsasaki@keio.jp)            #
#   Date:   2021/1/26                                   #
#                                                       #
#########################################################

import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

peak_threshold = 0.1    #ピーク検出のしきい値
N = 128                 #データ数
dt = 0.00078125         #時間間隔

raw = np.loadtxt('./data_kadai2.csv', delimiter=',')

fk = np.fft.fft(raw[:, 1])                                      #Fkを計算
spec = np.vectorize(lambda x: np.abs(x))(fk)                    #スペクトル(Fkの絶対値)を計算
phase = np.vectorize(lambda x: np.angle(x))(fk)                 #位相を計算
amp = np.vectorize(lambda x: np.abs(x)*2/N)(fk)                 #振幅を計算
freq = np.fft.fftfreq(raw[:, 1].size, dt)               #インデックス毎の周波数を計算
peaks, _ = find_peaks(amp[0:int(N/2)], height=peak_threshold)   #ピークのインデックス求める

#
#周波数スペクトルを描画
#
fig0 = plt.figure()
ax = fig0.add_subplot(111)
ax.set_xlabel('frequency / Hz')
ax.set_ylabel('Fk')
lists = sorted(zip(*[freq, spec]))
freq_sorted, spec_sorted = list(zip(*lists))
fig0.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
ax.plot(freq_sorted, spec_sorted, c="g", linewidth=0.5)
ax.scatter(freq_sorted, spec_sorted, s=10, c="g")
fig0.savefig("Spectrum.png")

#
#振幅の散布図を描画
#
fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_xlabel('frequency / Hz')
ax.set_ylabel('Amplitude')
ax.scatter(freq[peaks], amp[peaks], s=10, label="Peaks", c="darkorange")
ax.scatter(freq[0:int(N / 2)], amp[0:int(N / 2)], s=3, label="Amplitude", c="g")
fig1.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
plt.legend(loc='best')
fig1.savefig("Amplitude.png")

#
#位相の散布図を描画
#
fig2 = plt.figure()
ax = fig2.add_subplot(111)
ax.set_xlabel('frequency / Hz')
ax.set_ylabel('Phase / degrees')
ax.scatter(freq[peaks], [np.degrees(x) for x in phase[peaks]], s=10, label="Peaks", c="darkorange")
ax.scatter(freq[0:int(N / 2)], [np.degrees(x) for x in phase[0:int(N / 2)]], s=3, label="Phase", c="g")
fig2.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
plt.legend(loc='best')
fig2.savefig("Phase.png")

print("Nyquist frequency:", abs(freq[int(N/2)]))
print("peaks:")
for i in peaks:
    print("Frequency:", freq[i], "  Amplitude:", amp[i], "  Phase(degrees):", np.degrees(phase[i]))
