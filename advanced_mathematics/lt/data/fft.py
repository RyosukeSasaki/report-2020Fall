#########################################################
#                                                       #
#   Author: Ryosuke Sasaki (rsasaki@keio.jp)            #
#   Date:   2021/1/26                                   #
#                                                       #
#########################################################

import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import math

height_threshold = 0.1  # ピーク検出のしきい値
lower_threshold = 0.6  # 連続したピークを検出するための下限
N = 128  # データ数
dt = 0.0078125  # 時間間隔
num_sample = 1000  # 再合成のサンプル数
outdir = "kadai3"  # 出力ディレクトリ
name = "kadai3"  # ファイル名

#
# データの格納
#
raw = np.loadtxt('./data_'+name+'.csv', delimiter=',')

#
# FFT
#
fk = np.fft.fft(raw[:, 1])  # Fkを計算
spec = np.vectorize(lambda x: np.abs(x))(fk)  # スペクトルを計算
phase = np.vectorize(lambda x: np.angle(x))(fk)  # 位相を計算
amp = np.vectorize(lambda x: np.abs(x)*2/N)(fk)  # 振幅を計算
freq = np.fft.fftfreq(raw[:, 1].size, dt)  # index毎の周波数を計算

# 極大値によるピークの検出
peaks, _ = find_peaks(amp[0:int(N / 2)], height=height_threshold)
# 連続したピークの検出
peaks_ = [int(i) for i in range(int(N/2)) if amp[i]>=0.6 and not i in peaks]
if len(peaks_) > 0:
    peaks = np.append(peaks, peaks_)


#
# 波形再合成
#


def f_resyn(x):
    ret = 0
    for i in peaks:
        ret += amp[i] * math.cos(2 * math.pi * freq[i] * x + phase[i])

    return ret


time = [dt*N/num_sample*i for i in range(num_sample)]
f = [f_resyn(dt*N/num_sample*i) for i in range(num_sample)]

#
# 周波数スペクトルを描画
#
fig0 = plt.figure()
ax = fig0.add_subplot(111)
ax.set_xlabel('frequency / Hz')
ax.set_ylabel('$|F_k|$')
lists = sorted(zip(*[freq, spec]))
freq_sorted, spec_sorted = list(zip(*lists))
fig0.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
ax.plot(freq_sorted, spec_sorted, c="g", linewidth=0.5)
ax.scatter(freq_sorted, spec_sorted, s=10, c="g")
fig0.savefig(outdir+"/Spectrum-"+name+".png")

#
# 振幅の散布図を描画
#
fig1 = plt.figure()
ax = fig1.add_subplot(111)
ax.set_xlabel('frequency / Hz')
ax.set_ylabel('Amplitude')
ax.scatter(freq[peaks], amp[peaks], s=10, label="Peaks", c="b")
ax.scatter(freq[0:int(N / 2)], amp[0:int(N / 2)],
           s=3, label="Amplitude", c="g")
ax.plot(freq[0:int(N / 2)], amp[0:int(N / 2)], c="g", linewidth=0.5)
fig1.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
plt.legend(loc='best')
fig1.savefig(outdir+"/Amplitude-"+name+".png")

#
# 位相の散布図を描画
#
fig2 = plt.figure()
ax = fig2.add_subplot(111)
ax.set_xlabel('frequency / Hz')
ax.set_ylabel('Phase / degrees')
ax.scatter(freq[peaks], [np.degrees(x)
                         for x in phase[peaks]], s=10, label="Peaks", c="b")
ax.scatter(freq[0:int(N / 2)], [np.degrees(x)
                                for x in phase[0:int(N / 2)]], s=3, label="Phase", c="g")
fig2.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
plt.legend(loc='best')
fig2.savefig(outdir+"/Phase-"+name+".png")

#
# 入力データ,再合成波形の描画
#
fig3 = plt.figure()
ax = fig3.add_subplot(111)
ax.set_xlabel("Time / sec")
ax.set_ylabel("Amplitude")
ax.plot(raw[:, 0], raw[:, 1], c="g", linewidth=0.5)
ax.scatter(raw[:, 0], raw[:, 1], c="g", s=3, label="Raw data")
ax.plot(time, f, c="b", linewidth=0.5, label="Resynthesis wave")
fig3.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
plt.legend(loc='best')
fig3.savefig(outdir+"/Wave-"+name+".png")

#
# 入力データの描画
#
fig4 = plt.figure()
ax = fig4.add_subplot(111)
ax.set_xlabel("Time / sec")
ax.set_ylabel("Amplitude")
ax.plot(raw[:, 0], raw[:, 1], c="g", linewidth=0.5)
ax.scatter(raw[:, 0], raw[:, 1], c="g", s=3, label="Raw data")
fig4.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
plt.legend(loc='best')
fig4.savefig(outdir+"/WaveOrigin-"+name+".png")

#
# 結果の出力
#
print("Nyquist frequency:", abs(freq[int(N/2)]))
print("peaks:")
for i in peaks:
    print("Frequency:", freq[i], "  Amplitude:",
          amp[i], "  Phase(degrees):", np.degrees(phase[i]))

with open(outdir + "/" + name + ".txt", 'w') as f:
    print("Nyquist frequency:", abs(freq[int(N / 2)]), file=f)
    print("peaks:", file=f)
    for i in peaks:
        print("Frequency:", freq[i], "  Amplitude:", amp[i],
              "  Phase(degrees):", np.degrees(phase[i]), file=f)
