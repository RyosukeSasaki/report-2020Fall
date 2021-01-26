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
num_sample = 1000  # 再合成のサンプル数
resample = 1  # ダウンサンプリングの回数
outdir = "kadai3-1"  # 出力ディレクトリ
name = "kadai3"  # ファイル名


#
# データの格納
#
raw = np.loadtxt('./data_' + name + '.csv', delimiter=',')
raw_x = raw[:, 0][::resample]
raw_y = raw[:, 1][::resample]
dt = raw_x[1] - raw_x[0]  # 時間間隔
N = len(raw_x)  # ダウンサンプル後のデータ数
N_raw = len(raw[:, 0])  # 元データ数

#
# FFT
#
fk = np.fft.fft(raw_y)  # Fkを計算
spec = np.vectorize(lambda x: np.abs(x))(fk)  # スペクトルを計算
phase = np.vectorize(lambda x: np.angle(x))(fk)  # 位相を計算
amp = np.vectorize(lambda x: np.abs(x)*2/N)(fk)  # 振幅を計算
freq = np.fft.fftfreq(raw_y.size, dt)  # index毎の周波数を計算
F_Nyq = abs(freq[int(N/2)])  # ナイキスト周波数

fk_raw = np.fft.fft(raw[:, 1])  # ダウンサンプル前のFk
spec_raw = np.vectorize(lambda x: np.abs(x))(fk_raw)
amp_raw = np.vectorize(lambda x: np.abs(x)*2/N_raw)(fk_raw)
freq_raw = np.fft.fftfreq(raw[:, 1].size, raw[1, 0] - raw[0, 0])

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
ax.scatter(freq[0:int(N / 2)], amp[0:int(N / 2)], s=3, label="Amplitude", c="g")
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
ax.scatter(freq[peaks], [np.degrees(x) for x in phase[peaks]], s=10, label="Peaks", c="b")
ax.scatter(freq[0:int(N / 2)], [np.degrees(x) for x in phase[0:int(N / 2)]], s=3, label="Phase", c="g")
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
ax.plot(raw_x, raw_y, c="g", linewidth=0.5)
ax.scatter(raw_x, raw_y, c="g", s=3, label="Raw data")
ax.plot(time, f, c="b", linewidth=0.5, label="Resynthesis wave")
fig3.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
#plt.legend(loc='best')
plt.legend(loc='upper right')
fig3.savefig(outdir+"/Wave-"+name+".png")

#
# 入力データの描画
#
fig4 = plt.figure()
ax = fig4.add_subplot(111)
ax.set_xlabel("Time / sec")
ax.set_ylabel("Amplitude")
ax.plot(raw_x, raw_y, c="g", linewidth=0.5)
ax.scatter(raw_x, raw_y, c="g", s=3, label="Raw data")
fig4.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
plt.legend(loc='best')
fig4.savefig(outdir+"/WaveOrigin-"+name+".png")

#
# 周波数スペクトル(比較用)
#
fig5 = plt.figure()
ax = fig5.add_subplot(111)
ax.set_xlabel('frequency / Hz')
ax.set_ylabel('$|X_k|$')
plt.ylim(0, 70)
fig5.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)

ax.annotate('', xy=[F_Nyq, 0], xytext=[F_Nyq, 3.5],
            arrowprops=dict(shrink=0, width=0.7, headwidth=4,
                            headlength=7, connectionstyle='arc3',
                            facecolor='gray', edgecolor='gray'))
ax.annotate('', xy=[-F_Nyq, 0], xytext=[-F_Nyq, 3.5],
            arrowprops=dict(shrink=0, width=0.7, headwidth=4,
                            headlength=7, connectionstyle='arc3',
                            facecolor='gray', edgecolor='gray'))
ax.text(F_Nyq, 4.2, "$F_N$", size=10, horizontalalignment="center")
ax.text(-F_Nyq, 4.2, "$-F_N$", size=10, horizontalalignment="center")

lists = sorted(zip(*[freq_raw, spec_raw]))
freq_raw_sorted, spec_raw_sorted = list(zip(*lists))
ax.plot(freq_raw_sorted, spec_raw_sorted, c="b", linewidth=0.5)
ax.scatter(freq_raw_sorted, spec_raw_sorted, s=10, c="b", label="Original")

if resample > 1:
    lists = sorted(zip(*[freq, spec]))
    freq_sorted, spec_sorted = list(zip(*lists))
    ax.plot(freq_sorted, spec_sorted, c="g", linewidth=0.5)
    ax.scatter(freq_sorted, spec_sorted, s=10, c="g", label="Downsampled")
    plt.legend(loc='best')

fig5.savefig(outdir+"/Spectrums-"+name+".png")

#
# 振幅(比較用)
#
fig6 = plt.figure()
ax = fig6.add_subplot(111)
ax.set_xlabel('frequency / Hz')
ax.set_ylabel('Amplitude')

ax.annotate('', xy=[F_Nyq, 0], xytext=[F_Nyq, 0.05],
            arrowprops=dict(shrink=0, width=0.7, headwidth=4,
                            headlength=7, connectionstyle='arc3',
                            facecolor='gray', edgecolor='gray'))
ax.annotate('', xy=[-F_Nyq, 0], xytext=[-F_Nyq, 0.05],
            arrowprops=dict(shrink=0, width=0.7, headwidth=4,
                            headlength=7, connectionstyle='arc3',
                            facecolor='gray', edgecolor='gray'))
ax.text(F_Nyq, 0.06, "$F_N$", size=10, horizontalalignment="center")
ax.text(-F_Nyq, 0.06, "$-F_N$", size=10, horizontalalignment="center")

ax.scatter(freq_raw[0:int(N_raw)], amp_raw[0:int(N_raw)], s=3, label="Original", c="b")
ax.plot(freq_raw[0:int(N_raw)], amp_raw[0:int(N_raw)], c="b", linewidth=0.5)

if resample > 1:
    ax.scatter(freq[0:int(N)], amp[0:int(N)], s=3, label="Downsampled", c="g")
    ax.plot(freq[0:int(N)], amp[0:int(N)], c="g", linewidth=0.5)
    plt.legend(loc='best')

fig6.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95)
fig6.savefig(outdir+"/Amplitudes-"+name+".png")

#
# ダウンサンプリングしたデータを保存
#
csv = np.stack([raw_x, raw_y], 1)
np.savetxt(outdir + "/" + name + "-" + "downsample" + ".csv", csv, fmt='%.10f', delimiter=",")

#
# 結果の出力
#
print("Nyquist frequency:", F_Nyq)
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
