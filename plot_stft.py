import numpy as np
import matplotlib.pyplot as plt

# Parameters
sample_rate = 8000.0
duration = 1.0
signal_len = int(sample_rate * duration)

hop = 128
win = 256
fftsize = 256
outsize = (signal_len - win) // hop + 1

# Load signal
time_data = np.loadtxt("signal.txt", comments="#")

# Load STFT binary data
data = np.fromfile("stft_out.bin", dtype=np.float32)
data = data.reshape((outsize, fftsize // 2, 2))
complex_data = data[..., 0] + 1j * data[..., 1]

# Magnitude spectrogram
mag = np.abs(complex_data)

# Plot time signal
plt.figure(figsize=(12, 4))
plt.plot(time_data[:, 0], time_data[:, 1])
plt.title("Input Signal")
plt.xlabel("Time (s)")
plt.ylabel("Amplitude")
plt.grid(True)

# Plot spectrogram
plt.figure(figsize=(12, 6))
plt.imshow(20 * np.log10(mag.T + 1e-6), origin='lower', aspect='auto',
           extent=[0, duration, 0, sample_rate / 2], cmap="magma")
plt.colorbar(label="Magnitude (dB)")
plt.title("STFT Magnitude Spectrogram")
plt.xlabel("Time (s)")
plt.ylabel("Frequency (Hz)")
plt.tight_layout()

plt.savefig("stft-plot.png")
