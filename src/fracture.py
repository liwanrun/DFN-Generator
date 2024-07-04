import math
import numpy as np
import matplotlib.pyplot as plt

class Fracture:
    def __init__(self, npts) -> None:
        self.numPts = npts
        self.Points = np.zeros(shape=(npts, 2), dtype=float)

    def translate(self, target):
        self.Points += target

    def rotate(self, degree):
        radian = math.radians(degree)
        matrix = np.array([[np.cos(radian), np.sin(radian)],
                           [-np.sin(radian), np.cos(radian)]])
        self.Points = np.dot(self.Points, matrix) 

    def scale(self, factor):
        self.Points *= factor

    def write_to_gmsh(self):
        pass


class SmoothFracture(Fracture):
    def __init__(self, npts=2) -> None:
        super().__init__(npts)

    def write_to_gmsh(self):
        pass


class RoughFracture(Fracture):
    '''
    Generate rough fractures using FFT transform
    '''
    def __init__(self, npts=100, d238=(0.0218, 0.0111, 0.0019)) -> None:
        super().__init__(npts)
        self.Points[:, 0] = np.linspace(-0.5, 0.5, npts) 
        amp, phi = self._calc_harmonic(d238)
        self._transform_IFFT(amp, phi)

    def _calc_harmonic(self, d238):
        num = self.numPts//2 + 1
        phi = np.random.uniform(-np.pi, np.pi, num)
        amp = np.zeros(num, dtype=float)
        amp[2], amp[3], amp[8] = d238

        alpha = -1.4042
        beta = -1.3489
        if amp[3] > 0.0:
            amp[4:8] = 2**(alpha*np.log2(np.arange(4, 8)/3) + np.log2(amp[3]))
        if amp[8] > 0.0:
            amp[9:] = 2**(beta*np.log2(np.arange(9, num)/8) + np.log2(amp[8]))
        return amp, phi 

    def transform_DFT(self, x, y):
        N = len(y)
        self.Points[:, 0] = x
        self.Points[:, 1] = y
        a_n = np.zeros(N//2 + 1, dtype='float')
        b_n = np.zeros(N//2 + 1, dtype='float')
        for k in range(N // 2 + 1):
            for n in range(N):
                a_n[k] += y[n] * np.cos(2 * np.pi * k * n / N)
                b_n[k] += y[n] * np.sin(2 * np.pi * k * n / N)
        a_n[0] /= (N)
        a_n[1:] /= (N/2)
        b_n[1:] /= (N/2)
        return a_n, b_n
    
    def transform_IDFT(self, d, p):
        '''
        [d]: input magnitudes of harmonics
        [p]: input phase angle of harmonics
        '''
        N = self.numPts  
        y = self.Points[:, 1]
        for n in range(N):
            y[n] = d[0]
            for k in range(1, N//2 + 1):
                y[n] += d[k] * np.cos(2 * np.pi * k * n / N + p[k])
                # y[n] += d_n[k] * np.cos(2 * np.pi * k * n / N) + p_n[k] * np.sin(2 * np.pi * k * n / N)
    
    def transform_FFT(self, x, y):
        N = len(y)
        self.Points[:, 0] = x
        self.Points[:, 1] = y
        dft = np.fft.rfft(y)
        a_n = dft.real * 2 / N
        b_n = dft.imag * 2 / N
        a_n[0] = dft[0].real / N
        if N % 2 == 0:
            a_n[-1] = dft[-1].real / N
            b_n[-1] = 0.0
        return a_n, b_n

    def _transform_IFFT(self, d, p):
        num = self.numPts  
        a_n = d * np.cos(p)
        b_n = d * np.sin(p)
        dft = (a_n - 1j * (b_n)) * num / 2
        dft[0] = a_n[0] * num + 1j * 0
        if num % 2 == 0:
            dft[-1] = a_n[-1] * num + 1j * 0
        self.Points[:, 1] = np.fft.irfft(dft, num)

    def write_to_gmsh(self, file):
        for i, point in enumerate(self.Points):
            file.write(f'Point({i}) = {{{point[0]}, {point[1]}, {point[2]}, 1.0}};')
        

# Unit test
if __name__ == '__main__':
    # x = np.linspace(0.0, 1.0, num=100)
    # freq = 1.
    # y = 3*np.sin(2*np.pi*freq*x)

    # freq = 4
    # y += np.sin(2*np.pi*freq*x)

    # freq = 7   
    # y += 0.5* np.sin(2*np.pi*freq*x)
    
    # # 计算DFT
    # frac = RoughFracture()
    # A, B = frac.transform_FFT(x, y)  
    # # # 恢复原信号
    # A, B = np.sqrt(A**2 + B**2), np.arctan2(-B, A)
    # # frac.transform_IFFT(A, B, len(y))
    # frac.transform_IFFT(A, B)

    fractures = []
    rf_1 = RoughFracture(npts=100, d238=(np.random.uniform(0.0000, 0.0218), 
                                         np.random.uniform(0.0000, 0.0111),
                                         np.random.uniform(0.0000, 0.0019)))
    fractures.append(rf_1)

    rf_2 = RoughFracture(npts=100, d238=(np.random.uniform(0.0220, 0.0300), 
                                         np.random.uniform(0.0113, 0.0166),
                                         np.random.uniform(0.0020, 0.0028)))
    fractures.append(rf_2)

    rf_3 = RoughFracture(npts=100, d238=(np.random.uniform(0.0321, 0.0360), 
                                         np.random.uniform(0.0170, 0.0195),
                                         np.random.uniform(0.0029, 0.0039)))
    fractures.append(rf_3)

    rf_4 = RoughFracture(npts=100, d238=(np.random.uniform(0.0370, 0.0395), 
                                         np.random.uniform(0.0200, 0.0218),
                                         np.random.uniform(0.0040, 0.0050)))
    fractures.append(rf_4)

    rf_5 = RoughFracture(npts=100, d238=(np.random.uniform(0.0400, 0.0438), 
                                         np.random.uniform(0.0222, 0.0247),
                                         np.random.uniform(0.0051, 0.0061)))
    fractures.append(rf_5)

    rf_6 = RoughFracture(npts=100, d238=(np.random.uniform(0.0440, 0.0475), 
                                         np.random.uniform(0.0248, 0.0284),
                                         np.random.uniform(0.0062, 0.0072)))
    fractures.append(rf_6)

    rf_7 = RoughFracture(npts=100, d238=(np.random.uniform(0.0481, 0.0500), 
                                         np.random.uniform(0.0288, 0.0311),
                                         np.random.uniform(0.0073, 0.0085)))
    fractures.append(rf_7)

    rf_8 = RoughFracture(npts=100, d238=(np.random.uniform(0.0520, 0.0561), 
                                         np.random.uniform(0.0322, 0.0345),
                                         np.random.uniform(0.0086, 0.0100)))
    fractures.append(rf_8)

    rf_9 = RoughFracture(npts=100, d238=(np.random.uniform(0.0600, 0.0755), 
                                         np.random.uniform(0.0380, 0.0440),
                                         np.random.uniform(0.0100, 0.0110)))
    fractures.append(rf_9)

    rf_10 = RoughFracture(npts=100, d238=(np.random.uniform(0.0782, 0.0800), 
                                          np.random.uniform(0.0473, 0.0538),
                                          np.random.uniform(0.0110, 0.0130)))
    fractures.append(rf_10)

    fig, axs = plt.subplots(nrows=10, ncols=1, figsize=(6, 8), sharex=True,layout='tight')
    for i in range(10):
        fractures[i].translate((0.5, 0))
        fractures[i].scale(100.0)
        x = fractures[i].Points[:, 0]
        y = fractures[i].Points[:, 1]
        axs[i].plot(x, y, 'k', label=f'JRC = {2*i}-{2*(i+1)}')
        axs[i].set_xlim([0.0, 150.0])
        axs[i].set_ylim([-50, 50.0])
        axs[i].legend()
    plt.savefig('./img/numerical_JRC.png', dpi=300)
    plt.show()