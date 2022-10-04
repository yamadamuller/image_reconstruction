import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.signal import convolve, gausspulse, hilbert
from framework import post_proc
import math

#rugo_data = np.load('D:/Documentos/AUSPEX/RUGO/c_matrix_f520_lc25_n6.xlsx.npy',allow_pickle=True)

fig, axi = plt.subplots()

Ts = 1/(160*1e6)
t_final = 12e-6
t = np.arange(0,t_final,Ts)
delay = 1e-6
fc = 3e6
bw = .6
s = gausspulse(t-delay, fc, bw)

dx = .1e-3
dz = dx

Nt = round(t_final / Ts)
#Nt = 2048+250
#Nx = rugo_data.shape[0]
#Nz = rugo_data.shape[1]
#Nt = 1600
Nz = 160
Nx = 160
#Nt = 1500
#Ne = 64
dtype = 'float32'  # float32 is the fastest

c_agua = (1740 * Ts / dx)
c_alum = (6000 * Ts / dx)
c = c_agua*np.ones((Nz,Nx), dtype=dtype)
#c = np.copy(rugo_data)
print('Velocidade adimensional: meio 1 = {} e meio 2 = {}'.format(c_agua,c_alum))

#rec_x = np.arange(0, Nx, 8)
rec_x = np.array([Nx//2])
rec_z = np.zeros_like(rec_x)
s_rec = np.zeros((len(rec_x), Nt))
idx_r = np.arange(len(rec_x))

range_x = np.arange(Nx)
range_z = np.arange(Nz)
xx,zz = np.meshgrid(range_x, range_z)

c[zz > 60] = c_alum #refletor plano

#xtransd = np.arange(0, Nx, 8)
#ztransd = np.linspace(0, 0, Ne)

# t = np.arange(Nt)
# sig = np.zeros((Nt,Nz,Nx))
# #source = np.zeros_like(t, dtype=dtype)
# #source[:len(s)] += s
# for i in range(len(rec_x)):
#     sig[:, rec_z[i],rec_x[i]] += s

#coeff = [-5269/1800, 5/3, -5/21, 5/126, -5 /
        # 1008, 1/3150]  # Suggested in Dablain1986

coeff = [-1077749 / 352800, 16 / 9, -14 / 45, 112 / 1485, -7 / 396, 112 / 32175, -2 / 3861, 16 / 315315,
             -1 / 411840]

def lap(u, y=np.zeros((Nz, Nx), dtype=dtype)):
    y[:] = 2 * coeff[0] * u
    for k, c in enumerate(coeff[1:], start=1):
        y[:-k, :] += c * u[k:, :]
        y[k:, :] += c * u[:-k, :]
        y[:, :-k] += c * u[:, k:]
        y[:, k:] += c * u[:, :-k]

    return y

def d2x(u, y=np.zeros((Nz, Nx), dtype=dtype)):
    y[:] = coeff[0]*u
    for k, c in enumerate(coeff[1:], start=1):
        y[:-k, :] += c*u[k:, :]
        y[k:, :] += c*u[:-k, :]

    return y


def d2z(u, y=np.zeros((Nz, Nx), dtype=dtype)):
    y[:] = coeff[0]*u
    for k, c in enumerate(coeff[1:], start=1):
        y[:, :-k] += c*u[:, k:]
        y[:, k:] += c*u[:, :-k]
    return y

# Laplacian
deriv_accuracy = 10
deriv_order = 2
deriv_n_coef = 2*np.floor((deriv_order+1)/2).astype('int')-1+deriv_accuracy
p = np.round((deriv_n_coef-1)/2).astype('int')
A = np.arange(-p, p+1)**np.arange(0, 2*p+1)[None].T
b = np.zeros(2*p+1)
b[deriv_order] = math.factorial(deriv_order)
hx = np.zeros((deriv_n_coef, deriv_n_coef))
hx[deriv_n_coef//2, :] = np.linalg.solve(A, b)
hy = hx.T
h = hx + hy

def sim_yield():
    u_2 = np.zeros((Nz, Nx), dtype=dtype)
    u_1 = np.zeros((Nz, Nx), dtype=dtype)
    u_0 = np.zeros((Nz, Nx), dtype=dtype)
    u_s = np.zeros((Nx, Nt), dtype=dtype)
    y = np.zeros((Nz, Nx), dtype=dtype)

    for k in range(2, Nt):
        u_1, u_2 = u_0, u_1
        #lap(u_1, y)
        #y = lap(u_1)
        ax = d2x(u_1)
        az = d2z(u_1)
        a11x = np.zeros((Nz, Nx), dtype=dtype)
        a11x[:, 1:-1] = ax[:, 2:] + ax[:, :-2] - 2*ax[:, 1:-1]
        a11z = np.zeros((Nz, Nx), dtype=dtype)
        a11z[1:-1, :] = az[2:, :] + az[:-2, :] - 2*az[1:-1, :]
        a11zx = np.zeros((Nz, Nx), dtype=dtype)
        a11zx[1:-1, :] = ax[2:, :] + ax[:-2, :] - 2*ax[1:-1, :]
        a11xz = np.zeros((Nz, Nx), dtype=dtype)
        a11xz[:, 1:-1] = az[:, 2:] + az[:, :-2] - 2*az[:, 1:-1]
        u_0 = 2*u_1 - u_2 + c**2*(ax + az) - c**4 * \
            (a11x + a11z + a11xz + a11zx)/12  # y
        u_0[Nx//Nx,Nz//2] += s[k]

        axi.cla()
        axi.imshow(u_0)
        axi.set_xlabel('Nx')
        axi.set_ylabel('Nz')
        axi.set_title("Shot {}".format(k))
        plt.pause(0.001)

        s_rec[idx_r,k] = u_0[rec_z[idx_r],rec_x[idx_r]]

        yield u_0

    return u_s

def sim_yield_CRAS():
    # u_2 = np.zeros((Nz, Nx), dtype=dtype)
    u_1 = np.zeros((Nz, Nx), dtype=dtype)
    u_0 = np.zeros((Nz, Nx), dtype=dtype)
    y = np.zeros((Nz, Nx), dtype=dtype)
    c2 = c**2
    for k in range(2, Nt):
        u_1, u_2 = u_0, u_1
        lap(u_1, y)
        u_0 = -u_2 + 2 * u_1 + y*c2
        u_0[Nx//Nx,Nz//2] += s[k]

        axi.cla()
        axi.imshow(u_0)
        axi.set_xlabel('Nx')
        axi.set_ylabel('Nz')
        axi.set_title("Shot {}".format(k))
        plt.pause(0.001)

        yield u_0


start_time = time.time()
u_s2 = np.array([u[0, :] for u in sim_yield()])
end_time = time.time()
print(f"{end_time-start_time} seconds")
# plt.imshow(u_s2, aspect='auto')

# Mostra o sinal recebido nos elementos do transdutor
plt.figure(2)
gate_start = 400
gate_end = 1400
ymax = s_rec.ravel().max()
t = np.arange(Nt)
for ir in range(len(s_rec)):
     idx_t = (t*Ts)-1e-6
     plt.plot(idx_t[gate_start:]*1e6, (s_rec[idx_r, gate_start:].T + ymax * idx_r))
     plt.xlabel('Time (us)')
     plt.ylabel('Amplitude')
plt.show()

plt.figure(3)
plt.imshow(c)
# plt.figure(1)
# plt.imshow(u[0,:,:].T, aspect='auto')
# plt.figure(2)
# plt.imshow(u_s.T, aspect='auto')
#plt.figure(3)
#plt.imshow(u_s2, aspect='auto')

#Testando posição do refletor
# idx_max_r = np.argmax(post_proc.envelope(s_rec[idx_r, gate_start:gate_end].T + ymax * idx_r))
# t_fl = idx_t[gate_start+idx_max_r]
# z_r = t_fl*(c_agua*ad)/2
# print('Posição do refletor estimada = {} mm'.format(z_r*1000))
# print('Posição do refletor na simulação = {} mm'.format((np.where(c==c_agua)[0][-1]*dx)*1000,2))

# def fft_zp(signal, fft_size=2 ** 11):
#     import scipy.fftpack as ft
#     temp = np.zeros((fft_size,1))
#     temp[0:len(signal)] = signal
#     temp = np.abs(ft.fft(temp))
#     #normalized_spc = temp / temp.max()
#     return temp
#
# fft_size = 2**11
# fs = 1/Ts
# freqs = np.arange(0, fs, fs/fft_size)
# sinal = s_rec[idx_r, gate_start:].T
# spc = fft_zp(sinal, fft_size)









