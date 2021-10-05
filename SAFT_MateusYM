#Ritual de Iniciação - Implementação do SAFT
#Aluno: Mateus Yamada Muller
#Orientador: Thiago A. R. Passarin

import matplotlib.pyplot as plt
import numpy as np
import scipy.io
import math

file = "DadosEnsaio.mat"

mat = scipy.io.loadmat(file)

pi = 500  # Amostra inicial dos A-Scan
pf = 900  # Amostra final dos A-Scan
g = mat["ptAco40dB_1"]["AscanValues"][0][0][pi:pf]  # B-Scan
cl = mat["ptAco40dB_1"]["CscanData"][0][0]["Cl"][0][0][0][0]  # Velocidade
t = mat["ptAco40dB_1"]["timeScale"][0][0][pi:pf]*1e-6  # Tempo
T = t[1][0]-t[0][0]  # Período de amostragem
z = cl*t/2  # Conversação para posição /2->ida/volta
x = mat["ptAco40dB_1"]["CscanData"][0][0]["X"][0][0]*1e-3  # Posições transdut

def saft(g, x, z, cl, T, t):
    f = np.zeros_like(g)
    for j in range(31):
      for i in range(400):
        coord_x = j*1e-3  #coordenadas para eixo x da ROI e passo de 1mm = 1e-3m
        coord_z = ((t[i]*cl)/2) #[m]*[m/s] = [m] coordenadas para eixo z da ROI
        for transdutor in range(31):
          coord_transd_x = transdutor*1e-3 #coordenada para eixo x do transdutor e passo de 1mm
          coord_transd_z = 0 #transdutor se mantém em 0 para o eixo z  

          delta_x = coord_transd_x - coord_x #distancia no eixo x entre coordenada da ROI e o centro do transdutor
          delta_z = coord_transd_z - coord_z #distancia no eixo z entre coordenada da ROI e o centro do transdutor
          delta_x_pow2 = np.power(delta_x, 2)
          delta_z_pow2 = np.power(delta_z, 2)

          dist = np.sqrt(delta_x_pow2 + delta_z_pow2) #teorema de pitágoras dist = raiz_quadrada(delta_x^2 + delta_z^2)
          print('dist = {}' .format(dist))
  
          print('j = {}' .format(j))
          print('i = {}' .format(i))
          print('transdutor = {}' .format(transdutor))

          indice = np.argmin(np.abs(z-dist)) #índice relativo ao argumento mínimo entre a distancia convertida e a distancia calculada
          print('indice = {}' .format(indice))

          f[i,j] += g[indice,transdutor] #matriz da imagem composta pela amplitude do Ascan no índice e o valor anterior de amplitude já contido
          

    f = f/31
    print('final = {}' .format(f))
    return f

plt.figure()
plt.imshow(g, aspect="auto")
plt.title('B-Scan')

plt.figure()
f = saft(g, x, z, cl, T, t)
plt.imshow(f, aspect="auto")
plt.title('SAFT')
