#Método TFM para reconstrução de imagens (exclusivo para ensaio Full Matrix Capture) com loops de for
#Aluno: Mateus Yamada Muller
#Orientador: Thiago A.R. Passarin

import matplotlib.pyplot as plt
import numpy as np

file = 'D:\Documentos\AUSPEX\IMPLEMENTAÇÕES\TFM\TFM\dados_tfm.npy'

dados = np.load(file, allow_pickle=True).item()
ascans = dados.get('ascans')
speed_m_s = dados.get('speed_m_s')
f_sampling_MHz = dados.get('f_sampling_MHz')
samples_t_init_microsec = dados.get('samples_t_init_microsec')
elem_positions_mm = dados.get('elem_positions_mm')

#renomeando as variáveis
cl = speed_m_s #velocidade
f = f_sampling_MHz #frequência
t_init = samples_t_init_microsec #tempo
x = elem_positions_mm #posições transdutor

#conversões matemáticas
t = np.zeros((1858,1))
for linha in range(1858):
    t_i = (t_init + linha/f) #tempo das 1858 amostras em us
    t[linha] += t_i

z=(cl*t)*1e-3 #posições no eixo z em mm

#reduzindo a ROI para os 4 furos
#eixo x vai de -11 a 3 [mm]
#eixo z vai de 12 a 28 [mm]
x_min = np.argmin(np.abs(x-(-10)))
#print('x_min = {}' .format(x_min))
x_max = np.argmin(np.abs(x-(-5)))
#print('x_max = {}' .format(x_max))
z_min = np.argmin(np.abs(z-12))
#print('z_min = {}' .format(z_min))
z_max = np.argmin(np.abs(z-17.5))
#print('z_max = {}' .format(z_max))

#ascans necessários para algoritmo
g = ascans[:, :, :]
#print('g = {}' .format(g))
#print(np.shape(g))

#TFM para modelo FMC de ensaio
def tfm (g,x,z,t,cl):
    s_sum = 0
    f = np.zeros((869,25))
    for j in range(13,38):
        for i in range(36,905):
            coord_x = x[j]  #coordenadas para para o eixo x da ROI
            coord_z = z[i]/2    #coordenadas para o eixo z da ROI
            for emissores in range(64):
                for receptores in range(64):
                    coord_emissor_x = x[emissores]  #coordenadas eixo x do emissor
                    coord_emissor_z = 0
                    coord_receptor_x = x[receptores] #coordenadas eixo x do receptor
                    coord_receptor_z = 0

                    delta_e_x = coord_emissor_x - coord_x
                    delta_e_z = coord_emissor_z - coord_z
                    delta_r_x = coord_receptor_x - coord_x
                    delta_r_z = coord_receptor_z - coord_z

                    delta_e_x_pow2 = np.power(delta_e_x, 2)
                    delta_e_z_pow2 = np.power(delta_e_z, 2)
                    delta_r_x_pow2 = np.power(delta_r_x, 2)
                    delta_r_z_pow2 = np.power(delta_r_z, 2)

                    dist1 = np.sqrt(delta_e_x_pow2 + delta_e_z_pow2)
                    dist2 = np.sqrt(delta_r_x_pow2 + delta_r_z_pow2)
                    dist = dist1 + dist2
                    t_e_r = (dist/cl)*1e3
                    #print(dist)

                    #print('j = {}'.format(j))
                    #print('i = {}'.format(i))
                    #print('emissor = {}'.format(emissores))
                    #print('receptor = {}'.format(receptores))

                    indice1 = np.argmin(np.abs(z - dist))
                    indice2 = np.argmin(np.abs(t - t_e_r))
                    #print('indice1 = {}'.format(indice1))
                    #print('indice2 = {}'.format(indice2))

                    f[i-36, j-13] += g[indice2, emissores, receptores]
                    s_sum += 1 #registra, ao final, o número total de operações de soma envolvidas

    print('sinais somados = {}' .format(s_sum))
    f = np.abs(f) #retorna o valor absoluto (módulo) dos valores de f
    print('final = {}'.format(f))
    print(np.shape(f))
    return f


plt.figure()
f = tfm(g, x, z, t, cl)
plt.imshow(f, aspect="auto")
plt.title('TFM')
plt.show()
