#Implementação do algoritmo do SAFT para ensaio FMC
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

#variáveis de análise
t = np.zeros((1858,1))
for linha in range(1858):
    t_i = (t_init + linha/f) #tempo das 1858 amostras em microsegundos
    t[linha] += t_i

z=(cl*t/2)*1e-3 #posições no eixo z em mm

#reduzindo a ROI para os 4 furos
#eixo x vai de -11 a 3 [mm]
#eixo z vai de 12 a 28 [mm]
x_min = np.argmin(np.abs(x-(-11)))
print('x_min = {}' .format(x_min))
x_max = np.argmin(np.abs(x-3))
print('x_max = {}' .format(x_max))
z_min = np.argmin(np.abs(z-12))
print('z_min = {}' .format(z_min))
z_max = np.argmin(np.abs(z-28))
print('z_max = {}' .format(z_max))

#conversões de tempo
t = np.zeros((869,1))
for linha in range(36,905):
    t_i = (t_init + linha/f) #tempo das amostras reduzidas em microsegundos
    t[linha-36] += t_i

z_reduz = (cl*t/2)*1e-3 #posições reduzidas no eixo z em mm

x_reduz = np.zeros((25,1)) #posições reduzidas no eixo x em mm
for coluna in range(13,38):
    x_reduz[coluna-13] = x[coluna]
print('x_reduz = {}' .format(x_reduz))

#ascans relativos a região de interesse
g = np.zeros((869,25))
for elemento in range(13,38):
    for eixo_z in range(36,905):
        g[eixo_z-36, elemento-13] = ascans[eixo_z, elemento, elemento]
#print('g = {}' .format(g))
#print(np.shape(g))

#SAFT para modelo FMC de ensaio 
def saft (g,x_reduz,z_reduz):
    f = np.zeros_like(g)
    for j in range(25):
        for i in range(869):
            coord_x = x_reduz[j]
            coord_z = z_reduz[i]
            for transdutor in range(25):
                coord_transd_x = x_reduz[transdutor]
                coord_transd_z = 0

                delta_x = coord_transd_x - coord_x
                delta_z = coord_transd_z - coord_z
                delta_x_pow2 = np.power(delta_x, 2)
                delta_z_pow2 = np.power(delta_z, 2)
                dist = np.sqrt(delta_x_pow2 + delta_z_pow2)

                print('j = {}'.format(j))
                print('i = {}'.format(i))
                print('transdutor = {}'.format(transdutor))

                indice = np.argmin(np.abs(z_reduz - dist))
                print('indice = {}'.format(indice))

                f[i, j] += g[indice,transdutor]

    print('final = {}'.format(f))
    return f

plt.figure()
f = saft(g, x, z_reduz)
plt.imshow(f, aspect="auto")
plt.title('SAFT')
plt.show()
