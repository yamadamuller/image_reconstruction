#implementação inicial do método Newthon/Raphson

import matplotlib.pyplot as plt
import numpy as np

#derivada numérica - diferenças finitas
def derivada(f,x,delta):
    return (f(x+delta) - f(x))/delta

def f(x):
    return (-5.5*x**2 + 4*x + 10.5)

def min_max_funcao():
    lista = []
    x = np.linspace(0,1,100)
    delta = x[1] - x[0]
    for x_est in x:
        f_linha = derivada(f,x_est,delta)
        lista.append(f_linha)

    absolutos = np.abs(lista) #lista de todos os valores que o loop retorna
    minimo = min(absolutos)
    indice = np.argmin(np.abs(minimo-absolutos)) #índice do valor absoluto mínimo
    x_est_aprox = indice*delta #x_est ligado ao índice do mínimo
    return round(x_est_aprox,2)

print('O valor x de mínimo/máximo da função é aproximadamente = {}' .format(min_max_funcao()))

def newton_raphson():
    lista = []
    x = np.linspace(1.5,2,100)
    delta = x[1] - x[0]
    for x_est in x:
        f_linha = derivada(f, x_est, delta)
        iteracao = (f(x_est)/f_linha)
        iteracao_abs = np.abs(iteracao)
        if iteracao_abs <= 0.001:
            break
            
    return x_est

print('Valor aproximado de x que zera a função = {}'.format(np.round(newton_raphson(),3)))

#plot dos gráficos
x_graf =np.linspace(-2,2,100)
y=f(x_graf)
f_linha = derivada(f,x_graf,0.02)
plt.plot(x_graf,y,label='f(x)')
plt.plot(x_graf,f_linha,label="f'(x)")
plt.legend()
plt.grid(True)
plt.show()
