import numpy as np
import matplotlib.pyplot as plt

def dn(x, interval):
    '''
    Возвращает количество элементов в интервале
    '''
    x1, x2 = interval
    dn = 0
    for i in x:
        if x1 <= i and i < x2:
            dn +=1
    return dn

def function(x, arr, interval):
    '''
    Гистограмма, являющаяся приближением функции плотности вероятности
    Функция возвращает один её столбик в интервале от x до x + dv
    '''
    dv = (max(arr) - min(arr))/interval
    
    return dn(arr, (x, x + dv)) / len(arr) / dv

def speedSubplot(vx, m):
    t_x = np.linspace(min(vx),max(vx),m)
    f_x = [function(i, vx, m) for i in t_x]
    plt.bar(t_x, f_x, width = t_x[1]-t_x[0])
    plt.title('Распределение скоростей', fontsize=7)
    plt.xlabel('Скорость, м/с', fontsize=7)
    plt.ylabel('Вероятность', fontsize=7)

filename = input('Введите имя папки с данными, по которой нужно построить распределение скоростей:\n')

f = open('saves/' + filename + '/_speedAbs.txt', 'r')

v = []
for line in f:
    v += [float(i) for i in line.split()]

n = len(v)
n = int(n**(1/2))

plt.subplot(121)
speedSubplot(v, n)

plt.subplot(122)

t = np.linspace(min(v), max(v), n)
f = [function(i, v, n) for i in t]

dv = (max(v) - min(v)) / n
v2 = [(i + 1/2 * dv)**2 for i in t]
lnfv = [np.log(f[i] / v2[i]) for i in range(len(f))]

del v2[0:40]
del v2[int(40/75 * n):n]
del lnfv[0:40]
del lnfv[int(40/75 * n):n]

plt.plot(v2, lnfv, 'o')

p, v = np.polyfit(v2, lnfv, deg=1, cov=True)
p_f = np.poly1d(p)
plt.plot(v2, p_f(v2))

m = 6.63 * (10 ** (-26))
k = 1.38 * (10 ** (-23))

T = (-m) / (2 * k * p[0])
print(T, ' - температура в кельвинах через наклон графика')

T2 = m / (2 * np.pi * k * (np.exp(p[1] * 2 / 3)))
print(T2, ' - температура через сдвиг графика')

plt.show()

