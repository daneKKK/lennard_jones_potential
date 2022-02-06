import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.stats import maxwell

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

def MSD(t):
    mean = np.mean(t)
    for i in t:
        sum_square = (i - mean) ** 2
    return np.sqrt(sum_square / len(t))

def plot(x, arr, interval):
    dv = (max(arr) - min(arr))/interval
    

f=open('saves/'+input('Введите имя папки с данными, по которой нужно построить график энергий:\n') + '/_speed.txt','r')

m = 17

#t - набор полученных скоростей

vx = []
vy = []
vz = []


for line in f:
    vx_0, vy_0, vz_0 = line.split()
    vx += [float(vx_0)]
    vy += [float(vy_0)]
    vz += [float(vz_0)]

def speedSubplot(vx, m, axis):
    t_x = np.linspace(min(vx),max(vx),m)
    f_x = [function(i, vx, m) for i in t_x]
    plt.bar(t_x, f_x, width = t_x[1]-t_x[0])
    plt.title('Распределение компонент скорости по оси ' + axis, fontsize=7)
    plt.xlabel('Скорость, м/с', fontsize=7)
    plt.ylabel('Вероятность', fontsize=7)

sp = plt.subplot(131)
speedSubplot(vx, m, 'X')

sp = plt.subplot(132)
speedSubplot(vy, m, 'Y')

plt.subplot(133)
speedSubplot(vz, m, 'Z')

sigma = np.mean([np.std(vx), np.std(vy), np.std(vz)])

T = sigma ** 2 * 6.63 * (10 ** (-26)) / (1.38 * (10 ** (-23)))

print(T, 'К - температура через отклонение компонент скорости от среднего')

plt.show()
