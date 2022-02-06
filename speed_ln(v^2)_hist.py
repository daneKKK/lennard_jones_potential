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


f=open('saves/'+input('Введите имя папки с данными, по которой нужно построить график энергий:\n') + '/_speed.txt','r')

m = 17

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
    #t2_x = [t ** 2 for t in t_x]
    lnf_x = []
    t2_x = []
    for i in f_x:
        if i == 0:
            t_x = np.delete(t_x, f_x.index(i))
            f_x.remove(i)
        else:
            t2_x += [t_x[f_x.index(i)] ** 2]
            lnf_x += [np.log(i)]
        #print(len(t_x))
        #print(len(lnf_x))
    #print(len(t_x))
    #print(len(f_x))
    #t2_x = [t ** 2 for t in t_x]
    #print(len(t2_x))
    #print(len(lnf_x))
    plt.plot(t2_x, lnf_x, 'o')
    p_x, v_x = np.polyfit(t2_x, lnf_x, deg=1, cov=True)
    p_f_x = np.poly1d(p_x)
    plt.plot([i ** 2 for i in t_x], p_f_x([i ** 2 for i in t_x]))
    plt.title('Распределение компонент скорости по оси ' + axis, fontsize=7)
    plt.xlabel('Квадрат скорости, (м/с)^2', fontsize=7)
    plt.ylabel(r'ln(f(v^2))', fontsize=7)
    return p_x[0]

sp = plt.subplot(131)
p_x = speedSubplot(vx, m, 'X')

sp = plt.subplot(132)
p_y =speedSubplot(vy, m, 'Y')

sp = plt.subplot(133)
p_z = speedSubplot(vz, m, 'Z')



m = 6.63 * (10 ** (-26))
k = 1.38 * (10 ** (-23))

T = (-m) / (k * 2 * p_x)

print(T, 'T_x')

T = (-m) / (k * 2 * p_y)

print(T, 'T_y')

T = (-m) / (k * 2 * p_z)

print(T, 'T_z')

plt.show()
