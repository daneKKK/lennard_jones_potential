import numpy as np
import matplotlib.pyplot as plt

filename = input('Введите имя папки с данными, по которой нужно построить график диффузии от времени:\n')

f = open('saves/' + filename + '/_diffusion.txt', 'r')

t = []
r = []

distUnit = 0
timeUnit = 0
unitsRecieved = False

for line in f:
    if not unitsRecieved:
        distUnit, timeUnit = line.split()
        distUnit = float(distUnit)
        timeUnit = float(timeUnit)
        unitsRecieved = True
        continue
    t0, r0 = line.split()
    if len(t) >= 1:
        t += [float(t0) * timeUnit - t[0]]
        r += [float(r0) * distUnit * distUnit]
    else:
        print(t0)
        t += [float(t0) * timeUnit]
        r += [float(r0) * distUnit * distUnit]

t[0] = 0

plt.plot(t, r, 'o', markersize = 0.5)

p, v = np.polyfit(t, r, deg=1, cov=True)
p_f = np.poly1d(p)
plt.plot(t, p_f(t))

print('D = ', p[0], 'м^2/c')

plt.title('Зависимость среднего квадрата отклонения от времени')
plt.xlabel('Время, с')
plt.ylabel('Средний квадрат отклонения, м')

plt.show()
