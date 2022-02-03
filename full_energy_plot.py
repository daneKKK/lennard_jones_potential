import numpy as np
import matplotlib.pyplot as plt

filename = input('Введите имя папки с данными, по которой нужно построить график энергий:\n')

f=open('saves/'+filename + '/_energy.txt','r')
t=[]
sumE=[]
sumKE=[]
sumU=[]
x0=0
y0=0
z0=0
k = 1
for line in f:
    x0,T,U=line.split()
    t=t+[float(x0)]
    sumE=sumE+[float(T)+k*float(U)]
    sumKE=sumKE+[float(T)]
    sumU=sumU+[k*float(U)]
f.close()



plt.plot(t,sumE, 'b.', label='Сумма энергий')
plt.plot(t,sumKE,'r.', label='Кинетическая энергия')
plt.plot(t,sumU,'y.', label='Потенциальная энергия')
print('Синий - сумма энергий \n' +
      'Красный - кинетическая энергия \n' +
      'Жёлтый - потенциальная энергия \n')



del sumE[0:1500]
    
print(abs(np.std(sumE)/np.mean(sumE)*100), ' % - относительное отклонение энергии')

del sumKE[0:1500]

KE = sum(sumKE)/len(sumKE)/300

print(KE, 'Дж - средняя кинетическая энергия')

print(KE * 2 / (3 * 1.38 * (10 ** (-23))), ' К - температура через среднюю кинетическую энергию')

print(np.std(sumKE)/np.mean(sumKE)*100, ' % - относительное отклонение кинетической энергии')

plt.title('Графики энергий. Название файла: ' + filename)
plt.ylabel('Энергия в джоулях')
plt.xlabel('Шаг симуляции')



plt.show()

#физика:
#t = 1.094 * 10^-12 с
#m0 = 6.63 * 10^-26 кг
#epsilon = 120 K * boltzmann_constant
#sigma = 0.34 * 10^-9 м
#v_числ = 0.003217 * v [м/с]
