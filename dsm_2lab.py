from random import random
from math import log, exp
from matplotlib import pyplot as plt
import numpy as np
from scipy.stats.distributions import chi2

N = 10000 #количетсво потоков
W = [0, 5] #omega
source_array = []
temp = [] #каждый элемент - поток с k событиями

def intens(t):
	return exp(-t)
	
#################### Generator
for i in range(0, N):
	i_flow = []
	t_prev = 0;
	while t_prev < W[1]:
		t_curr = t_prev - log(random())
		yi = random() # [0, 1]
		if yi < intens(t_curr) and t_curr < W[1]:
			i_flow.append(t_curr)
			temp.append(t_curr)
		t_prev = t_curr
	source_array.append(i_flow)
######################################################
max_elem_array = []
min_elem_array = []
len_array = []
for i in source_array:
	len_array.append(len(i))
	if i:
		max_elem_array.append(max(i))
		min_elem_array.append(min(i))
'''
print(f'max len: {max(len_array)}')
print(f'min len: {min(len_array)}')
print(f'min elem: {min(min_elem_array)}')
print(f'max elem: {max(max_elem_array)}')
'''
#####################################################
#Интенсивность
hist, bins = np.histogram(temp, 30, (0, 5))
s = bins[1] - bins[0]
Iexp = hist / (N * s)
print(Iexp/N)
print(hist)
centrs = list()
for i in range(0, len(bins) - 1):
	centrs.append(bins[i] + s / 2)

Iteor = list()
for i in range(0, len(bins) - 1):
	Iteor.append(intens(centrs[i]))
	
fig = plt.figure()
plt.plot(Iteor, centrs, c='red')
plt.plot(Iexp, centrs)
plt.title('Интенсивность')

####################################################
#корреляционная
enum_Iteor = enumerate(Iteor)
ind = 0
maxI = 0
for i, value in enum_Iteor:
	if value > maxI:
		maxI = value
		ind = i

K = np.zeros((N, len(centrs)))

for i in range(0, len(source_array)):
	hist, bins = np.histogram(source_array[i], 30, (0, 5))
	s = bins[1] - bins[0]
	
	K[i] = hist*hist[ind] / (s * s)
	K[i, ind] = hist[ind] * (hist[ind] - 1) / (s * s)

cmax = centrs[ind]
summ = np.zeros((30))
for i in range(0, len(source_array)):
	summ = K[i] + summ

korExp = summ/N

korTeor = []
for i in range(0, len(bins) - 1):
	korTeor.append(intens(centrs[i]) * intens(cmax))


fig = plt.figure()
plt.plot(bins[0:30], korTeor, c='red')
plt.plot(bins[0:30], korExp)
plt.title('Корелл-ая ф-ция')

###################################################
#веротность отдельного числа событий
max_length = max(len_array)

P = []

for i in range(0, max_length):
	P.append(((1**(i)) / np.math.factorial(i)) * exp(-1))
	
f = []
for i in range(0, max_length + 2):
	f.append(i)

Nk, bins = np.histogram(len_array, max_length, (0, max_length))
Rk = Nk / N

fig = plt.figure()
plt.plot(f[0:max_length], P)
plt.plot(f[0:max_length], Rk)
plt.title('Вероятность n-го числа событий')
plt.show()

###################################################]
#chi2
NP, mu = [], []

trigger, counter = 0, 0
mmm = 0

for i in range(0, max_length):
	if N * P[i] >= 5:
		NP.append(N * P[i])
		mu.append(Nk[i])
		counter += 1
		trigger = 0
		mmm += 1
	else:
		if trigger == 0:
			NP.append(N * P[i])
			mu.append(Nk[i])
			counter += 1
		else:
			NP[-1] = NP[-1] + (N * P[i])
			mu[-1] = mu[-1] + (Nk[i])
			counter += 1

#чтобы NP было равно N
mm = sum(NP[0:mmm])
mm_end = N - mm
NP[-1] += mm_end

mu = np.array(mu)
NP = np.array(NP)
mem = (mu - NP)**2/NP

chi2EXP = sum(mem)
print(chi2EXP)

chi2Pred = chi2.ppf(0.95, max_length - 1)
print(chi2Pred)

