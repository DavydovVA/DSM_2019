from random import random
from math import log, exp
from matplotlib import pyplot as plt
import numpy as np
from scipy.stats.distributions import chi2
from time import time
import pprint

N = 500000 #количетсво потоков
a = 1
b = 2
w_A = [0, a]


source_array = []
temp = []


def intens(t, tau):
	return exp(-b * (t - tau))

	
def moment1(t):
	return (1 / 2) * (1 - exp( -2 * t))
	

def moment2(t1, t2):
	first = 0
	if t1 < t2:
		first = exp(-2 * t1) * exp( -2 * t2) * exp(4 * t1) / 4 - exp( -2 * t1) * exp( -2 * t2) / 4
	else:
		first = exp(-2 * t1) * exp( -2 * t2) * exp(4 * t2) / 4 - exp( -2 * t1) * exp( -2 * t2) / 4
		
	second = exp(-2 * t1) * exp( -2 * t2) * (exp(2 * t1) / 2 - 1 / 2) * (exp(2 * t2) / 2 - 1 / 2)
	
	return first + second


#################### Generator
for i in range(0, N):
	i_flow_B = []
	t_prev_A = w_A[0];
	while t_prev_A < w_A[1]:
		t_curr_A = t_prev_A - log(random())
		t_prev_B = t_curr_A
		flows_B = []
		while t_prev_B < w_A[1]:
			t_curr_B = t_prev_B - log(random())
			yi = random()
			if yi < intens(t_curr_B, t_curr_A) and t_curr_B < w_A[1]:
				flows_B.append(t_curr_B)
				temp.append(t_curr_B)
			
			t_prev_B = t_curr_B
		
		i_flow_B.append(flows_B)
		t_prev_A = t_curr_A
	
	source_array.append(i_flow_B)

sum_of_flows = []
for i in source_array:
	temp_flow = []
	for j in i:
		for k in j:
			temp_flow.append(k)
	sum_of_flows.append(temp_flow)
####################################################
#Интенсивность
hist, bins = np.histogram(temp, 25, (0, w_A[1]))
s = bins[1] - bins[0]
Iexp = hist / (N * s)

centrs = list()
for i in range(0, len(bins) - 1):
	centrs.append(bins[i] + s / 2)

Iteor = list()
for i in range(0, len(bins) - 1):
	Iteor.append(moment1(centrs[i]))
	
fig = plt.figure()
plt.plot(centrs, Iteor, c='red')
plt.plot(centrs, Iexp)
plt.title('Интенсивность')
#################################################### 
#корреляционная
K = np.zeros((N, len(centrs)))
ind = 17
for i in range(0, len(sum_of_flows)):
	hist, bins = np.histogram(sum_of_flows[i], 25, (0, w_A[1]))
	
	K[i] = hist*hist[ind] / (s * s)
	K[i, ind] = hist[ind] * (hist[ind] - 1) / (s * s)

summ = np.zeros((25))
for i in range(0, len(sum_of_flows)):
	summ = K[i] + summ

korExp = summ/(N)

cmax = bins[17]
korTeor = []
for i in range(0, len(bins) - 1):
	korTeor.append(moment2(bins[i], cmax))

fig = plt.figure()
plt.plot(bins[0:25], korTeor, c='red')
plt.plot(bins[0:25], korExp)
plt.title('Корелл-ая ф-ция')
plt.show()
