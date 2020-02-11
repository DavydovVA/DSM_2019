from math import exp, log, sqrt
from random import random, gauss
import numpy as np
from matplotlib import pyplot as plt


W = [0, 20]
N = 10000

f1 = 1
sigma = 1 / sqrt(2 * f1)
alpha = 1 / 2
lambd = 1


def boze_einstein(m, n):
	return (m ** n) / ((m + 1) ** (n + 1))
	
	
def expression_p82(prev_1, prev_2):
	first = prev_1 ** 2
	second = prev_2 ** 2
	
	return first + second


def f2(t):
	global f1
	
	return f1 * f1 * (1 + exp(-t / 2))

	
def generator():
	source_array = [] #flows
	
	for i in range(0, N):
		i_flow = []
		t_prev = 0;
		while t_prev < W[1]:
			t_curr = t_prev - log(random()) / 9 # 9 - batta, max ksi
			yi = random() * 9
			if t_curr < W[1]:
				i_flow.append(t_curr)

			t_prev = t_curr
			
		source_array.append(i_flow)
		
	return source_array
		
		
		
if __name__ == '__main__':
	####################################### step1
	print('Step 1...')
	source_array = generator()
	
	max_elem_array = []
	min_elem_array = []
	len_array = []
	for i in source_array:
		len_array.append(len(i))
		if i:
			max_elem_array.append(max(i))
			min_elem_array.append(min(i))
	
	print(f'\nmax len: {max(len_array)}')
	print(f'min len: {min(len_array)}')
	print(f'min elem: {min(min_elem_array)}')
	print(f'max elem: {max(max_elem_array)}\n')
	####################################### step2
	print('Step 2...')
	for i, j in enumerate(source_array):
		source_array[i] = sorted(j)

	step_2 = []
	for sub in source_array:
		step_2_arr = []
		for i in range(len(sub)):
			if i == 0:
				curr_1 = sigma * gauss(0, 1)
				curr_2 = sigma * gauss(0, 1)
				step_2_arr.append(expression_p82(curr_1, curr_2))
				prev_1 = curr_1
				prev_2 = curr_2
			else:
				curr_1 = prev_1 * exp(-alpha * (sub[i] - sub[i - 1]) / 2) + sigma * gauss(0, 1) * sqrt(1 - exp(-2 * alpha * (sub[i] - sub[i - 1]) / 2))
				curr_2 = prev_2 * exp(-alpha * (sub[i] - sub[i - 1]) / 2) + sigma * gauss(0, 1) * sqrt(1 - exp(-2 * alpha * (sub[i] - sub[i - 1]) / 2))
				
				step_2_arr.append(expression_p82(prev_1, prev_2))	
				
				prev_1 = curr_1
				prev_2 = curr_2
				
		step_2.append(step_2_arr)
	####################################### step3
	print('Step 3...')
	step_3 = []
	
	for i, j in enumerate(step_2):
		step_3_arr = []
		for k, z in enumerate(j):
			rndm = 9 * random()
			if rndm < z:
				step_3_arr.append(source_array[i][k])
		step_3.append(step_3_arr)
	
	step_3_arr = []
	for i in step_3:
		for j in i:
			step_3_arr.append(j)
	#######################################
	print('Graph. 1...')
	#Интенсивность
	hist, bins = np.histogram(step_3_arr, 25, (0, W[1]))
	s = bins[1] - bins[0]
	Iexp = hist / (N * s)
	
	centrs = list()
	for i in range(0, len(bins) - 1):
		centrs.append(bins[i] + s / 2)

	Iteor = list()
	for i in range(0, len(bins) - 1):
		Iteor.append(f1)
		
	fig = plt.figure()
	plt.plot(centrs, Iteor, c='red')
	plt.plot(centrs, Iexp)
	plt.title('Интенсивность')
	####################################################
	print('Graph. 2...')
	#корреляционная
	
	K = []
	for i in range(len(step_3)):
		hist, bins = np.histogram(step_3[i], 25, (0, W[1]))
		
		K_temp = []
		for j in range(len(hist)):
			if j == 0:
				b = [val - 1 for val in hist]
				c = [hist[k] * b[k] / (s * s) for k in range(len(hist))]
				K_temp.append(sum(c) / (len(hist)))
			else:
				a = hist[j:]
				b = hist[:len(hist) - j]
				c = [a[k] * b[k] / (s * s) for k in range(len(a))]
				K_temp.append(sum(c) / (len(hist) - j))
		K.append(K_temp)
	
	summ = np.zeros((25))
	for i in range(len(K)):
		summ = K[i] + summ
	
	korExp = summ/(N)

	korTeor = []
	for i in range(len(bins) - 1):
		korTeor.append(f2(bins[i]))

	fig = plt.figure()
	plt.plot(bins[0:25], korTeor, c='red')
	plt.plot(bins[0:25], korExp)
	plt.title('Корелл-ая ф-ция')
	##################################################
	print('Graph. 3...')
	#распределение вероятностей
	tauc = 2
	m = f1 * tauc / 10
	
	tauc_step3 = []
	for i in step_3:
		ns3 = []
		for j in i:
			if j <= m:
				ns3.append(j)
		tauc_step3.append(ns3)
	
	lengths = []
	for i in tauc_step3:
		lengths.append(len(i))
		
	Nk, bins = np.histogram(lengths, max(lengths), (0, max(lengths)))
	Rk = Nk / len(tauc_step3)
	K = len(bins) - 1
	
	Teor = []
	for i in range(K):
		Teor.append(boze_einstein(m, i))
		
	fig = plt.figure()
	plt.plot(bins[0:K], Rk)
	plt.plot(bins[0:K], Teor)
	plt.title('Распределение числа событий на инт-ле [0, tuac/10]')
	plt.show()
