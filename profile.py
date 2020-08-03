"""
t_n = t_(n-1) + v_(n-1)
v_n = v_(n-1) + a_(n-1)
a_n = a_(n-1) + j

t_n = v_0 + v_1 + ... + v_n-1
v_n = a_0 + a_1 + ... + a_n-1

Velocity profile
T1 quadratic / T2 linear / T3  quadratic / T4 flat /T5  quadratic / T6 linear / T7  quadratic 

given
	j
	a
	v
	d

goal
	find t_1, t_2, t_3, t_4, t_5, t_6, t_7	

concatinate different motions together and form a smooth motion profile

motion_0 (j_0, a_0, v_0, d_0): start ramp up and achieve the max speed, and exit this motion with max speed,
d_0 and zero accelaration

motion_1 (j_1, a_1, v_1, d_1): start ramp up and achieve the max speed, and exit this motion with max speed,
d_1 and zero accelaration 

...
last motion
motion_e (j_e, a_e, v_e, d_e): start ramp up and achieve the max speed, and exit this motion with max speed,
d_1 and zero accelaration 

type_2: s-curve with no ending tail
	t_1, t_2, t_3, t_4

	a_m = j * t_1
	v_m = a_m (t_1 + t_2) + v_0
	d_m = v_m * (t_1 + t_2 /2 + t_4) + v_0 * (t_2/2 + t_1)
	j_m = j * (d/d_m)

type_3: 
	t_1, t_2, t_3, t_4, t_5, t_6, t_7
	start with v_s
	get to v_m during t_1, t_2, t_3
	stay in v_m during t_4
	get to zero during t_5, t_6, t_7 
"""

import math
import matplotlib.pyplot as plt
import time
import numpy as np

def tick(jerk, ticks, v_0 = 0, a_0 = 0):
	q = [0]
	v = [v_0]
	a = [a_0]
	
	for i in range(len(jerk)):
		for j in range(ticks[i]):
			q.append(q[-1] + v[-1])
			v.append(v[-1] + a[-1])
			a.append(a[-1] + jerk[i])
	return [a[1:], v[1:], q[1:]]

"""
remove intervals with no ticks
"""
def tick_jerk(t_j):
	data = []
	for tj in t_j:
		if tj[0] > 0:
			data.append({"tick": tj[0], "jerk": tj[1]})
	return data

"""
add segment order to each tick_jerk
0: first interval
1: middle interval
2: last interval
3: first and last interval
"""
def tick_jerk_segment_order(t_j):
	l = length(t_j)
	for i in range(length(t_j)):
		t_j[i]["segment_order"] = 1 
		if i == l-1 and i == 0:
			t_j[i]["segment_order"] = 3
		elif i == 0:
			t_j[i]["segment_order"] = 0
		elif i == l-1:
			t_j[i]["segment_order"] = 2 

	return t_j

"""
type_1: full s-curve
	t_1, t_2, t_3, t_4, t_5, t_6, t_7
	 
	a_m = j * t_1 < a
	v_m = a_m * (t_1 + t_2) + v_0= j * t_1 * (t_1 + t_2) + v_0 < v   
	t_m = 4*t_1 + 2*t_2 + t_4 >= t
	d_m = (v_m - v_0)*(2*t_1 + t_2 + t_4) + v_0 * (4*t_1 + 2*t_2 + t_4) = 
	j * t_1 * (t_1 + t_2) * (2*t_1 + t_2 + t_4) + v_0 * (4*t_1 + 2*t_2 + t_4) < d
	
input: j, a, v, d, v_0 
return: j_m, a_m, v_m, t_1, t_2, t_4
"""
def profile_1(j, a, v, d, v_0 = 0):
	
	t_1_a = math.floor(a / j)
	t_1_v = math.floor((abs(v-v_0)/j)**(1/2))
	# t_1_d
	roots = [ np.real(x) for x in np.roots([ np.sign(v-v_0) *2*j, 0, 4*v_0, -d]) if np.isreal(x)]
	roots.sort()
	if v >= v_0:
		# positive j
		t_1_d = math.floor(max(roots))
	else:
		# negative j
		if len(roots) <= 1:
			t_1_d = math.floor(((2*v_0)/(3*j))**(1/2))
		else:
			t_1_d = math.floor(roots[1]) 
	# find t_1_a
	t_1 = min(t_1_a, t_1_v, t_1_d)
	t_1 = max(0, t_1)
	t_2 = 0
	
	if t_1 == 0:
		if v_0 > 0: 
			t_4 = math.ceil(d/v_0)
			t_4 = max(0, t_4)
		else: 
			t_4 = 0
		j_m = 0

	else:
		if t_1_a <= min(t_1_v, t_1_d):

			t_2_v = math.floor((abs(v-v_0) / (j*t_1)) - t_1)
			#t_2_d
			roots = [ np.real(x) for x in np.roots([np.sign(v-v_0) *j*t_1, np.sign(v-v_0)*3*j*(t_1**2)+2*v_0, np.sign(v-v_0)*2*j*(t_1**3)+4*v_0*t_1 -d]) if np.isreal(x)]
			roots.sort()
			if v >= v_0:
				# positive j
				t_2_d = math.floor(max(roots))  
			else:
				# negative j
				if len(roots) <= 1:
					t_2_d = math.floor(((-3*j*(t_1**2) + 2*v_0)/(2*j*t_1))**(1/2))
				else:
					t_2_d = math.floor(roots[0])
			t_2 = min(t_2_v, t_2_d)
			t_2 = max(t_2, 0)
	
		t_4 = math.ceil((d - (2*t_1 + t_2)*(2*v_0 + np.sign(v-v_0) * j*t_1*(t_1+t_2)))/(np.sign(v-v_0) *j*t_1*(t_1+t_2)+v_0)) 
		t_4 = max(t_4, 0)
		j_m = (d- v_0 * (4*t_1 + 2*t_2 + t_4))/(t_1 * (t_1 + t_2) * (2*t_1 + t_2 + t_4))

	# a_m, v_m, j_m
	a_m = j_m*t_1
	v_m = a_m * (t_1 + t_2) + v_0
	d_m = (v_m - v_0)*(2*t_1 + t_2 + t_4) + v_0 * (4*t_1 + 2*t_2 + t_4)
	t_j = [[t_1, j_m], [t_2, 0], [t_1, -j_m], [t_4, 0], [t_1,-j_m], [t_2, 0], [t_1, j_m]]
	return {"tick_jerk": t_j, "v_e": v_0, "d": d_m}

"""
type_2: s-curve with no ending tail
	t_1, t_2, t_3, t_4

	a_m = j * t_1 < a
	v_m = a_m (t_1 + t_2) + v_0 = j * t_1 * (t_1 + t_2) + v_0 < v
	d_m = (v_m - v_0)* (t_1 + t_2 /2 + t_4 - 1) + v_0 * (2*t_1 + t_2 + t_4) < d
		j * t_1 * (t_1 + t_2) * (t_1 + t_2/2 + t_4 - 1) + v_0 * (2*t_1 + t_2 + t_4) < d
input: j, a, v, d, v_0 
return: j_m, a_m, v_m, t_1, t_2, t_4
"""
def profile_2(j, a, v, d, v_0 = 0):
	
	t_1_a = math.floor(a / j)
	t_1_v = math.floor((abs(v-v_0)/j)**(1/2))
	# t_1_d
	roots = [ np.real(x) for x in np.roots([ np.sign(v-v_0) *j, np.sign(v-v_0) *-j, 2*v_0, -d]) if np.isreal(x)]
	roots.sort()
	if v >= v_0:
		# positive j
		t_1_d = math.floor(max(roots))
	else:
		# negative j
		if len(roots) <= 1:
			t_1_d = math.floor(1.5*(1+(1+6*v_0/j)**0.5))
		else:
			t_1_d = math.floor(roots[1]) 
	# find t_1_a
	t_1 = min(t_1_a, t_1_v, t_1_d)
	t_1 = max(0, t_1)
	t_2 = 0

	if t_1 == 0:
		if v_0 > 0: 
			t_4 = math.ceil(d/v_0)
			t_4 = max(0, t_4)
		else: 
			t_4 = 0 
		j_m = 0

	else:
		if t_1_a <= min(t_1_v, t_1_d):

			t_2_v = math.floor((abs(v-v_0) / (j*t_1)) - t_1)
			#t_2_d
			roots = [ np.real(x) for x in np.roots([np.sign(v-v_0) *j*t_1 / 2, np.sign(v-v_0)*3*j*(t_1**2)/2-np.sign(v-v_0)*j*t_1+v_0, np.sign(v-v_0) *j*(t_1**3)-np.sign(v-v_0) *j*(t_1**2)+2*v_0*t_1 -d]) if np.isreal(x)]
			roots.sort()
			if v >= v_0:
				# positive j
				t_2_d = math.floor(max(roots))  
			else:
				# negative j
				if len(roots) <= 1:
					t_2_d = math.floor(1+(v_0/(j*t_1))-t_1*3/2)
				else:
					t_2_d = math.floor(roots[0])

			t_2 = min(t_2_v, t_2_d)
			t_2 = max(0, t_2)
	
		t_4 = math.ceil((d-(v_0*(2*t_1+t_2)+np.sign(v-v_0)*j*t_1*(t_1+t_2)*(t_1-1+t_2/2)))/(np.sign(v-v_0) *j*t_1*(t_1+t_2)+v_0))
		t_4 = max(0, t_4)
		j_m = (d - v_0 * (2*t_1 + t_2 + t_4)) / (t_1 * (t_1 + t_2) * (t_1 + t_2/2 + t_4 - 1)) 		
	
	# a_m, v_m, j_m
	a_m = j_m*t_1
	v_m = a_m * (t_1 + t_2) + v_0
	d_m = (v_m - v_0)* (t_1 + t_2 /2 + t_4 - 1) + v_0 * (2*t_1 + t_2 + t_4)

	t_j = [[t_1, j_m], [t_2, 0], [t_1, -j_m], [t_4, 0]]
	#return {"tick_jerk": tick_jerk(t_j), "v_e": v_m, "d": d_m}
	return {"tick_jerk": t_j, "v_e": v_m, "d": d_m}

"""
type_3: only ending tail. Similar to halt
	t_4, t_5 = t_1, t_6 = t_2, t_7 = t_1

	The idea is to get from v_0 to 0 in distance d.
	if d is not achievable then we get to the closest d.
	we also return d_m (acheived d)

	if d is equal to 0 this acts as halt

	fit is True: we adjsut the j to satisfy the d constarint
	fit is false: The d constraint is not necessarily, and keep the j to its limit.

	a_m = j * t_1 < a
	v_s = a_m (t_1 + t_2) = v_0
	d_m = j * t_1 * (t_1 + t_2)* (t_1 + t_2 /2 + t_4 + 1) = d
	v_s*(t_1 + t_2 /2 + t_4 + 1) = d
input: j, a, v_0, d
return: j_m, a_m, d_m, t_1, t_2, t_4
""" 
def profile_3(j, a, d, v_0):
	t_1_a = math.floor(a / j)
	t_1_v = math.floor((v_0/j)**(1/2))
	t_1_d = math.floor((d/v_0)-1)

	# find t_1 and set t_2
	t_1 = min(t_1_a, t_1_v, t_1_d)
	t_1 = max(0, t_1)
	t_2 = 0

	if t_1 == 0:
		if v_0 > 0: 
			t_4 = math.ceil(d/v_0) - 1
			t_4 = max(0, t_4)
		else: 
			t_4 = 0 
		j_m = 0
	else:
		if t_1_a <= min(t_1_v, t_1_d):
			t_2_v = math.floor((v_0 / (j*t_1)) - t_1)
			t_2_d = 2*(t_1_d - t_1)

			# t_2
			t_2 = min(t_2_v, t_2_d)
			t_2 = max(0, t_2)
		
		# t4
		t_4 = math.ceil((d/v_0) - t_1 - t_2/2 - 1)
		t_4 = max(0, t_4)
		j_m = v_0 / (t_1 * (t_1 + t_2))
	#j_m
	a_m = j_m*t_1
	v_m = a_m * (t_1 + t_2)
	d_m = j_m* t_1 * (t_1 + t_2)* (t_1 + t_2 /2 + t_4 + 1) 
	t_j = [[t_4, 0], [t_1, -j_m], [t_2, 0], [t_1, j_m]]
	return {"tick_jerk": t_j, "v_e": 0, "d": d_m}


"""
factor = t2 / t1
"""
def profile_no_jerk(a, vm, d, factor = 1):
	# adjust d
	d = d+ 0.01
	# find n1 = 2*t1 + t2
	n1 = math.floor(vm/a)
	# find t4
	t4 = math.floor(d/vm) - n1	
	# adjust vm
	if vm * n1 > d:
		vm = d / n1
	# find t1, t2
	# make sure t1 is not zero
	t1 = math.floor(n1/(2+factor))+ 1
	t2 = math.floor(factor*t1)
	"""
	vm = t1*(t1+t2)*jm
	d = jm* t1*(t1+t2) * (2*t1 + t2 + t4)
	jm = d/ (t1*(t1+t2) * (2*t1 + t2 + t4)) 	
	"""
	# find jerk
	jm = d/ (t1*(t1+t2) * (2*t1 + t2 + t4))
	# tick and jerk
	t_j = [[t1, jm], [t2, 0], [t1, -jm], [t4, 0], [t1,-jm], [t2, 0], [t1, jm]]
	
	return {"tick_jerk": t_j, "d": d, "v_e": 0}


"""
given the number of ticks and intial values, find final a, v, q

def tavq(t, j0, a0 = 0, v0 = 0, q0 = 0):
	q0 + = t*v0 + t*(t-1)*a0/2 + t*(t-1)(t-2)*j0/6
	v0 + = t*a0 + t*(t-1)*j0/2
	a0 + = t*j0
	return a0, v0, q0
"""

"""
given initial j and a, and final a, find ticks and jerk that get you there
"""
def a_inv(j0, a0 , am):
	if am == a0:
		t = 0
		j = j0
	else:
		t = 1 + math.floor((am - a0) / jm)
		j = (am - a0) / t

	return t, j


"""
type_4: 
	start from a0, v0 and get to v while not exceeding a and j, and total distance d
	it will return tick jerk and final accel, and v (ae, ve)
"""	

"""
def profile_4(jm, am, vm, d, a0 = 0, v0 = 0):
	# do all the initial filtering
	if vm > v0 and am >= a0:
		
		t3, j3 = a_inv(jm, a0, 0)
		a3, v3, q3 = tavq (t3, j3, a0, v0, 0)
		
		# decrease a to zero
		if vm <= v3:
			# achieve d before changing j
			if q3 > d
				# find t that achieve d
				# t*v0 + t*(t-1)*a0/2 + t*(t-1)(t-2)*j3/6 == d
				return True
			# touch a==0 
			return True

		
		# increasing a to am
		t1, j1 = a_inv(jm, a0, am)
		a1, v1, q1 = tavq (t1, j1, a0, v0, 0)
		# decreasing a to 0
		t3, j3 = a_inv(jm, a1, 0)
		a3, v3, q3 = tavq (t3, j3, a1, v1, q1)
		
		# find t1 and t3 s.t it will achieve vm
		if vm <= v3:
			beta = -a0/j3
			alpha = -j1/j3

			# find t1 s.t
			# solve(alpha+1)*j1/2 (t1**2) + (alpha+1)*a0*t1 + (v0-vm+a0/2+(beta*a0/2)) ==0
			# equations
			v1 = v0 + t1*a0 + (t1*(t1-1))*j1/2 
			v3 = 
			# asign value # solve based on achieving to zero a
			t1 = math.floor(t1+2)
			j1 = 2*(v1-v0 - t1*a0)/(t1*(t1-1))
			a1 = a0 + t1*j1  
			t3 = math.floor(2*(vm - v1 - a1/2)/a1)
		# find t2 s.t it will achieve vm
"""


def plot(a, v, t, c = "r"):
	plt.figure(1)
	
	plt.subplot(311)
	plt.plot(a, c+'o-')
	plt.title("a_end "+ str(a[-1]) )
	
	plt.subplot(312)
	plt.plot(v, c+'o-')
	plt.title("v_end " + str(v[-1]))

	plt.subplot(313)
	plt.plot(t, c+'o-')
	plt.title("t_end "+ str(t[-1]))

	plt.show()


def main_test_sadegh_code():
	javn = []
	filepath = 'test_500.txt'
	with open(filepath) as fp:
		for line in fp:
			x = line.split()
			#print(x)
			x = [float(y) for y in x]
			try:
				print(x[5]/x[0], x[6]/x[1], x[7]/x[2])
			except Exception as e:
				pass
			javn.append({"j": x[0], "a": x[1], "v": x[2], "n": int(x[4])})
	
	ticks = [x["n"] for x in javn]
	jerk = [x["j"] for x in javn]
	a, v, q = tick(jerk, ticks, 0)
	#plot(a, v, q)

if __name__ == '__main__':
	
	j = 2.0e-10
	a = 17.0e-6
	v = 0.01
	d = 10
	v_0 =0
	factor = 1

	result = profile_no_jerk(a, v, d, factor)
	#result = profile_1(j, a, v, d)
	print(result)
	jerk = [x[1] for x in result["tick_jerk"]]
	ticks = [x[0] for x in result["tick_jerk"]]
	print("number of ticks: ", sum(ticks))
	a, v, q = tick(jerk, ticks, v_0)
	plot(a, v, q)
	print(len(q))
	