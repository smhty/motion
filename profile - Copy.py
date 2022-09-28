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

concatenate different motions together and form a smooth motion profile

motion_0 (j_0, a_0, v_0, d_0): start ramp up and achieve the max speed, and exit this motion with max speed,
d_0 and zero laceration

motion_1 (j_1, a_1, v_1, d_1): start ramp up and achieve the max speed, and exit this motion with max speed,
d_1 and zero acceleration 

...
last motion
motion_e (j_e, a_e, v_e, d_e): start ramp up and achieve the max speed, and exit this motion with max speed,
d_1 and zero acceleration 

type_2: s-curve with no ending tail
	t_1, t_2, t_3, t_4

type_3: 
	t_1, t_2, t_3, t_4, t_5, t_6, t_7
	start with v_s
	get to v_m during t_1, t_2, t_3
	stay in v_m during t_4
	get to zero during t_5, t_6, t_7 
	
	case_1: go up
	t1: j
	a_m = t_1*j
	v = v_0 + t_1*(t_1-1)*j/2
	q = t_1*v_0

	t2: 0
	a = t_1*j
	v = v_0 + t_1*(t_1-1)*j/2 + t_2*t_1*j
	q = 

	t3: -j
	a = 0
	v_m = v_0 + a_m*(t_1+t_2)
	q = (v_m-v_0)*(t_1+t_2/2-1)+v_0*(2*t_1+t_2)

	t4: 0
	a = 0
	v = v_0 + t_1*(t_1+t_2)*j
	q = (v_m-v_0)*(t_1+t_2/2+t_4-1)+v_0*(2*t_1+t_2+t_4)

	t5: -j
	a = -t_5*j
	v = v_0 + t_1*(t_1+t_2)*j - t_5*(t_5-1)*j/2
	q = 

	t6: 0
	a = -t_5*j
	v = v_0 + t_1*(t_1+t_2)*j -t_5*(t_5-1)*j/2 -t_5*t_6*j
	q = 

	t7: j
	a = 0
	v = v_0 + t_1*(t_1+t_2)*j - t_5*(t_5+t_6)*j 
	q = 
"""

import math
import matplotlib.pyplot as plt
import time
import numpy as np



def tick(jerk, ticks, a_0 = 0, v_0 = 0, q_0 = 0):
	q = [q_0]
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
	d_m = (v_m - v_0)*(2*t_1 + t_2 + t_4) + v_0 * (4*t_1 + 2*t_2 + t_4) = 
	j * t_1 * (t_1 + t_2) * (2*t_1 + t_2 + t_4) + v_0 * (4*t_1 + 2*t_2 + t_4) < d
	
input: j, a, v, d, v_0 
return: j_m, a_m, v_m, t_1, t_2, t_4
"""
class type_1(object):
	"""docstring for type_1"""
	def __init__(self):
		super(type_1, self).__init__()
		
	def prm_j_a_v(self, j, a, v, d, v_0 = 0):
		
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
	v: a number between 0-1
	a: a number between 0-1
	t: total number of ticks
	d: total distance

	v=0 min v_m
	v=1 max v_m
	
	a = 0 min a_m
	a = 1 max a_m

	t_4 = math.floor(t*(-1 + 2/(1+v)))
	# v_m = 2 * (d - t*v_0) / (t + t_4)

	t_1 = math.floor((t-t_4)/2 * (1-1/(1+a)))
	t_2 = math.floor((t - t_4 - 4*t_1)/2)
	t_4 = math.floor(t - 4*t_1 - 2*t_2) # update t_4

	j = d /  (t_1 * (t_1 + t_2) * (2*t_1 + t_2 + t_4))
	a_m = j * t_1
	v_m = a_m *(t_1 + t_2)
	"""
	def prm_t_v_a(self, d, t, v, a):
		t_4 = math.floor(t*(-1 + 2/(1+v)))
		# v_m = 2 * (d - t*v_0) / (t + t_4)

		t_1 = math.floor((t-t_4)/2 * (1-1/(1+a)))
		t_2 = math.floor((t - t_4 - 4*t_1)/2)
		t_4 = math.floor(t - 4*t_1 - 2*t_2) # update t_4

		j_m = d /  (t_1 * (t_1 + t_2) * (2*t_1 + t_2 + t_4))
		a_m = j_m * t_1
		v_m = a_m *(t_1 + t_2)
		d_m = v_m *(2*t_1 + t_2 + t_4)

		t_j = [[t_1, j_m], [t_2, 0], [t_1, -j_m], [t_4, 0], [t_1,-j_m], [t_2, 0], [t_1, j_m]]
		return {"tick_jerk": t_j, "v_e": 0, "d": d_m}

	"""
	j: jerk
	d: travel distance
	t_d: total number of ticks
	p_a: portion of td that the robot is accelerating and decelerating
	v_0: initial velocity  

	t_d = 4*t_1 + 2*t_2 + t_4
	p_a = (4*t_1 + 2*t_2) / t_d

	t_4 = t_d * (1-p_a)

	v_m = (d - v0 * t_d*p_a/2) / (t_d*(1-p_a/2))
	"""
	def prm_time(self, j, d, t_d, p_a=0.5, v_0 = 0):
		m = p_a*t_d

		# t_4
		t_4 = t_d - m

		# v_m from d
		v_m = (d - v_0 * m /2)/(t_4 + m / 2)

		r = (v_m - v_0)/j
		if abs(r) >= m**2 / 16: # increase j
			t_1 = m/4 # t_1 >=1
		else:
			t_1 = m/4 - (1/2)* math.sqrt(m**2 / 4 - 4*abs(r)) # pick the smaller solution
		
		# t_1
		t_1 = math.floor(t_1)

		# t_2
		t_2 = math.floor(m/2 - 2*t_1)

		#update t_4 based on t_d
		t_4 = math.floor(t_d - 4* t_1 - 2*t_2)

		# j based on d
		j_m = (d - v_0 * (4*t_1 + 2*t_2 + t_4))/ (t_1 * (t_1 + t_2) * (2*t_1 + t_2 + t_4))
		
		t_j = [[t_1, j_m], [t_2, 0], [t_1, -j_m], [t_4, 0], [t_1,-j_m], [t_2, 0], [t_1, j_m]]
		return {"tick_jerk": t_j, "v_e": v_0, "d": d}



"""
type_2: s-curve with no ending tail
	t_1, t_2, t_3, t_4

	a_m = j * t_1 < a
	v_m = a_m (t_1 + t_2) + v_0 = j * t_1 * (t_1 + t_2) + v_0 < v
	d_m = (v_m - v_0)* (t_1 + t_2 /2 + t_4 - 1) + v_0 * (2*t_1 + t_2 + t_4) < d
		j * t_1 * (t_1 + t_2) * (t_1 + t_2/2 + t_4 - 1) + v_0 * (2*t_1 + t_2 + t_4) < d
"""
class type_2(object):
	"""docstring for type_2"""
	def __init__(self):
		super(type_2, self).__init__()

	"""
	v: a number between 0-1
	a: a number between 0-1
	t: total number of ticks
	d: total distance

	v=0 min v_m
	v=1 max v_m
	
	a = 0 min a_m
	a = 1 max a_m

	
	if d-v_0*t >=0:
		alpha = (t-1)*(t-2) / (t*(1+v)-2)
		t_4 = math.floor(2*alpha -t + 1)

	if d-v_0*t < 0:
		denom = 1/2 - v_0 * (t-1) / (v*(d-v_0))
		nom = ((d*(t-1) -v_0*(t-1)*(t/2 + 1)) / (v*(d-v_0))) -t/2 + 1
		t_4 = math.floor(nom/denom)

	t_1 = math.floor((t-t_4)*(1-1/(1+a)))
	t_2 = math.floor((t - t_4 - 2*t_1))
	t_4 = math.floor(t - 2*t_1 - t_2) # update t_4 

	"""
	def prm_t_v_a(self, d, t, v, a, v_0=0):
		if d-v_0*t >=0:
			alpha = (t-1)*(t-2) / (t*(1+v)-2)
			t_4 = math.floor(2*alpha -t + 1)

		if d-v_0*t < 0:
			v_min = max(0, d-v_0*(1+t/2))
			print("v_min: ", v_min)
			v_max = (d-v_0) / (t-1)
			t_4 = 2*(d-v_0*t) / (v_min + v*(v_max-v_min) -v_0) -t+2 
			t_4 = math.floor(t_4)

		t_1 = math.floor((t-t_4)*(1-1/(1+a)))
		t_2 = math.floor((t - t_4 - 2*t_1))
		t_4 = math.floor(t - 2*t_1 - t_2) # update t_4 

		j_m = (d - v_0 * (2*t_1 + t_2 + t_4)) / (t_1 * (t_1 + t_2) * (t_1 + t_2/2 + t_4 - 1))
		a_m = j_m * t_1
		v_m = a_m * (t_1 + t_2) + v_0 

		t_j = [[t_1, j_m], [t_2, 0], [t_1, -j_m], [t_4, 0]]
		return {"tick_jerk": t_j, "v_e": v_m, "d": d}


	"""
	j: jerk
	d: travel distance
	t_d: total number of ticks
	p_a: portion of td that the robot is accelerating and decelerating
	v_0: initial velocity 

	t_d =  2*t_1 + t_2 + t_4
	p_a = (2*t_1 + t_2) / t_d
	"""
	def prm_time(self, j, d, t_d, p_a=0.5, v_0 = 0):
		m = p_a * t_d
		
		# check v_m > 0
		if d - v_0 * (1 + m /2) <= 0:
			m = 1.8 * (d/v_0 - 1)
			p_a = m / t_d

 		# t_4
		t_4 = t_d - m
		v_m = (d - v_0 * (1 + m /2))/(t_4 + m / 2 - 1)

		# ratio
		r = (v_m - v_0)/j
		if abs(r) >= m**2 / 4: # increase j
			t_1 = m/2 # t_1 >=1
		else:
			t_1 = m/2 - (1/2)* math.sqrt(m**2 - 4*abs(r)) # pick the smaller solution
		
		# t_1
		t_1 = math.floor(t_1)
			
		# t_2
		t_2 = math.floor(m- 2*t_1)
		
		#update t_4 based on t_d
		t_4 = math.floor(t_d - 2* t_1 - t_2)

		# j based on d
		j_m = (d - v_0 * (2*t_1 + t_2 + t_4))/ (t_1 * (t_1 + t_2) * (t_1 + t_2/2 + t_4 - 1))
		
		v_m = (d - v_0 * (1 + m /2))/(t_4 + m / 2 - 1) 
		
		t_j = [[t_1, j_m], [t_2, 0], [t_1, -j_m], [t_4, 0]]
		return {"tick_jerk": t_j, "v_e": v_m, "d": d, "t_d": t_d, "p_a": p_a}


	"""
	input: j, a, v, d, v_0 
	return: j_m, a_m, v_m, t_1, t_2, t_4
	"""
	def prm_j_a_v(self, j, a, v, d, v_0 = 0):
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
	we also return d_m (achieved d)

	if d is equal to 0 this acts as halt

	fit is True: we adjust the j to satisfy the d constraint
	fit is false: The d constraint is not necessarily, and keep the j to its limit.

	a_m = j * t_1 < a
	v_s = a_m (t_1 + t_2) = v_0
	d_m = j * t_1 * (t_1 + t_2)* (t_1 + t_2 /2 + t_4 + 1) = d
	v_s*(t_1 + t_2 /2 + t_4 + 1) = d
input: j, a, v_0, d
return: j_m, a_m, d_m, t_1, t_2, t_4
""" 
class type_3(object):
	"""docstring for type_3"""
	def __init__(self):
		super(type_3, self).__init__()

	"""
	run the motion in positive (or negative), zero and negative accel depending on the initial parameter
	the max accel on positive and negative are the same 
	"""
	def prm_t_v_a(self, d, t, v, a, v_0=0):
		t_tmp = t

		if v_0 * t < 0.8 * d:
			t_tmp = math.floor(1.2*d / v_0)

		# case 1 t is smaller		
		elif v_0*t/2 > d:
			t_tmp = math.floor(2*d / v_0)

		# make sure t_4 >= 0
		t_4 = math.floor(2*(d/v_0 -1) - t_tmp)
		if t_4 < 0:
			t_4 = 0
		
		# m = 2*t_1 + t_2
		m = t_tmp - t_4
		t_1 = math.floor(m/(1+a) - m/2)
		t_2 = math.floor(m-2*t_1)
		t_4 = math.floor(t_tmp - 2*t_1 - t_2)

		j_m = d / (t_1 * (t_1 + t_2)* (t_1 + t_2 /2 + t_4 + 1))
		a_m = j_m * t_1
		v_e = 0 

		t_j = [[t_4, 0], [t_1, -j_m], [t_2, 0], [t_1, j_m]]
		return {"tick_jerk": t_j, "v_e": v_e, "d": d, "t_e": t_tmp, "t": t}


	def prm_t_v_a_copy(self, d, t, v, a, v_0=0):
		# case 1
		if v_0*t - d <= 0:
			v_m_min = d/t
			theta_y = math.sqrt((2*v_0*t-4*d)+4*(v_0*t)**2) + (2*v_0*t - 4*d)  
			theta_x = 2*v_0**2
			tan_theta = theta_y / theta_x
			v_m_max = (v_0 * tan_theta + t) / (2*tan_theta)

		if d/t >= v_0/2:
			v_m_min = d/t
			v_m_max = v_0
			v_m = v_min + v*(v_max-v_min)

			# t_4
			t_4 = math.floor(t*(2*v_m_min - v_0)/(2*v_m - v_0))
			T_5 = math.floor(2*t_4*(v_m/v_0)*(v_m - v_m_min)/(2*v_m_min - v_0))
			T_1 = math.floor(t-t_4-T_5)
		else:
			v_m_min = 0
			v_m_max = d/t

		t_1 = math.floor((t-t_4)*(1-1/(1+a)))
		t_2 = math.floor((t - t_4 - 2*t_1))
		t_4 = math.floor(t - 2*t_1 - t_2) # update t_4 

		j_m = (d - v_0 * (2*t_1 + t_2 + t_4)) / (t_1 * (t_1 + t_2) * (t_1 + t_2/2 + t_4 - 1))
		a_m = j_m * t_1
		v_m = a_m * (t_1 + t_2) + v_0 

		t_j = [[t_1, j_m], [t_2, 0], [t_1, -j_m], [t_4, 0]]
		return {"tick_jerk": t_j, "v_e": v_m, "d": d}

		
	def prm_j_a_v(self, j, a, d, v_0):
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
	j: jerk
	d: travel distance
	td: total number of ticks
	v_0: initial velocity 

	
	We use two parts, t1 and t2, with -j1 and j2 jerk respectively
	starting from vo, at the end of t1 we have
	q[t_1 - 1] = t_1*v_0 - t_1*(t_1-1)(t_1-2)*j_1/6
	v[t_1 - 1] = v_0 - t_1*(t_1-1)*j_1/2
	a[t_1 - 1] = -t_1*j_1

	if t_d * v_0 <= d:

	"""
	def prm_time(self, j, d, t_d, p_a=0.5, v_0 = 0):
		if t_d * v_0 <= d:
			t_1 = math.floor(t_d/4)
	
			v_m = 2*(d - (v_0* t_d/4))/t_d
			j_1 = (v_m - v_0) / (t_1**2)
			j_2 = v_m / (t_1**2)

		else:
			v_m = v_0 / 2
			t_1 = 123		  

			# can not make the time
			return {"tick_jerk": [], "error": 1}
		m = t_d - d/v_0 +1
		t_4 = math.floor(t_d - 2*m)
		print(t_4)
		if v_0 / j >= m**2:
			t_1 = m
		else:
			t_1 = m - math.sqrt(m**2 - v_0/j)

		t_1 = math.floor(t_1)

		t_2 = math.floor(t_d - t_4 - 2*t_1)

		j_m = d / (t_1*(t_1+t_2)*(t_1 + t_2 /2 + t_4 + 1))

		t_j = [[t_4, 0], [t_1, -j_m], [t_2, 0], [t_1, j_m]]
		return {"tick_jerk": t_j, "v_e": 0, "d": d}


"""
given the number of ticks and initial values, find final a, v, q
a[t-1] = a_0 + t*j0

"""
def tavq(t, j0, a0 = 0, v0 = 0, q0 = 0):
	q0 += t*v0 + t*(t-1)*a0/2 + t*(t-1)(t-2)*j0/6
	v0 += t*a0 + t*(t-1)*j0/2
	a0 += t*j0
	return a0, v0, q0


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


def plot(a, v, t, c = "r"):
	plt.figure(1)
	
	plt.subplot(311)
	plt.plot(a, c+'-')
	plt.title("a_end "+ str(a[-1]) )
	"""
	#plt.title("Accelaration" )
	#plt.xticks([])
	#plt.yticks([])
	"""

	plt.subplot(312)
	plt.plot(v, c+'-')
	plt.title("v_end " + str(v[-1]))
	"""
	plt.title("Velocity" )
	plt.xticks([])
	plt.yticks([])
	"""


	plt.subplot(313)
	plt.plot(t, c+'-')
	plt.title("t_end "+ str(t[-1]))
	"""
	plt.title("Distance" )
	plt.xticks([])
	plt.yticks([])
	"""
	plt.show()

def plot_2(x, a, v, t, c = "ro"):
	#plt.figure(1)
	
	plt.subplot(311)
	plt.plot(x,a, c)
	plt.title("a_end "+ str(a[-1]) )
	"""
	#plt.title("Accelaration" )
	#plt.xticks([])
	#plt.yticks([])
	"""

	plt.subplot(312)
	plt.plot(x, v, c)
	plt.title("v_end " + str(v[-1]))
	"""
	plt.title("Velocity" )
	plt.xticks([])
	plt.yticks([])
	"""


	plt.subplot(313)
	plt.plot(x, t, c)
	plt.title("t_end "+ str(t[-1]))
	"""
	plt.title("Distance" )
	plt.xticks([])
	plt.yticks([])
	"""
	#plt.show()	




def main_1():
	j = 1.0e-15
	a = 1.0e-12
	v =1.0e-7
	v_0 =1.0e-2
	d = 300
	t_d = 5.0e4
	p_a= 1

	"""
	x = time.time()
	p = profile()
	slc = p.slice(t, j, a, v_0, d)
	print(time.time() - x)

	print(p.q_n(slc["tick_jerk"][0][0], slc["tick_jerk"][0][1], a, v_0), d)
	"""
	factor = 1

	#result = profile_no_jerk(a, v, d, factor)
	#result = type_3().prm_time(j, d, t_d, p_a, v_0)
	#result = type_2().prm_j_a_v(j, a, v, d, v_0)
	result = type_3().prm_t_v_a(d, t_d, 0.8, 0.5, v_0)
	print(result)
	jerk = [x[1] for x in result["tick_jerk"]]
	ticks = [x[0] for x in result["tick_jerk"]]
	a, v, q = tick(jerk, ticks, v_0=v_0)
	plot(a, v, q)


def main_cont():
	a0 = 0
	v0 = 0
	q0 = 0
	a = []
	v = []
	q = []
	d = 0
	
	path_list = [
		{"jm": 7.29e-12, "am": 2.43e-07, "vm": 0.0045, "d":   47.7231, "color": "r--"},
		{"jm": 7.29e-12, "am": 2.43e-07, "vm": 0.0045, "d":   790.8024, "color": "b--"},
		{"jm": 7.29e-12, "am": 2.43e-07, "vm": 0.002, "d":   47.7231, "color": "r--"},
		{"jm": 7.29e-12, "am": 2.43e-07, "vm": 0.002, "d":   790.8024, "color": "b--"},
		{"jm": 7.29e-12, "am": 2.43e-07, "vm": 0, "d":   790.8024, "color": "r--"},
		]

	"""
	path_list = [
		{"jm": 1.0e-10, "am": 5.0e-6, "vm": 0.1, "d": 500, "color": "r--"},
		{"jm": 5.0e-10, "am": 5.0e-6, "vm": 0.3, "d": 6000, "color": "b--"},
		{"jm": 1.0e-10, "am": 5.0e-6, "vm": 0.05, "d": 3000, "color": "r--"},
		{"jm": 1.0e-10, "am": 5.0e-6, "vm": 0.01, "d": 8000, "color": "b--"},
		{"jm": 1.0e-10, "am": 5.0e-6, "vm": 0, "d": 8000, "color": "r--"},
	]
	"""
	plt.figure(1)
	p = profile()
	for path in path_list:
		d += path["d"]  
		result = p.cont(a0, v0, path["jm"], path["am"], path["vm"], d - q0)
		print("result: ", result)
		p.depth = 0
		jerk = [x[1] for x in result["tick_jerk"]]
		ticks = [x[0] for x in result["tick_jerk"]]
		_a, _v, _q = tick(jerk, ticks, a0, v0, q0)

		# append
		x = list(range(len(a), len(a) + len(_a)))
		a += _a
		v += _v
		q += _q

		# update
		a0 = result["ae"]
		v0 = result["ve"]
		q0 += result["de"]

		plot_2(x, _a, _v, _q, path["color"])
	
	plt.show()


def main_test():
	a0 = -3.8098097394908223e-06
	v0 = 0.19427840220378942
	q0 = 0
	a = []
	v = []
	q = []
	d = 0

	path_list = [
		{"jm": 5.0e-10, "am": 5.0e-5, "vm":0, "d": 5100.179543354534, "color": "r--"},
		]
	plt.figure(1)
	p = profile()
	for path in path_list:
		d += path["d"]  
		result = p.cont(a0, v0, path["jm"], path["am"], path["vm"], d - q0)
		print(result)
		p.depth = 0
		jerk = [x[1] for x in result["tick_jerk"]]
		ticks = [x[0] for x in result["tick_jerk"]]
		_a, _v, _q = tick(jerk, ticks, a0, v0, q0)

		# append
		x = list(range(len(a), len(a) + len(_a)))
		a += _a
		v += _v
		q += _q

		# update
		a0 = result["ae"]
		v0 = result["ve"]
		q0 += result["de"]

		plot_2(x, _a, _v, _q, path["color"])
	
	plt.show()
if __name__ == '__main__':
	#main_cont()
	#main_test()
	main_1()