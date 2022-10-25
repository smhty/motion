import math
import time
import numpy as np
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

"""
[t_0, 0] -> d_0 flat
	T_1 -> accel a v_0 to v_m 
	T_4 -> flat 
	T_1 -> accel a V_M to v_e
[t_m, 0] -> d_2 flat	

if total time t is not enough then we increase it
if a_avg is not within the range then we change it 

conditions 
a_avg > 0

D = v_m*(t_4+t_1+t_2/2+t_5+t_6/2) + v_0*(t_1+t_2/2+1)+v_e*(t_5+t_6/2-1)
"""
def profile(d_list, t, a_avg, v_0, v_e, j):
	# threshold
	thr = 0.0001

	# t_0
	if v_0 == 0:
		t_0 = 0
	else:
		t_0 = math.floor(d_list[0]/v_0)

	# t_e
	if v_e == 0:
		t_e = 0
	else:
		t_e = math.floor(d_list[2]/v_e)
	

	# adjust D: from v_0 to v_e
	D = d_list[1] + (d_list[0]-t_0*v_0) + (d_list[2]-t_e*v_e)

	# adjust t
	t_tmp = math.floor(t-t_0-t_e)

	# time is to short
	if t_tmp <=0:
		t_tmp = 4

	"""
	adjust a_avg to a feasible accel
	"""
	# minimum value of a_avg
	_a_avg = max(a_avg, abs(v_0-v_e)/t_tmp)

	# d_max
	v_max = (_a_avg*t_tmp+v_0+v_e)/2
	_T_1 = (v_max-v_0)/(2*_a_avg)
	_T_3 = (v_max-v_e)/(2*_a_avg) 
	d_max = v_max*(_T_1+_T_3)+v_0*(_T_1+1)+v_e*(_T_3-1)

	# d_min 
	v_min = (-_a_avg*t_tmp+v_0+v_e)/2
	if v_min < 0:
		v_min = 0
		_T_1 = v_0/(2*_a_avg)
		_T_3 = v_e/(2*_a_avg)
	else:
		_T_1 = (v_0-v_min)/(2*_a_avg)
		_T_3 = (v_e-v_min)/(2*_a_avg) 
	d_min = v_min*(_T_1+_T_3)+v_0*(_T_1+1)+v_e*(_T_3-1)

	# check D and adjust _a_avg based on that
	if D > d_max:
		# increase _a_avg
		A = (t_tmp**2)/4
		B = t_tmp*(v_e+v_0)/2 + v_0-v_e-(1+thr)*D
		C = -(v_0-v_e)**2/4
		delta = B**2-4*A*C
		_a_avg = (-B+math.sqrt(delta))/(2*A) # pick the largest root


	elif D < d_min:
		# increase _a_avg
		A = -(t_tmp**2)/4
		B = t_tmp*(v_e+v_0)/2 + v_0-v_e-(1-thr)*D
		C = (v_0-v_e)**2/4
		delta = B**2-4*A*C
		_a_avg = (-B-math.sqrt(delta))/(2*A) # pick the largest root 

		# check if v_min is negative or not
		if (v_0+v_e - _a_avg*t_tmp)/2 <0:
			_a_avg = (v_0**2 + v_e**2)/(2*((1-thr)*D+v_e-v_0))

	# init v_m
	v_m = None
	
	"""
	case 1 v_m >= max
	"""
	while v_m == None:
		A = -1/_a_avg
		B = (v_0 + v_e)/_a_avg + t_tmp
		C = -(v_0**2 + v_e**2)/(2*_a_avg) + v_0-v_e-D   
		delta = B**2-4*A*C

		# check delta
		if delta < 0:
			break
		
		# pick the max root
		roots = [(-B+math.sqrt(delta))/(2*A),(-B-math.sqrt(delta))/(2*A)]
		
		# check if v_m and t_4 is valid
		for root in roots:
			if root >= max(v_0, v_e): # v_m condition
				if t_tmp - (2*root-v_e-v_0)/_a_avg >= 0:
					v_m = root
		
		# break the loop
		break

	"""
	case 2 v_m <= min
	"""
	while v_m == None:
		A = 1/_a_avg
		B = -(v_0 + v_e)/_a_avg + t_tmp
		C = (v_0**2 + v_e**2)/(2*_a_avg) + v_0-v_e-D   
		delta = B**2-4*A*C
		# check delta
		if delta < 0:
			break
		
		# all possible solutions
		sols = [(-B+math.sqrt(delta))/(2*A),(-B-math.sqrt(delta))/(2*A)]
		for sol in sols:
			if sol >=0 and sol <= min(v_0, v_e):
				v_m = sol
				break

		# break the loop
		break

	"""
	case 3 v_e <= v_m <= v_0
	"""
	while v_0 >= v_e and v_m == None:
		A = (-v_0 + v_e)/_a_avg + t_tmp
		B = (v_0**2 - v_e**2)/(2*_a_avg) + v_0-v_e-D   
		
		sol = -B/A
		
		if sol >=0 and sol <= v_0 and sol >= v_e:
			v_m = sol
		
		# break the loop
		break

	"""
	case 4 v_0 <= v_m <= v_e
	"""
	while v_e >= v_0 and v_m == None:
		A = (v_0 - v_e)/_a_avg + t_tmp
		B = (-v_0**2 + v_e**2)/(2*_a_avg) + v_0-v_e-D   
		
		sol = -B/A
		
		if sol >=0 and sol <= v_e and sol >= v_0:
			v_m = sol
		
		# break the loop
		break

	# check if v_m exists
	if v_m == None:
		return False

	# t_1 and t_2
	if v_0 == v_m:
		t_1 = math.floor(0)
		t_2 = math.floor(0)
	else:
		T_1 = abs(v_m-v_0)/(2*_a_avg) # t_1+t_2/2 
		#t_1 = max(1, math.floor(2*T_1*(1-1/(1+a))))
		t_1 = max(1, math.floor(T_1 + j*(1-T_1)))
		t_2 = max(0, math.floor(2*T_1 - 2*t_1))

	# t_5 and t_6
	if v_e == v_m and t_1 != 0 and t_2 != 0:
		t_5 = math.floor(0)
		t_6 = math.floor(0)
	else:
		T_3 = abs(v_m-v_e)/(2*_a_avg) # t_5+t_6/2
		#t_5 = max(1, math.floor(2*T_3*(1-1/(1+a))))
		t_5 = max(1, math.floor(T_3 + j*(1-T_3)))
		t_6 = max(0, math.floor(2*T_3 - 2*t_5))

	# t_4
	t_4 = t_tmp - 2*t_1 - t_2 - 2*t_5 - t_6

	# adjust v_m
	v_m = (D - v_0*(t_1+(t_2/2)+1) - v_e*(t_5+(t_6/2)-1)) / (t_tmp - t_1-(t_2/2) - t_5-(t_6/2))
	
	# j_1
	if t_1 == 0:
		j_1 = 0
	else:
		j_1 = (v_m - v_0)/ (t_1 * (t_1 + t_2))
	
	# j_2
	if t_5 == 0:
		j_2 = 0
	else:
		j_2 = (v_e - v_m)/ (t_5 * (t_5 + t_6))
	
	d_m = v_m*(t_1 + t_2/2 + t_5 + t_6/2+ t_4) + v_0*(t_1 + t_2/2 + 1) + v_e*(t_5 + t_6/2 - 1)
	t_m = t_0+2*t_1+t_2+t_4+2*t_5+t_6+t_e
	t_j = [[t_0, 0], [t_1, j_1], [t_2, 0], [t_1, -j_1], [t_4, 0], [t_5,j_2], [t_6, 0], [t_5, -j_2], [t_e, 0]]

	return {"tick_jerk": t_j, "v_0": v_0, "v_e": v_e, "d": sum(d_list), "d_e": d_m+v_0*t_0+v_e*t_e, "a_avg": a_avg, "a_avg_e": _a_avg, "t": t, "t_m": t_m}


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
def type_1_j_a_v(j, a, v, d, v_0 = 0):
	
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

def type_1_t_v_a(d, t, v, a):
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
type_2: s-curve with no ending tail
	t_1, t_2, t_3, t_4

	a_m = j * t_1 < a
	v_m = a_m (t_1 + t_2) + v_0 = j * t_1 * (t_1 + t_2) + v_0 < v
	d_m = (v_m - v_0)* (t_1 + t_2 /2 + t_4 - 1) + v_0 * (2*t_1 + t_2 + t_4) < d
		j * t_1 * (t_1 + t_2) * (t_1 + t_2/2 + t_4 - 1) + v_0 * (2*t_1 + t_2 + t_4) < d
"""

"""
v: a number between 0-1
a: a number between 0-1
t: total number of ticks
d: total distance list [d_0, d_1, d_2]

d_0: no accel -> [t_0, 0]
d_1: accel and final vel -> [t_1, j_m], [t_2, 0], [t_1, -j_m], [t_4, 0]
d_2: no accel -> [t_5, 0]

tick_ jerck: [[t_0, 0], [t_1, j_m], [t_2, 0], [t_1, -j_m], [t_4, 0], [t_5, 0]]
"""
def type_2_t_v_a(d_list, t, v, a, v_0=0):
	# case
	sgn = 1
	if t*v_0 > sum(d_list):
		sgn = -1
	elif t*v_0 == sum(d_list):
		return {"tick_jerk": [[t,0]], "v_e": v_0, "d": sum(d_list), "t_e": t, "t_init": t}

	#t_0
	t_0 = math.floor(d_list[0]/v_0)

	# remaining t
	t_tmp = t - t_0

	# update d
	#d_list[1] += d_list[0]-(t_0*v_0)

	# v_m, and find the largest root
	b = -(t_tmp+2)*v_0+d_list[2]+2*d_list[1]
	root_1 = (b+math.sqrt(b**2+4*(t_tmp-2)*v_0*d_list[2]))/(2*(t_tmp-2))
	root_2 = (b-math.sqrt(b**2+4*(t_tmp-2)*v_0*d_list[2]))/(2*(t_tmp-2))
	beta = max(0, root_1, root_2)
	if sgn > 0:
		v_m_min = max(0, v_0, (d_list[2]+d_list[1]-v_0)/(t_tmp-1))
		v_m_max = beta
	else:
		v_m_min = beta
		v_m_max = min(v_0, (d_list[2]+d_list[1]-v_0)/(t_tmp-1))

	if v_m_min > v_m_max:
		return False

	v_m = v_m_min + v*(v_m_max - v_m_min)
	t_5 = math.floor(d_list[2] / v_m)
	t_4 = math.floor(((2-t_tmp)*v_m-(t_tmp+2)*v_0+d_list[2]+2*d_list[1] + v_0*d_list[2]/v_m)/(v_m-v_0))
	M_2 = math.floor(t_tmp-t_5-t_4)

	# a_m
	t_1 = math.ceil((M_2)*(1-1/(1+a)))
	t_2 = math.floor(M_2-2*t_1)

	j_m = (sum(d_list)-v_0*t)/(t_1*(t_1+t_2)*(t_1+t_2/2+t_4+t_5- 1))
	a_m = j_m * t_1
	v_m = a_m * (t_1 + t_2) + v_0 

	t_j = [[t_0, 0], [t_1, j_m], [t_2, 0], [t_1, -j_m], [t_4, 0], [t_5, 0]]
	return {"tick_jerk": t_j, "v_e": v_m, "d": sum(d_list), "t_e": t_0 + 2*t_1+t_2+t_4+t_5, "t_init": t}


"""
input: j, a, v, d, v_0 
return: j_m, a_m, v_m, t_1, t_2, t_4
"""
def type_2_j_a_v(j, a, v, d, v_0 = 0):
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

"""
run the motion in positive (or negative), zero and negative accel depending on the initial parameter
the max accel on positive and negative are the same 
"""
"""

t_0 with v_0 and d_0

d_1 in [t_1, 0] [t_2, -j], [t_3, j]
"""
def type_3_t_v_a(d_list, t, v_0, a=0.5):
	# t_0
	t_0 = math.floor(d_list[0]/v_0)
	
	# update d
	d_list[1] += d_list[0]-(t_0*v_0) 

	# update t
	t_tmp = t - t_0
	if d_list[1] < (1/2)*v_0*t_tmp: # case 1: d < (1/2) * v_0*t => down
		M = d_list[1]/v_0 - 1 # M = t_1 + t_2 /2
		t_1 = math.floor(2*M*(1-1/(1+a)))
		t_2 = math.floor(2*M - 2*t_1)
		
		j_m = v_0 / (t_1 * (t_1 + t_2))
		v_e = 0 
		d_m = v_0*(t_1 + t_2 /2 + t_0 + 1)

		t_j = [[t_0, 0], [t_1, -j_m], [t_2, 0], [t_1, j_m]]
		return {"tick_jerk": t_j, "v_e": v_e, "d": d_m, "t_e": t_0+2*t_1+t_2, "t_init": t}


	elif 1.2 * d_list[1] > v_0 * t_tmp: # case 2: t is short flat down
		t_tmp = math.floor(1.2*d_list[1] / v_0)
	
	# M = t_1 + t_2/2
	M = t_tmp - (d_list[1]/v_0) +1
	t_4 = math.floor(d_list[1]/v_0 - 1 - M )
	t_1 = math.floor(2*M*(1-1/(1+a)))
	t_2 = math.floor(2*M - 2*t_1)

	j_m = v_0 / (t_1 * (t_1 + t_2))
	v_e = 0 
	d_m = v_0*(t_1 + t_2 /2 + t_4 + 1)+t_0*v_0

	t_j = [[t_0, 0], [t_4, 0], [t_1, -j_m], [t_2, 0], [t_1, j_m]]
	return {"tick_jerk": t_j, "v_e": v_e, "d": d_m, "t_e": t_0+2*t_1+t_2+t_4, "t_init": t}

		
def type_3_j_a_v(j, a, d, v_0):
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
given the number of ticks and initial values, find final a, v, q
a[t-1] = a_0 + t*j0

"""
def tavq(t, j0, a0 = 0, v0 = 0, q0 = 0):
	q0 += t*v0 + t*(t-1)*a0/2 + t*(t-1)*(t-2)*j0/6
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
"""

def plot(a, v, t, c = "r"):
	plt.figure(1)
	
	plt.subplot(311)
	plt.plot(a, c+'-')
	plt.title("a_end "+ str(a[-1]) )

	plt.subplot(312)
	plt.plot(v, c+'-')
	plt.title("v_end " + str(v[-1]))


	plt.subplot(313)
	plt.plot(t, c+'-')
	plt.title("t_end "+ str(t[-1]))
	plt.show()


def main_1():
	d_list = [[0, 1e5, 0]]
	t = 1e5
	a_avg = 0
	v_0 = [1]
	v_e = [1]
	j = 0
	print("prm: ", "d_list: ", d_list, " t: ",t, " a_avg: ", a_avg, " v_0: ", v_0, " v_e: ", v_e, " j: ", j)

	tick_jerk = []
	for i in range(len(d_list)): 
		result = profile(d_list[i], t, a_avg, v_0[i], v_e[i], j) #profile(d_list, t, a_avg, v_0, v_e, a)
		tick_jerk += result["tick_jerk"]

	jerk = [x[1] for x in tick_jerk]
	ticks = [x[0] for x in tick_jerk]
	a, v, q = tick(jerk, ticks, v_0=v_0[0])
	plot(a, v, q)

if __name__ == '__main__':
	import matplotlib.pyplot as plt
	main_1()
"""