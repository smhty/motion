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

import sys
sys.setrecursionlimit(1500)


"""
continous
"""
class profile(object):
	"""docstring for profile"""
	def __init__(self, cont = 1):
		super(profile, self).__init__()
		self.depth = 0

	def q_n(self, t, j, a0, v0, q0 = 0):
		return q0 + t*v0 + t*(t-1)/2 * a0 + t*(t-1)*(t-2)/6 * j  

	def q_n_prime(self,t,j,a0,v0):
		return v0 - a0 / 2 + t * a0 + j * (3*t*t - 6*t + 2) /6  #(t-1)*(t-2)/6 * j + t*(t-2)/6 * j  + t*(t-1)/6 * j 

	def v_n(self, t, j, a0, v0):
		return v0 + t*a0 + t*(t-1)/2 * j  

	def a_n(self, t, j, a0):
		return a0 + t*j  

	def sign(self, x):
		if x >= 0:
			return 1
		return -1
	
	# ax^2 + bx +c =0
	def solver_quadratic(self, a, b, c):
		delta = b**2 - 4*a*c
		if delta < 0:
			return False
		delta = math.sqrt(delta)
		return [(-b-delta)/(2*a), (-b+delta)/(2*a)]  	


	"""
	given t, j, a0, v0, d
	we travel the path with t and j and return the remaining
	"""
	def slice(self,t, j0, a0, v0, d):
		print("xxx12: ", t, j0, a0, v0, d)
		q_t = self.q_n(t, j0 , a0 , v0 )
		print("xxx13: ", q_t)
		if q_t <= d:
			print("xxx18: ", q_t, )
			return {"tick_jerk": [[t,j0]], "ae": self.a_n(t,j0,a0) , "ve": self.v_n(t,j0,a0,v0) , "de": q_t, "finish": 1}
		
		if j0 ==0:
			result = self.solver_quadratic(a0/2, v0 - a0/2, -d)
			print("xxx81: ", result)
			if result:
				x = math.floor(min(result))
				if x >= 0:
					return {"tick_jerk": [[x,j0]], "ae": self.a_n(x,j0,a0), "ve": self.v_n(x,j0,a0,v0) , "de": self.q_n(x,j0,a0,v0), "finish": 0}
			else:
				return {"tick_jerk": [[0,j0]], "ae": self.a_n(0,j0,a0), "ve": self.v_n(0,j0,a0,v0) , "de": self.q_n(0,j0,a0,v0), "finish": 0}
		

		#using newton method:
		x_u = t
		q_u = q_t

		x_l = 0
		q_l = 0

		# find the right initial point
		if j0 > 0:
			x = max(0, min(((6*d)/j0)**(1/3), t))
		else:
			x = t
		
		print("xxx14: ", x, t)
		q_n = self.q_n(x , j0 , a0 , v0 , 0)

		i = 0
		while i < 100 and math.floor(x_u) > math.ceil(x_l):
			x = x - (q_n - d) / self.q_n_prime(x , j0 , a0 , v0);
			# update q_n
			q_n = self.q_n(x , j0 , a0 , v0 , 0)

			print("xxx15: ", x, q_n)
			if x >= x_l and x <= x_u:
				# update x_u, x_l
				if q_n > d and q_n < q_u:
					x_u = x
					q_u = q_n
				elif q_n < d and q_n > q_l:
					x_l = x
					q_l = q_n
				else:
					print("xxx")
					x_u = x
					q_u = q_n
					x_l = x
					q_l = q_n

			print("xxx16: ", x, q_n, q_u, q_l)
			i = i+1

		# integer solution
		x = math.floor(x_u)		
		print("xxx17: ", x, i)
		return {"tick_jerk": [[x,j0]], "ae": self.a_n(x,j0,a0), "ve": self.v_n(x,j0,a0,v0) , "de": self.q_n(x,j0,a0,v0), "finish": 0}

	"""
	assumption
	jm, am, vm are always positive
	a0 <= v0
	j <= jm
	??? 
	a0 < 0 choose right t s.t v is not negative
	"""
	def cont(self, a0, v0, jm, am, vm, d):
		print("xxx1: ", a0, v0, jm, vm,  d)
		rtn = {
			"tick_jerk": [],
			"ae": a0,
			"ve": v0,
			"de": 0
		}
		# check that d is valid
		if d < v0:
			return rtn
		  
		if vm > 0: # middle
			
			if a0 == 0: # initial accel is zero

				"""
				v0 == vm
				"""
				if v0 == vm: # stick to this velocity
					t4 = math.floor(d/vm)
					rtn["tick_jerk"] = [[t4, 0]]
					rtn["de"] = self.q_n(t4, 0, 0, v0)
					print("xxx8: ", a0, v0, jm, d)
					return rtn 
				
				"""
				v0 != vm
				"""
				# find t1_
				t1 = math.floor(am/jm) +1 
				if t1 * am >= abs(vm - v0):
					# t1**2 * jm = vm- v0
					t1 = math.floor(math.sqrt(abs(vm-v0)/jm)) + 1
					j = (vm - v0) / (t1**2) 
					a = j * t1
					t2 = 0
					print("xxx9: ", a0, v0, jm, d)
				else:
					t2 = math.floor(abs(vm - v0) / am) - t1 + 1
					a = (vm- v0)/ (t1 + t2)
					j = a / t1
					print("xxx10: ", t1,  t2, a, j, vm, v0)

				# 3 motions: [t1,sign* jm], [t2, 0], [t1, -sign*jm]
				rtn = {
					"tick_jerk": [],
					"ae": a0,
					"ve": v0,
					"de": 0
				}
				print("xxx20: ", t1, rtn)
				slc= self.slice(t1, j , a0, v0, d)
				print("xxx11: ", slc)
				rtn["tick_jerk"] += list(slc["tick_jerk"])
				rtn["ae"] = slc["ae"]
				rtn["ve"] = slc["ve"]
				rtn["de"] += slc["de"] 
				if slc["finish"]:
					print("xxx21: ", t1)
					slc = self.slice(t2, 0, a, rtn["ve"] , d - rtn["de"])
					rtn["tick_jerk"] += list(slc["tick_jerk"])
					rtn["ae"] = slc["ae"]
					rtn["ve"] = slc["ve"]
					rtn["de"] += slc["de"] 					
					if slc["finish"]:
						print("xxx22: ", t1)
						slc = self.slice(t1, -j, a, rtn["ve"] , d - rtn["de"])
						rtn["tick_jerk"] += list(slc["tick_jerk"])
						rtn["ae"] = slc["ae"]
						rtn["ve"] = slc["ve"]
						rtn["de"] += slc["de"] 					
						if slc["finish"]:
							print("xxx23: ", t1)
							slc = self.cont(0, rtn["ve"], jm, am, vm, d - rtn["de"])
							rtn["tick_jerk"] += list(slc["tick_jerk"])
							rtn["ae"] = slc["ae"]
							rtn["ve"] = slc["ve"]
							rtn["de"] += slc["de"] 					

				return rtn						 
			# a0 != 0
			# find t1 and t3
			result = self.solver_quadratic(self.sign(a0)*jm, 2*a0, v0 - a0**2/(2*jm) + a0/2 - vm )
			if result:
				if a0 >= 0:
					t1 = min(result)
				else:
					t1 = max(result)
				#t1 = min(result)
				if t1 >= 1:
					t3 = t1 + abs(a0)/jm
					t2 = 0
					
					# adjust t1, t3
					a1 = self.a_n(t1, self.sign(a0)*jm, a0)
					v1 = self.v_n(t1, self.sign(a0)*jm, a0, v0)


					# go toward am: j = self.sign(a0)* self.sign(am - abs(a0))* jm
					t1_adj = math.floor(abs((min(abs(a1), am) - abs(a0))) / jm)
					v1_adj = self.v_n(t1_adj, self.sign(a0)* self.sign(am - abs(a0))* jm, a0, v0)
					a1_adj = self.a_n(t1_adj, self.sign(a0)* self.sign(am - abs(a0))* jm, a0)

					# zero jerk: j = 0
					t_tmp = abs(a1 - a1_adj) / jm				
					v_tmp = self.v_n(t_tmp, - self.sign(a0)*jm, a1, v1)
					t2_adj = math.floor((v_tmp - v1_adj) / (self.sign(a0)*am)) 
					v2_adj = self.v_n(t2_adj, 0, a1_adj, v1_adj)

					# go toward a = 0: j = - self.sign(a0)*jm  => adjust j3
					t3_adj = math.floor(abs(a1_adj) / jm)
					j3 = -a1_adj / t3_adj

					# 3 motions: [t1,sign* jm], [t2, 0], [t1, -sign*jm]				
				


					print("xxx24: ", t1_adj)
					slc = self.slice(t1_adj, self.sign(a0)* self.sign(am - abs(a0))* jm , a0, v0, d)
					rtn["tick_jerk"] += list(slc["tick_jerk"])
					rtn["ae"] = slc["ae"]
					rtn["ve"] = slc["ve"]
					rtn["de"] += slc["de"] 				
					if slc["finish"]:
						print("xxx25: ", t2_adj, rtn)
						slc = self.slice(t2_adj, 0, a1_adj, v1_adj , d - rtn["de"])
						rtn["tick_jerk"] += list(slc["tick_jerk"])
						rtn["ae"] = slc["ae"]
						rtn["ve"] = slc["ve"]
						rtn["de"] += slc["de"]
						if slc["finish"]:
							print("xxx26: ", t3_adj, rtn)
							slc = self.slice(t3_adj, j3, a1_adj, v2_adj , d - rtn["de"])
							rtn["tick_jerk"] += list(slc["tick_jerk"])
							rtn["ae"] = slc["ae"]
							rtn["ve"] = slc["ve"]
							rtn["de"] += slc["de"]
							if slc["finish"]:
								slc = self.cont(0, rtn ["ve"], jm, am, vm, d - rtn["de"])
								rtn["tick_jerk"] += list(slc["tick_jerk"])
								rtn["ae"] = slc["ae"]
								rtn["ve"] = slc["ve"]
								rtn["de"] += slc["de"]

					return rtn						 

			# go toward a= 0 as fast as possible
			if a0 < -v0:

				a0 = -v0			
			t3 = math.floor(abs(a0)/jm) + 1
			j3 = -a0 / t3
			v3 = self.v_n(t3, j3, a0, v0)  
			if v3 < 0:
				t3 = math.floor(2*v0/abs(a0)) - 1
				j3 = -a0 / t3

			print("xxx27: ", t3)
			slc = self.slice(t3, j3 , a0, v0, d)
			print("xxx37: ", slc)
			rtn["tick_jerk"] += list(slc["tick_jerk"])
			rtn["ae"] = slc["ae"]
			rtn["ve"] = slc["ve"]
			rtn["de"] += slc["de"]
			print("xxx2: ", a0, v0, jm, d, slc["ve"])
			if slc["finish"]:
				print("xxx3: ", a0, v0, jm, d, rtn ["ve"])
				slc = self.cont(0, rtn ["ve"], jm, am, vm, d - rtn["de"])
				print("xxx7: ", slc["de"])
				rtn["tick_jerk"] += list(slc["tick_jerk"])
				rtn["ae"] = slc["ae"]
				rtn["ve"] = slc["ve"]
				rtn["de"] += slc["de"]
			return rtn						 

		else: # stop
			"""
			a0 == 0
			"""
			print("xxx30: ", a0, v0, jm, vm,  d)
			if a0 == 0:
				# sadegh add
				t1_a = math.floor(am / jm)
				t1_v = math.floor((v0/jm)**(1/2))
				t1_d = math.floor((d/v0)-1)

				# find t_1 and set t_2
				t1 = min(t1_a, t1_v, t1_d)
				t1 = max(0, t1)
				t2 = 0

				if t1 == 0:
					if v0 > 0: 
						t4 = math.ceil(d/v0) - 1
						t4 = max(0, t4)
					else: 
						t4 = 0 
					j_t = 0
				else:
					if t1_a <= min(t1_v, t1_d):
						t2_v = math.floor((v0 / (jm*t1)) - t1)
						t2_d = 2*(t1_d - t1)

						# t2
						t2 = min(t2_v, t2_d)
						t2 = max(0, t2)
					
					# t4
					t4 = math.ceil((d/v0) - t1 - t2/2 - 1)
					t4 = max(0, t4)
					j_t = v0 / (t1 * (t1 + t2))
				#j_m
				rtn["tick_jerk"] += [[t4, 0], [t1, -j_t], [t2, 0], [t1, j_t]]
				rtn["ae"] = 0
				rtn["ve"] = 0
				rtn["de"] = j_t* t1 * (t1 + t2)* (t1 + t2 /2 + t4 + 1) 

				return rtn
				# sadegh add
				# sadegh remove
				"""
				#v0*(t_1 + t_2 /2 + t_4 + 1) = d
				t1 = math.floor(d/v0) - 1
				j = v0 / (t1**2)
				rtn["tick_jerk"] += [[t1, -j], [t1, j]]
				rtn["ae"] = 0
				rtn["ve"] = 0
				rtn["de"] = v0 * (t1+1)	
				print("xxx4: ", a0, v0, jm, d)					
				return rtn
				"""
				# sadegh remove


			"""
			d1 = v0*t + a0/3 * (t**2 - 1)
			v1 = v0 + a0 (t+1) /2
			a0 = j*t1
			"""
			d1 = d/2
			result = self.solver_quadratic(a0/3, v0, -a0/3 - d1)
			if result:
				t1 = math.floor(min(result))
				if t1 < 0:
					t1 = math.floor(max(result))
				
				print("xxx34: ", result)
				# adjust t1 for negative a0
				if a0 < 0:
					if t1 > -2*(v0/a0) - 1:
						t1 = math.floor(-(v0/a0) - 0.5)
						print("xxx32: ", t1)


				j1 = -a0/t1
				v1 = self.v_n(t1, j1, a0, v0)
				d1 = self.q_n(t1, j1, a0, v0)

				rtn["tick_jerk"] += list([[t1, j1]])
				print("xxx33: ", d, d1)
				slc = self.cont(0, v1, jm, am, vm, d - d1)
				rtn["tick_jerk"] += list(slc["tick_jerk"])
				rtn["ae"] = slc["ae"]
				rtn["ve"] = slc["ve"]
				rtn["de"] += slc["de"]
				print("xxx5: ", a0, v0, jm, d)
				return rtn

		print("xxx6: ", a0, v0, jm, d)
		return rtn				


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

def profile_2_t(j, a, t, d, v=1, v_0 = 0):
	return 0
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
solve the profile based on time t
solve the curve based on a maximum velocity v
if t_m >= t 
	then we already have a solution
else
	then t_4 += (t-t_m) 
	update j_m based on the d and the new t_4
"""
def profile_t(j, a, t, d, _type=1, v=1, v_0 = 0):
	return True

"""
factor = t2 / t1
"""
def profile_no_jerk(a, vm, d, factor = 1):
	# find n1 = 2*t1 + t2
	n1 = math.floor(vm/a)
	# adjust vm and a
	if vm * n1 > d:
		vm = d / n1
	# find t4
	t4 = math.floor(d/vm) - n1	
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
given the number of ticks and initial values, find final a, v, q

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
	plt.plot(a, c+'o')
	plt.title("a_end "+ str(a[-1]) )
	"""
	#plt.title("Accelaration" )
	#plt.xticks([])
	#plt.yticks([])
	"""

	plt.subplot(312)
	plt.plot(v, c+'o')
	plt.title("v_end " + str(v[-1]))
	"""
	plt.title("Velocity" )
	plt.xticks([])
	plt.yticks([])
	"""


	plt.subplot(313)
	plt.plot(t, c+'o')
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


def main_1():
	j = 1.0e-10
	a = 1.0e-6
	v = 0.08
	v_0 =0
	d = 10000
	t = 10000000000

	"""
	x = time.time()
	p = profile()
	slc = p.slice(t, j, a, v_0, d)
	print(time.time() - x)

	print(p.q_n(slc["tick_jerk"][0][0], slc["tick_jerk"][0][1], a, v_0), d)
	"""
	factor = 1

	#result = profile_no_jerk(a, v, d, factor)
	result = profile_1(j, a, v, d, v_0)
	jerk = [x[1] for x in result["tick_jerk"]]
	ticks = [x[0] for x in result["tick_jerk"]]
	a, v, q = tick(jerk, ticks, v_0)
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
	main_cont()
	#main_test()
	