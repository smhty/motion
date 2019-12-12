import numpy as np
import matplotlib.pyplot as plt
from profile import *
"""
prm: 
	p_1: start
	p_2: end
return: A, B 
p(t) = A +q(t)*B 
"""
def path_line_prm(p_1, p_2, d):
	A = p_1
	if d == 0:
		B = (p_2-p_1)
	else:
		B =  (p_2-p_1)/d		
	return {"A": A, "B": B}

def path_line_q(A, B, q):
	return A+q*B 

"""
prm: dictionary
	profile
	path
	"A"
	"B"
	...
	"p_1", "p_2", "d", "r_inv"
	"j", "a", "v", "d", "v_0"	
	"id"
"""
def path_tick(prm):
	# path
	if prm["path"] == "line":
		# {"A", "B"}
		path = path_line_prm(prm["p_1"], prm["p_2"], prm["d"]) 

	"""
	profile
	find ticks for the given d
	{"tick_jerk": [{"t", "j"},...], "v_e": v_0, "d_m": d_m}
	"""
	if prm["profile"] == "profile_1":
		profile = profile_1(prm["j"], prm["a"], prm["v"], prm["d"], prm["v_0"] )
	elif prm["profile"] == "profile_2":
		profile = profile_2(prm["j"], prm["a"], prm["v"], line_prm["d"], prm["v_0"] )
	else prm["profile"] == "profile_3":
		profile = profile_3(prm["j"], prm["a"], line_prm["d"], prm["v_0"] )


	# segment order
	profile["tick_jerk"] = tick_jerk_segment_order(profile["tick_jerk"])
	
	# profile to tick_jerk
	for t_j in profile["tick_jerk"]:
		tj["id"] = prm["id"]
		tj["path"] = prm["path"]
		tj["A"] = prm["A"]
		tj["B"] = prm["B"]

	return profile

"""
prm:
	p_1: start
	C: center
	V: vector perpendicular to p_1-C and length r: sets the move direction and circle plane  
return: [A, B, C]
p(t) = A +B*cons(theta)+C*sin(theta)  
"""
def circle_prm(p_1, C, r):
	if p_1 == C:
		return 
	return [{"type": "circle", "A": C, "B": p_1-C, "C": V, "r_inv": r_inv}]	

"""
prm:
	p_1: center
	p_2: start
	p_3: vector perpendicular to p_2-p_1 and length r: sets the move direction and circle plane  
	d: distance
return: [A, B, C]
p(t) = A +B*cons(theta)+C*sin(theta)  
"""
def circle_prm_copy(p_1, p_2, p_3, d, r_inv):
	return [{"type": "circle", "A": p_1, "B": p_2-p_1, "C": p_3, "d": d, "r_inv": r_inv}]	



def circle_q(a, b, c,r_inv, q):
	theta = q*r_inv
	return a + b*np.cos(theta)+c*np.sin(theta)

"""
prm: dictionary
	profile
	path
	"A"
	"B"
	...
	"p_1", "p_2", "d", "r_inv"
	"j", "a", "v", "d", "v_0"	


"""
def circle_tick(prm):
	# find line prm
	# [{"type": "circl", "A": a, "B": b, "d": d}]
	if prm["path"] == "circle":
		path = circle_prm("p_1", "p_2", "p_3", d, r_inv)[0]

	# find ticks for the given d
	if prm["profile"] == "profile_1":
		# {"tick_jerk": [{"t", "j"},...], "v_e": v_0, "d": d_m}
		profile = profile_1(prm["j"], prm["a"], prm["v"], prm["d"], prm["v_0"] )
	elif prm["profile"] == "profile_2":
		profile = profile_2(prm["j"], prm["a"], prm["v"], line_prm["d"], prm["v_0"] )
	elif prm["profile"] == "profile_3":
		profile = profile_3(prm["j"], prm["a"], line_prm["d"], prm["v_0"] )
	else:
		return False

	# add segment order

	result["tick_jerk"] = tick_jerk_segment_order(result["tick_jerk"])
	for tj in result["tick_jerk"]:
		tj["id"] = prm["id"]
		tj["type"] = prm["type"]
		tj["A"] = prm["A"]
		tj["B"] = prm["B"]
	
	return result


"""
prm:
	p_1: start
	p_2: middle
	p_3: end 
	r: corner radius
return:
	[line, circle, line]
"""
def line_circle_line_prm(p_1, p_2, p_3, r):
	# vector
	v_1 = np.subtract(p_1, p_2)
	v_2 = np.subtract(p_3, p_2)

	# norm
	v_1_n = np.linalg.norm(v_1)
	v_2_n = np.linalg.norm(v_2)

	if v_1_n == 0 or v_2_n == 0:
		return path_line(p_1, p_3)

	# the angle forms by p_1, p_2, p_3
	theta = np.arccos(np.inner(v_1, v_2)/(v_1_n * v_2_n))
	
	# check for 180 degree
	cos_t_2 = np.cos(theta/2)
	if cos_t_2 == 0:
		return path_line(p_1, p_3)  

	# tan theta/2
	tan_t_2 = np.sin(theta/2) / cos_t_2

	# modify r
	r = min(r, v_1_n*tan_t_2*0.9, v_2_n*tan_t_2*0.9)

	# point where v_1 touches the circle
	A = p_2 + r/tan_t_2 * v_1/v_1_n 
	
	# point where v_2 touches the circle
	B = p_2 + r/tan_t_2 * v_2/v_2_n 	

	# center of the circle
	C = p_2 + (r/np.sin(theta)) *(v_1/v_1_n + v_2/v_2_n)
	
	return line_prm(p_1, A) +  circle_prm(C, A, - r*v_1/v_1_n, r*(np.pi - theta), 1/r) + line_prm(B, p_3)  


if __name__ == '__main__':
	
	p_3 = [10,0]
	p_2 = [5,4]
	p_1 = [22,12]
	r = 30

	path = line_circle_line_prm(p_1, p_2, p_3, r)

	# Create plots with pre-defined labels.
	fig, ax = plt.subplots()	
	
	for p in path:
		x = []
		y = []
		if p["type"] == "circle":	
			t = 0
			while t <= p["distance"]:
				point = circle_q(p["prm"][0], p["prm"][1], p["prm"][2], p["r_inv"], t)
				x.append(point[0])
				y.append(point[1])
				t += 1
			print(p["prm"][0], p["prm"][1], p["prm"][2], p["r_inv"])		 
		elif p["type"] == "line":	
			t = 0
			while t <= p["distance"]:
				point = line_q(p["prm"][0], p["prm"][1], t)
				x.append(point[0])
				y.append(point[1])
				t += 1			

		ax.plot(x, y, "o-")
	plt.show()	