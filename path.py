import numpy as np
import matplotlib.pyplot as plt
from profile import *
"""
prm: 
	p_1: start
	p_2: end
return: [A, B, d] 
p(t) = A +q(t)*B 
"""
def line_prm(p_1, p_2):
	d = np.linalg.norm(p_2-p_1)
	a = p_1
	if d == 0:
		b = (p_2-p_1)
	else:
		b =  (p_2-p_1)/d		
	return [{"type": "line", "A": a, "B": b, "d": d}]

def line_q(a, b, q):
	return a+q*b 

"""
prm: dictionary

"""
def line_tick(prm):
	# find line prm
	# [{"type": "line", "A": a, "B": b, "d": d}]
	line_prm = line_prm(prm["p_1"], prm["p_2"])[0]

	# find ticks for the given d
	if prm["type"] == "type_1":
		# {"tick_jerk": [{"t", "j"},...], "v_e": v_0, "d": d_m}
		result = type_1(prm["j"], prm["a"], prm["v"], line_prm["d"], prm["v_0"] )
	elif prm["type"] == "type_2":
		result = type_2(prm["j"], prm["a"], prm["v"], line_prm["d"], prm["v_0"] )
	elif prm["type"] == "type_3":
		result = type_3(prm["j"], prm["a"], line_prm["d"], prm["v_0"] )
	else:
		return False

	data = []
	for tj in result["tick_jerk"]:
		



"""
prm:
	p_1: center
	p_2: start
	p_3: vector perpendicular to p_2-p_1 and length r: sets the move direction and circle plane  
	d: distance
return: [A, B, C]
p(t) = A +B*cons(theta)+C*sin(theta)  
"""
def circle_prm(p_1, p_2, p_3, d, r_inv):
	return [{"type": "circle", "prm":[p_1, p_2-p_1, p_3], "distance": d, "r_inv": r_inv}]	

def circle_q(a, b, c,r_inv, q):
	theta = q*r_inv
	return a + b*np.cos(theta)+c*np.sin(theta)


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