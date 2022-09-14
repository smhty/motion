import numpy as np
import math
import matplotlib.pyplot as plt
import profile
import kinematic
import json
import traverse




def joint(jnt_s, jnt_e, j, a, t, profile_type=profile_1_t):
	jnt_s = np.array(jnt_s)
	jnt_e = np.array(jnt_e)

	# traverse
	trv = traverse.linear(jnt_s, jnt_e)

	# find profile for the path distance
	prf = profile.profile_1(j, a, v, path.d, v_0)
	
	# tick jerk
	tick_jerk_l = profile.tick_jerk(result["tick_jerk"])

	# sampling path
	return path_sampling(tick_jerk_l, path, False)	


def joint_copy(jnt_s, jnt_e, j, a, v, v_0 = 0):
	jnt_s = np.array(jnt_s)
	jnt_e = np.array(jnt_e)

	#path
	path = linear(jnt_s, jnt_e)

	# find profile for the path distance
	result = profile.profile_1(j, a, v, path.d, v_0)
	
	# tick jerk
	tick_jerk_l = profile.tick_jerk(result["tick_jerk"])

	# sampling path
	return path_sampling(tick_jerk_l, path, False)	
		

def circle(xyz_s, xyz_m, xyz_e, j, a, v, v_0 = 0, turn = 0):
	xyz_s = np.array(xyz_s)
	xyz_e = np.array(xyz_e)

	#path
	path = circular(xyz_s, xyz_m, xyz_e, turn)

	print("length: ", path.d)
	# find profile for the path distance
	result = profile.profile_1(j, a, v, path.d, v_0)
	
	# tick jerk
	tick_jerk_l = profile.tick_jerk(result["tick_jerk"])

	# sampling path
	return path_sampling(tick_jerk_l, path, True)


def line(xyz_s, xyz_e, j, a, v, v_0 = 0):
	xyz_s = np.array(xyz_s)
	xyz_e = np.array(xyz_e)

	#path
	path = linear(xyz_s, xyz_e)

	# find profile for the path distance
	result = profile.profile_1(j, a, v, path.d, v_0)
	
	# tick jerk
	tick_jerk_l = profile.tick_jerk(result["tick_jerk"])

	# sampling path
	return path_sampling(tick_jerk_l, path, True)


# start from line_s, go to line_e, after that go on a circle touches circle_m, and stop at circle_e
def line_circle(line_s, line_e, circle_m, circle_e, j, a, v, v_0 = 0, turn = 0):
	return True

def point_to_motor(point, invk):
	point = np.array(point)
	if invk: # run inverse kinematic
		k = kinematic.kinematic()
		point = np.array(k.inverse(point)[0])
	m = kinematic.motor()
	return np.array(m.joint_to_motor(point))	


def motor_cmd(motor_s, motor_e, t, njerk, naccel, nvel):
	motor_path = linear(motor_s, motor_e)

	return {
			"t": t,
			"jerk": njerk* motor_path.d  ,
			"accel": naccel * motor_path.d,
			"vel": nvel * motor_path.d,
			"q": 0,  
			"a0": motor_path.a[0],
			"a1": motor_path.a[1],
			"a2": motor_path.a[2],
			"a3": motor_path.a[3],
			"a4": motor_path.a[4],
			"a5": 0, 
			"a6": 0,
			"a7": 0,   
			"b0": motor_path.b[0],
			"b1": motor_path.b[1],
			"b2": motor_path.b[2],
			"b3": motor_path.b[3],
			"b4": motor_path.b[4],
			"b5": 0, 
			"b6": 0,
			"b7": 0,
			"m0": motor_e[0],
			"m1": motor_e[1],
			"m2": motor_e[2],
			"m3": motor_e[3],
			"m4": motor_e[4],
			"m5": 0,
			"m6": 0,
			"m7": 0,
			"njerk": njerk, 
			"naccel": naccel,
			"nvel": nvel,			
			}


def path_sampling(tick_jerk_l, path, invk, a_0 = 0, v_0 = 0, d_0 = 0, t_s = 250):
	motor = []

	# initialization
	motor_e = point_to_motor(path.traverse(d_0), invk)
	
	# segments
	for tick_jerk in tick_jerk_l:
		t_s_total = 0
		l = math.floor(tick_jerk["tick"]/t_s)
		for j in range(l):				
			cmd = {}
			
			# tick
			if j < l-1:
				t = t_s
			else: # j == l-1
				t = tick_jerk["tick"] - j*t_s
			
			# initial jerk
			j_0 = tick_jerk["jerk"]

			# final d, v, a
			t_t_1 = t*(t-1)
			t_t_1_t_2 = t_t_1*(t-2)
			d = t * v_0 + t_t_1 * a_0/2 + t_t_1_t_2 * j_0/6 
			
			# normalized 
			njerk = j_0 / d    
			naccel = a_0 / d
			nvel = v_0 / d

			d_0 += d
			v_0 += t * a_0 + t_t_1 * j_0/2 
			a_0 += t * j_0

			# motor start and end
			motor_s = np.array(motor_e)
			motor_e = point_to_motor(path.traverse(d_0), invk)
		
			# motor parameter
			motor.append(motor_cmd(motor_s, motor_e, t, njerk, naccel, nvel))				

	return motor	


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
	elif prm["profile"] == "profile_3":
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

def main_copy():
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

def write_to_file(cmds, path = "cmds.txt"): 
	f = open(path, "w")
	for cmd in cmds:
		f.write(json.dumps(cmd) + "\n")
	f.close()

def main_line():

	j = 1.0e-12
	a = 1.0e-7
	v = 0.001
	v_0 = 0

	jnts = [
			[0, 0, 0, 0, 0],
			[0, 50, 0, 0, 0],
			[0, 0, 0, 0, 0],
			[0, 50, 0, 0, 0],
			[0, 0, 0, 0, 0],]

	k = kinematic.kinematic()
	result = []
	for i in range(len(jnts) - 1):
		result += line(k.forward(jnts[i]) , k.forward(jnts[i+1]) , j, a, v, v_0)

	# write to file
	write_to_file(result)


def main_joint():

	j = 1.0e-12
	a = 1.0e-7
	v = 0.001
	v_0 = 0

	jnts = [
			[0, 0, 0, 0, 0],
			[0, 50, 0, 0, 0],
			]

	result = []
	for i in range(len(jnts) - 1):
		result += joint(jnts[i],jnts[i+1] , j, a, v, v_0)

	# write to file
	write_to_file(result)	


def main_circle():
	j = 1.0e-11
	a = 1.0e-6
	v = 0.002
	v_0 = 0
	turn = 0
	joint_0 = [0, 44, -111, 67, 0]

	k = kinematic.kinematic()
	xyz_s = k.forward(joint_0) 

	xyz_m = np.array(xyz_s)
	xyz_m[0] += 100

	xyz_e = np.array(xyz_s)
	xyz_e[0] += 60
	xyz_e[1] += 60
	xyz_e[2] += 60


	xyz_s = [400.0, 0.0 , 180.0, 0.0, 0.0]
	xyz_m = [300.0, 0.0 , 180.0, 0.0, 0.0]
	xyz_e = [351.0, -50.0 , 180.0, 0.0, 0.0]
	# list_tvaj
	result = []
	result += circle(xyz_s, xyz_m, xyz_e, j, a, v, v_0, turn)
	#result += circle(xyz_e, xyz_m, xyz_s, j, a, v, v_0, turn)

	write_to_file(result)

	"""
	# 3d plot
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	# For each set of style and range settings, plot n random points in the box
	# defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
	for r in result:
		xs = r["b0"]
		ys = r["b1"]
		zs = r["b2"]
		ax.scatter(xs, ys, zs, c="b", marker="o")
	plt.show()
	"""


def round_line(p0, p1, p2, r):
	return True


def main_line_circle():
	#give two lines, and round corner

	j = 1.0e-12
	a = 1.0e-7
	v = 0.001
	v_0 = 0
	r 

	jnts = [
			[0, 0, 0, 0, 0],
			[0, 45, -135, 90, 0],
			[0, 90, -90, 0, 0],
			]

	k = kinematic.kinematic()
	
	pnts = [k.forward(jnt) for jnt in jnts]
	#######
	# vector
	v_1 = pnts[0] - pnts[1]
	v_2 = pnts[2] - pnts[1]

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
	main_circle()