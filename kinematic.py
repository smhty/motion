import math
import profile
import numpy as np

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import json
class motor(object):
	"""docstring for motor"""
	def __init__(self, micro_step = 4000):
		super(motor, self).__init__()
		self.micro_step = micro_step

	def joint_to_motor(self, joint):
		return [
			-joint[0] * 16.0 * self.micro_step / 360.0,
			-joint[1] * 20.0 * self.micro_step / 360.0,
			-joint[2] * 20.0 * self.micro_step / 360.0,
			((-joint[3] + joint[4]) / 360.0) * (self.micro_step * 4.0),
			((joint[3] + joint[4]) / 360.0) * (self.micro_step * 4.0),
			0.0,
			0.0,
			0.0,			
		]

	def motor_to_joint(self, motor):
		return [
			-motor[0] *360.0/( 16.0 * self.micro_step),
			-motor[1] *360.0/( 20.0 * self.micro_step),
			-motor[2] *360.0/( 20.0 * self.micro_step),
			(-motor[3] + motor[4]) * 360.0 / (2.0 * self.micro_step * 4.0),
			(motor[3] + motor[4]) * 360.0 / (2.0 * self.micro_step * 4.0),
			0.0,
			0.0,
			0.0,			
		]
			
"""
use kinematic.jpg 
all dimensions are in mm
l = [
		[l_00 = 95.53023, 0, l_02 = 200.09099],
		[l_10 = 203.2, 0, 0],
		[l_20 = 152.4, 0, 0],
		[0, 0, 0],
		[0,  0, l_42 = 48.92446]
	]
l = [
		[95.53023, 0, 200.09099],
		[203.2, 0, 0],
		[152.4, 0, 0],
		[0, 0, 0],
		[0,  0, 48.92446]
	]	
"""
class kinematic(object):
	"""docstring for kinematic"""
	def __init__(self, 
		l = [
		[95.53023, 0, 200.09099],
		[203.2, 0, 0],
		[152.4, 0, 0],
		[0, 0, 0],
		[0,  0, 48.92446]]
	):
		super(kinematic, self).__init__()
		# l matrix
		self.l = l
		self.e = 0.001
	

	"""
	forward kinematics: joint to xyz
	"""
	def forward(self, joint):
		# joint to radian
		joint = [math.radians(j) for j in joint]

		# first we find x, y, z assuming base rotation is zero (j_0 = 0). Then we rotate everything
		# then we rotate the robot around z axis for j_0
		a = math.degrees(joint[1] + joint[2] + joint[3])
		b = math.degrees(joint[4])
		tmp_d = self.l[0][0] + self.l[1][0] * math.cos(joint[1]) + self.l[2][0] * math.cos(joint[1] + joint[2]) + self.l[4][2] * math.cos(joint[1] + joint[2] + joint[3])
		x = tmp_d * math.cos(joint[0])
		y = tmp_d * math.sin(joint[0])
		z = self.l[0][2] + self.l[1][0] * math.sin(joint[1]) + self.l[2][0] * math.sin(joint[1] + joint[2]) + self.l[4][2] * math.sin(joint[1] + joint[2] + joint[3])

		return [x, y, z, a, b]


	"""
	inverse kinematics: xyz to joint
	"""
	def inverse(self, xyz):
		x = xyz[0]
		y = xyz[1]
		z = xyz[2]
		a = math.radians(xyz[3])
		b = math.radians(xyz[4])
		
		joint = []	
		
		try:
			# first we find the base rotation
			theta_0 = math.atan2(y, x)


			# next we assume base is not rotated and everything lives in x-z plane
			x = math.sqrt(x ** 2 + y ** 2)

			# next we update x and z based on base dimensions and hand orientation
			x -= (self.l[0][0] + self.l[4][2] * math.cos(a))
			z -= (self.l[0][2] + self.l[4][2] * math.sin(a))

			# at this point x and z are the summation of two vectors one from lower arm and one from upper arm of lengths l1 and l2
			# let L be the length of the overall vector
			# we can calculate the angle between l1 , l2 and L
			L = math.sqrt(x ** 2 + z ** 2)
			
			if L >= self.l[1][0] + self.l[2][0]:
				theta_l1_L = 0  # l1 angle to L
				theta_l1_l2 = math.pi  # l1 angle to l2

			elif L <= self.l[1][0] - self.l[2][0]:
				theta_l1_L = 0  # l1 angle to L
				theta_l1_l2 = 0  # l1 angle to l2

			else:
				theta_l1_L = math.acos((self.l[1][0] ** 2 + L ** 2 - self.l[2][0] ** 2) / (2 * self.l[1][0] * L))  # l1 angle to L
				theta_l1_l2 = math.acos((self.l[1][0] ** 2 + self.l[2][0] ** 2 - L ** 2) / (2 * self.l[1][0] * self.l[2][0]))  # l1 angle to l2

			theta_L_x = math.atan2(z, x)  # L angle to x axis
			theta_1_0 = theta_L_x + theta_l1_L
			theta_1_1 = theta_L_x - theta_l1_L

			theta_2_0 = theta_l1_l2 - math.pi
			theta_2_1 = -(theta_l1_l2 - math.pi)

			theta_3_0 = a - theta_1_0 - theta_2_0
			theta_3_1 = a - theta_1_1 - theta_2_1

			joint = [[theta_0, theta_1_0, theta_2_0, theta_3_0, b],[theta_0, theta_1_1, theta_2_1, theta_3_1, b]]
			joint = [[math.degrees(j) for j in joint[0]], [math.degrees(j) for j in joint[1]]]

		except:
			pass

		return joint


"""
given a line do the following:
- find its motion profile
- sample it and convert it to joint
- find v_0, a_0, j_0 associate to each profile 

The sampleling is done based on parameter T, and on the following points 
Each segment has T points, except the last segment
the last segment has more than or equal to T points
[0, T-1]
[T, 2T-1]
[2T, 3T-1]
...
[kT, T_end] -> last part
k = floor((number_of_ticks / T))

initailly v[-1] is zero

in each segment do the following:

main equation: 
q[n] =n * v + n*(n-1) * a/2 + n*(n-1)*(n-2)*j/6
v[n] =        v +     n * a +      n*(n-1)* j/2 	
a[n] =                        a +       n * j 	



q(n) = q(n-1)+v(n-1)
v(n) = v(n-1) + a(n-1)
a(n) = a(n-1) + j

a, v, j 
*v, a, j are the last parameter, *[T-1] from previous segment

we assume that v is continous from segment to segment, so v[-1] comes from

v = v[T-1] from the previous segment calculation
v[T-1] = distance from p[0] in next segment -  p[T-1] in current segment  

"""
class path(object):
	"""docstring for path"""
	def __init__(self):
		super(path, self).__init__()
		self.kinematic = kinematic()
		self.motor = motor()
	

	"""
	[{"tick": 10, "jerk": 0.0001},...]
	"""
	def tick_gen_circle(self, xyz_s, xyz_e, p0, r1, r4, d, r, tick_jerk_l , t_s=250, d_0 = 0, v_0 = 0, a_0 = 0, j_0 = 0):
		rtn = []
		for tick_jerk in tick_jerk_l:
			t_s_total = 0
			for j in range(math.floor(tick_jerk["tick"]/t_s)-1):				
				rtn.append({
					"tick": t_s,
					"jerk": tick_jerk["jerk"],
					})
				t_s_total += t_s
			rtn.append({
				"tick": tick_jerk["tick"] - t_s_total,
				"jerk": tick_jerk["jerk"] 
				})

		_xyz_s = np.array(xyz_s)
		_motor_s = np.array(self.kinematic.inverse(_xyz_s)[0])
		_motor_s = np.array(self.motor.joint_to_motor(_motor_s))

		for i in range(len(rtn)):
			rtn[i]["xyz_s"] = np.array(_xyz_s)
			rtn[i]["d_0"] = d_0
			rtn[i]["accel"] = a_0
			rtn[i]["vel"] = v_0
			rtn[i]["d"] = rtn[i]["tick"] * v_0 + rtn[i]["tick"]*(rtn[i]["tick"]-1) * a_0/2 + rtn[i]["tick"]*(rtn[i]["tick"]-1)*(rtn[i]["tick"]-2)*rtn[i]["jerk"]/6 
			
			# next itteration
			d_0 += rtn[i]["d"]
			v_0 +=  rtn[i]["tick"] * a_0 +      rtn[i]["tick"]*(rtn[i]["tick"]-1)* rtn[i]["jerk"]/2 
			a_0 += rtn[i]["tick"]  * rtn[i]["jerk"]

			t_0 = d_0 / r
			xyz = p0 + math.cos(t_0)*r1 + math.sin(t_0)*r4
			ab = xyz_s[3:] + (d_0/d)* (xyz_e[3:] - xyz_s[3:]) 
			rtn[i]["xyz_e"] =  np.append( xyz , ab )
			_xyz_s = np.array(rtn[i]["xyz_e"])			

			# motor
			rtn[i]["motor_s"] = np.array(_motor_s)
			rtn[i]["motor_e"] = np.array(self.kinematic.inverse(rtn[i]["xyz_e"])[0])
			rtn[i]["motor_e"] = np.array(self.motor.joint_to_motor(rtn[i]["motor_e"]))
			_motor_s = np.array(rtn[i]["motor_e"])

			# motor parameter
			rtn[i]["d_motor"] = np.linalg.norm(rtn[i]["motor_s"]-rtn[i]["motor_e"])
			rtn[i]["a_m"] = np.array(rtn[i]["motor_s"])
			rtn[i]["b_m"] = (np.array(rtn[i]["motor_e"]) - np.array(rtn[i]["motor_s"]))/rtn[i]["d_motor"]
			rtn[i]["jerk_m"] = rtn[i]["jerk"] * rtn[i]["d_motor"] / rtn[i]["d"]    
			rtn[i]["accel_m"] = rtn[i]["accel"] * rtn[i]["d_motor"] / rtn[i]["d"]
			rtn[i]["vel_m"] = rtn[i]["vel"] * rtn[i]["d_motor"] / rtn[i]["d"]
			rtn[i]["njerk"] = rtn[i]["jerk"] / rtn[i]["d"]    
			rtn[i]["naccel"] = rtn[i]["accel"] / rtn[i]["d"]
			rtn[i]["nvel"] = rtn[i]["vel"]/ rtn[i]["d"]

		return rtn



	"""
	[{"tick": 10, "jerk": 0.0001},...]
	"""
	def tick_gen_line(self, xyz_s, xyz_e, tick_jerk_l , t_s=250, d_0 = 0, v_0 = 0, a_0 = 0, j_0 = 0):
		rtn = []
		for tick_jerk in tick_jerk_l:
			t_s_total = 0
			for j in range(math.floor(tick_jerk["tick"]/t_s)-1):				
				rtn.append({
					"tick": t_s,
					"jerk": tick_jerk["jerk"],
					})
				t_s_total += t_s
			rtn.append({
				"tick": tick_jerk["tick"] - t_s_total,
				"jerk": tick_jerk["jerk"] 
				})

		d = np.linalg.norm(xyz_e-xyz_s)

		for i in range(len(rtn)):
			rtn[i]["xyz_s"] = np.array(xyz_s + d_0*(xyz_e - xyz_s)/d)
			#print(rtn[i]["xyz_s"] )
			rtn[i]["d_0"] = d_0
			rtn[i]["accel"] = a_0
			rtn[i]["vel"] = v_0
			rtn[i]["d"] = rtn[i]["tick"] * v_0 + rtn[i]["tick"]*(rtn[i]["tick"]-1) * a_0/2 + rtn[i]["tick"]*(rtn[i]["tick"]-1)*(rtn[i]["tick"]-2)*rtn[i]["jerk"]/6 
			
			# next itteration
			d_0 += rtn[i]["d"]
			v_0 +=  rtn[i]["tick"] * a_0 +      rtn[i]["tick"]*(rtn[i]["tick"]-1)* rtn[i]["jerk"]/2 
			a_0 += rtn[i]["tick"]  * rtn[i]["jerk"]
			
			rtn[i]["xyz_e"] = np.array(xyz_s + d_0*(xyz_e - xyz_s)/d) 
			#print(rtn[i]["xyz_e"], rtn[i]["xyz_s"]  ) 
			# motor
			rtn[i]["motor_s"] = np.array(self.kinematic.inverse(rtn[i]["xyz_s"])[0])
			rtn[i]["motor_s"] = np.array(self.motor.joint_to_motor(rtn[i]["motor_s"]))
			rtn[i]["motor_e"] = np.array(self.kinematic.inverse(rtn[i]["xyz_e"])[0])
			rtn[i]["motor_e"] = np.array(self.motor.joint_to_motor(rtn[i]["motor_e"]))
			
			# motor parameter
			rtn[i]["d_motor"] = np.linalg.norm(rtn[i]["motor_s"]-rtn[i]["motor_e"])
			rtn[i]["a_m"] = np.array(rtn[i]["motor_s"])
			rtn[i]["b_m"] = (np.array(rtn[i]["motor_e"]) - np.array(rtn[i]["motor_s"]))/rtn[i]["d_motor"]
			rtn[i]["jerk_m"] = rtn[i]["jerk"] * rtn[i]["d_motor"] / rtn[i]["d"]    
			rtn[i]["accel_m"] = rtn[i]["accel"] * rtn[i]["d_motor"] / rtn[i]["d"]
			rtn[i]["vel_m"] = rtn[i]["vel"] * rtn[i]["d_motor"] / rtn[i]["d"]
			rtn[i]["njerk"] = rtn[i]["jerk"] / rtn[i]["d"]    
			rtn[i]["naccel"] = rtn[i]["accel"] / rtn[i]["d"]
			rtn[i]["nvel"] = rtn[i]["vel"]/ rtn[i]["d"]

		return rtn

	def cmd_to_dorna(self, cmd):
		dorna = []
		for c in cmd:
			dorna.append({
				"t": c["tick"],
				"jerk": c["jerk_m"],
				"accel": c["accel_m"],
				"vel": c["vel_m"],
				"q": 0,  
				"a0": c["a_m"][0],
				"a1": c["a_m"][1],
				"a2": c["a_m"][2],
				"a3": c["a_m"][3],
				"a4": c["a_m"][4],
				"a5": 0, 
				"a6": 0,
				"a7": 0,   
				"b0": c["b_m"][0],
				"b1": c["b_m"][1],
				"b2": c["b_m"][2],
				"b3": c["b_m"][3],
				"b4": c["b_m"][4],
				"b5": 0, 
				"b6": 0,
				"b7": 0,
				"m0": c["motor_e"][0],
				"m1": c["motor_e"][1],
				"m2": c["motor_e"][2],
				"m3": c["motor_e"][3],
				"m4": c["motor_e"][4],
				"m5": 0,
				"m6": 0,
				"m7": 0,
				"njerk": c["njerk"], 
				"naccel": c["naccel"],
				"nvel": c["nvel"],			
				})
		return dorna

	"""
	give start and end in xyz plane
	d_s is the sampeling distance, default is 1mm

	start sampeling in xyz and convert it to joint space
	sum all the dsitance in joint space and find the associated profile

	cmd: qmove
	a0, .. ,a7
	b0,..., b7
	q_max,
	j0, a0, v0
	"""
	def line(self, xyz_s, xyz_e, j, a, v, v_0 = 0):
		xyz_s = np.array(xyz_s)
		xyz_e = np.array(xyz_e)

		# distance
		d = np.linalg.norm(xyz_e-xyz_s)

		# find profile for this distance
		result = profile.profile_1(j, a, v, d, v_0)
		#jerk_l = [x[1] for x in result["tick_jerk"]]
		#tick_l = [x[0] for x in result["tick_jerk"]]
		#a, v, q = profile.tick(jerk_l, tick_l, v_0)
		
		# tick jerk
		tick_jerk_l = profile.tick_jerk(result["tick_jerk"])

		cmd = self.tick_gen_line(xyz_s, xyz_e, tick_jerk_l)
		dorna = self.cmd_to_dorna(cmd)
		
		return dorna
		

	"""
	xyz_s, 
	"""
	def circle(self, xyz_s, xyz_m, xyz_e, j, _a, v, v_0 = 0, turn = 0):
		# make a copy
		xyz_s = np.array(xyz_s)
		xyz_m = np.array(xyz_m)
		xyz_e = np.array(xyz_e)
		
		# xyz space
		A = np.array(xyz_s[0:3])
		B = np.array(xyz_e[0:3])
		C = np.array(xyz_m[0:3])

		a = A - C
		b = B - C

		# center of circle
		_a_b_cross = np.cross(a,b)
		p0 =  (np.linalg.norm(a)**2)*b - (np.linalg.norm(b)**2)*a
		p0 = np.cross(p0,_a_b_cross) / (2*np.linalg.norm(_a_b_cross)**2)
		p0 = p0 + C
		# radius
		r = np.linalg.norm(A-p0)

		
		# traverse vector
		r1 = A - p0
		r2 = B - p0
		r3 = r2 -  (r1* np.inner(r1, r2)/(r**2)) 
		r4 = r3 *(r/np.linalg.norm(r3))
		r5 = C - p0

		# find t_b and t_c
		cos_t_b = min(max(-1, np.inner(r1,r2) / r**2), 1) 
		t_b = np.arccos(cos_t_b)
		
		# adjust r4 and t_b
		if np.inner(r5,r4) < 0:
			r4 = -r4
			t_b = 2*math.pi - t_b
		elif np.arccos(min(max(-1, np.inner(r1,r5) / r**2), 1)) > t_b:
			r4 = -r4
			t_b = 2*math.pi - t_b

		# distance
		d = r * (2*turn*math.pi + t_b)

		print(p0, r1, r4, d, r)
		# profile
		result = profile.profile_1(j, _a, v, d, v_0)

		# tick jerk
		tick_jerk_l = profile.tick_jerk(result["tick_jerk"])

		cmd = self.tick_gen_circle(xyz_s, xyz_e, p0, r1, r4, d, r, tick_jerk_l)


		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')

		# For each set of style and range settings, plot n random points in the box
		# defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
		for r in cmd:
			xs = r["xyz_s"][0]
			ys = r["xyz_s"][1]
			zs = r["xyz_s"][2]
			ax.scatter(xs, ys, zs, c="b", marker="o")
		plt.show()

		dorna = self.cmd_to_dorna(cmd)
		
		return dorna


	def line_type_3(self, xyz_s, xyz_e, j, a, v, v_0 = 0, d_s = 1):
		xyz_s = np.array(xyz_s)
		xyz_e = np.array(xyz_e)

		# distance
		d = np.linalg.norm(xyz_e-xyz_s)

		# get the motors
		motor_sample = []
		joint_sample = []
		for i in range(math.floor(d/d_s)+1):
			xyz = xyz_s + i*d_s*(xyz_e - xyz_s)/d
			joint_xyz = np.array(self.kinematic.inverse(xyz)[0])
			motor_xyz =np.array(self.motor.joint_to_motor(joint_xyz))
			motor_sample.append(motor_xyz)
			joint_sample.append(joint_xyz)

		# add laste sample
		motor_sample.pop(-1)
		xyz = xyz_e	
		joint_xyz = np.array(self.kinematic.inverse(xyz)[0])
		motor_xyz =np.array(self.motor.joint_to_motor(joint_xyz))
		motor_sample.append(motor_xyz)


		# find total distance
		d_motor = [np.linalg.norm(motor_sample[i+1]-motor_sample[i]) for i in range(len(motor_sample) - 1)]
		# find profile for this distance
		result = profile.profile_1(j, a, v, sum(d_motor), v_0)
		list_jerk = [x[1] for x in result["tick_jerk"]]
		list_tick = [x[0] for x in result["tick_jerk"]]
		_a, _v, q = profile.tick(list_jerk, list_tick, v_0)
		profile.plot(_a, _v, q)


		tick_jerk = profile.tick_jerk(result["tick_jerk"])
		t_j = []
		t = 0 
		for i in range(len(tick_jerk)):
			#t += tick_jerk[i]["tick"]
			jerk = tick_jerk[i]["jerk"]
			q_max = q[t]
			cmd = { 
				"q_max": q_max,
				"jerk": jerk,
				}			
			t_j.append(cmd)
			#print(json.dumps(cmd))			
			t += tick_jerk[i]["tick"]
		t_j[0]["q"] = 0

		# {"a": , "b":, "q": } => p(t) = A + q(t)B
		a_b = []
		_a_b = []
		q_max = 0
		for i in range(len(d_motor)):
			A = motor_sample[i] - q_max*(motor_sample[i+1] - motor_sample[i])/d_motor[i]
			B = (motor_sample[i+1] - motor_sample[i])/d_motor[i]
			cmd = {
				"q_max": q_max, 
				"a0": A[0],
				"a1": A[1],
				"a2": A[2],
				"a3": A[3],
				"a4": A[4],
				"a5": 0, 
				"a6": 0,
				"a7": 0,   
				"b0": B[0],
				"b1": B[1],
				"b2": B[2],
				"b3": B[3],
				"b4": B[4],
				"b5": 0, 
				"b6": 0,
				"b7": 0,
			}
			_a_b.append(dict(cmd))
			a_b.append(dict(cmd))
			q_max += d_motor[i] 
			#print(json.dumps(cmd))

		rtn = []

		#print("###size: ", len(a_b), len(t_j))
		while t_j and a_b:				
			
			if t_j[0]["q_max"] < a_b[0]["q_max"]:
				rtn.append(t_j.pop(0))
			elif t_j[0]["q_max"] > a_b[0]["q_max"]:  
				rtn.append(a_b.pop(0))
			else:
				a_b[0]["jerk"] = t_j[0]["jerk"]
				rtn.append(a_b.pop(0))
				t_j.pop(0)
		
		#print("###a_b: ", a_b)
		#print("###t_j: ", t_j)
		while t_j:
			rtn.append(t_j.pop(0))

		while a_b:
			rtn.append(a_b.pop(0))

		#print("size: ", len(rtn))
		#print("size: ", len(a_b))
		#print("size: ", len(t_j))

		"""
		plot j0 velocity
		"""
		xyz = []
		joint_index = "1"
		j = 0
		try:
			for _q in q:
				j += 1
				if j == 100:
					j = 0
					print(j)
					for i in range(len(_a_b)-1,-1,-1):
						if _q >= _a_b[i]["q_max"]:
							_m = [_a_b[i]["a0"] + _q*_a_b[i]["b0"],
									_a_b[i]["a1"] + _q*_a_b[i]["b1"],
									_a_b[i]["a2"] + _q*_a_b[i]["b2"],
									_a_b[i]["a3"] + _q*_a_b[i]["b3"],
									_a_b[i]["a4"] + _q*_a_b[i]["b4"]
									]
							_j = self.motor.motor_to_joint(_m)
							_xyz = self.kinematic.forward(_j) 		 

							xyz.append(_xyz)
							break										

		except Exception as ex:
			pass
			print("error: ", ex)
			#print(len(_a_b), i)

		#print("a, b, q, m: ",_a + _qq*_b, _qq, _q)
		y = [xyz[i][1] for i in range(len(xyz))] 
		v_y = [y[i+1] - y[i] for i in range(len(y)-1) ]
		a_y = [v_y[i+1] - v_y[i] for i in range(len(v_y)-1) ]

		print(y)
		plt.figure(1)

		plt.subplot(311)
		plt.plot(a_y, 'o-')
		#plt.title("a_end "+ str(a[-1]) )
		
		plt.subplot(312)
		plt.plot(v_y, 'o-')
		#plt.title("v_end " + str(v[-1]))

		plt.subplot(313)
		plt.plot(y, 'o-')

		plt.show()



		m = []
		joint_index = "1"
		try:
			for _q in q:
				for i in range(len(_a_b)-1,-1,-1):
					if _q >= _a_b[i]["q_max"]:
						_a =_a_b[i]["a"+joint_index]
						_b = _a_b[i]["b"+joint_index]  
						_qq = _q
						m.append(_a + _q*_b )
						break										

		except Exception as ex:
			pass
			#print(ex, i)
			#print(len(_a_b), i)

		#print("a, b, q, m: ",_a + _qq*_b, _qq, _q)
		v = [m[i+1] - m[i] for i in range(len(m)-1) ]
		a = [v[i+1] - v[i] for i in range(len(v)-1) ]

		plt.figure(1)

		plt.subplot(411)
		plt.plot(a, 'o-')
		#plt.title("a_end "+ str(a[-1]) )
		
		plt.subplot(412)
		plt.plot(v, 'o-')
		#plt.title("v_end " + str(v[-1]))

		plt.subplot(413)
		plt.plot(m, 'o-')
		plt.title("m_end "+ str(m[-1]))

		plt.subplot(414)
		ms = [x[0] for x in motor_sample]
		plt.plot(ms, 'o-')
		plt.title("m_end "+ str(ms[-1]))
		plt.show()
		 


		"""
		"""

		for i in range(len(rtn)-1):
			rtn[i]["q_max"] = rtn[i+1]["q_max"] 
		rtn[-1]["q_max"] = sum(d_motor)

		rtn[0]["q"] = 0
		rtn[0]["vel"] = 0
		rtn[0]["accel"] = 0
		"""
		with open("cmds.txt", 'w') as out_file:
	
		out_file.write(parsed_line)		
		"""
		#print("aliali123")
		
		for r in rtn:
			print(json.dumps(r))
		

		"""
		for r in _a_b:
			print(r)
		"""
		return rtn



	def line_type_2(self, xyz_s, xyz_e, j, a, v, t_s = 100, v_0 = 0):
		xyz_s = np.array(xyz_s)
		xyz_e = np.array(xyz_e)

		# distance
		d = np.linalg.norm(xyz_e-xyz_s)

		"""
		profile
		return: {"tick_jerk": t_j, "v_e": v_0, "d": d_m}
		t_j = [[t_1, j_m], [t_2, 0], [t_1, -j_m], [t_4, 0], [t_1,-j_m], [t_2, 0], [t_1, j_m]]
		"""
		result = profile.profile_1(j, a, v, d)
		list_jerk = [x[1] for x in result["tick_jerk"]]
		list_tick = [x[0] for x in result["tick_jerk"]]
		a, v, q = profile.tick(list_jerk, list_tick, v_0)
		profile.plot(a, v, q)

		# number of segments
		number_of_segments = math.floor(len(q)/t_s)

		"""
		for each segment find
		t is the number of ticks, so ticks changes from 0 to t-1
		v_1, a_1, j_1 are the initiaal v, a and j 
		"""
		list_tvaj = [] 		
		
		# v[-1]
		v_s = 0

		# joint p[0] of next segment
		p_e = np.array(self.kinematic.inverse(xyz_s)[0])
		p_e = np.array(self.motor.joint_to_motor(p_e))

		total_distance = 0

		for i in range(number_of_segments):
			# update p_0 (p start)
			p_s = np.array(p_e)
			
			# update p_e and p_0_f
			if i < number_of_segments - 1:
				p_e = np.array(self.kinematic.inverse(xyz_s + q[(i+1)*t_s]*(xyz_e - xyz_s)/d)[0])
				p_e =np.array(self.motor.joint_to_motor(p_e))

				p_e_1 = np.array(self.kinematic.inverse(xyz_s + q[(i+1)*t_s - 1]*(xyz_e - xyz_s)/d)[0])
				p_e_1 =np.array(self.motor.joint_to_motor(p_e_1)) 
				t = t_s 
			else:
				p_e = np.array(self.kinematic.inverse(xyz_s + q[- 1]*(xyz_e - xyz_s)/d)[0])
				p_e =np.array(self.motor.joint_to_motor(p_e))

				p_e_1 = np.array(self.kinematic.inverse(xyz_s + q[- 2]*(xyz_e - xyz_s)/d)[0])
				p_e_1 =np.array(self.motor.joint_to_motor(p_e_1))
				t = len(q) - i*t_s				

			# q(t_1)
			d_e = np.linalg.norm(p_e - p_s)
			total_distance += d_e 
			
			#v[T-1]
			v_e_1 = np.linalg.norm(p_e - p_e_1)

			inv = np.linalg.inv([[t*(t-1)/2 , t*(t-1)*(t-2)/6], [t-1, (t-1)*(t-2)/2]])
			ans = np.matmul(inv, np.array([[d_e - t* v_s],[v_e_1 - v_s]]))

			# t, v, a, j
			tvaj = [t, v_s, ans[0][0], ans[1][0], p_s, p_e]
			list_tvaj.append(tvaj)
			
			cmd = {
				"cmd": "raw",
				"t": t,
				"v": v_s,
				"a": ans[0][0],
				"j": ans[1][0],
				"p0": p_s[0],
				"p1": p_s[1],
				"p2": p_s[2],
				"p3": p_s[3],
				"p4": p_s[4],
				"p5": 0,
				"p6": 0,
				"p7": 0,
				"q0": p_e[0],
				"q1": p_e[1],
				"q2": p_e[2],
				"q3": p_e[3],
				"q4": p_e[4],
				"q5": 0,
				"q6": 0,
				"q7": 0,
			}
			
			print(json.dumps(cmd))
			#print(tvaj[0], tvaj[1], tvaj[2], tvaj[3],  p_s[0], p_s[1], p_s[2], p_s[3], p_s[4], 0, 0,0, p_e[0], p_e[1], p_e[2], p_e[3], p_e[4], 0, 0, 0)


			# update v_s
			v_s += t*ans[0][0] + t*(t-1)* ans[1][0]/2

		return list_tvaj

def main_1():
	j = 1.0e-11
	a = 1.0e-6
	v = 0.05
	v_0 = 0
	d_s = 1
	joint_s = [0, 0, 0, 0, 0] 
	d = 150
	
	p = path()
	xyz_s = p.kinematic.forward(joint_s) 
	xyz_e = list(xyz_s)
	xyz_e[0] -= d

	# list_tvaj
	result = p.line_type_2(xyz_s, xyz_e, j, a, v, t_s, v_0)
	result = p.line_type_2(xyz_e, xyz_s, j, a, v, t_s, v_0)
	
	# plot

	v_sample = [r[2] for r in result]
	plt.figure(1)
	
	plt.subplot(111)
	plt.plot(v_sample, 'o-')
	plt.title("q_end "+ str(v_sample[-1]) )
	

	plt.show()	
	

	a_sample = []
	v_l = [0]
	a_l = [0]
	q_l = [0]
	for i in range(len(result)):
		t_m = result[i][0]
		a_l = [result[i][2]]
		j_m = result[i][3]
		for j in range(t_m):
			q_l.append(q_l[-1] + v_l[-1])
			v_l.append(v_l[-1] + a_l[-1])
			a_l.append(a_l[-1] + j_m)
		a_sample += a_l
	plt.figure(1)
	
	plt.subplot(211)
	plt.plot(a_sample, 'o-')
	plt.title("q_end "+ str(a_sample[-1]) )

	plt.subplot(212)
	plt.plot(v_l, 'o-')
	plt.title("v_end " + str(v_l[-1]))
	

	plt.show()

	"""
	# plot
	v_l = [r[1] for r in result]
	
	plt.figure(1)
	
	plt.subplot(111)
	plt.plot(v_l, 'o-')
	plt.title("v_end "+ str(v_l[-1]) )
	
	plt.show()	
	"""
def main_2():
	joint = [45, 45, 45, 45, 45]
	k = kinematic()
	xyz = k.forward(joint)

	print(xyz)
	print(k.inverse([400.05469,   0,      200.09099,   0,        0.     ]))


def main_line():
	j = 1.0e-12
	a = 1.0e-7
	v = 0.001
	v_0 = 0
	d_s = 5
	joint_0 = [0, 0, 0, 0, 0]
	joint_1 = [0, 60, -60, 0, 0]
	
	p = path()
	xyz_0 = p.kinematic.forward(joint_0) 
	xyz_1 = p.kinematic.forward(joint_1)

	# list_tvaj
	result = []
	result += p.line(xyz_0, xyz_1, j, a, v, v_0)
	result += p.line(xyz_1, xyz_0, j, a, v, v_0)

	f = open("cmds.txt", "w")
	for r in result:
		f.write(json.dumps(r) + "\n")

	f.close()	
	#result = p.line(xyz_p1, xyz_p2, j, a, v, v_0, d_s)
	#result = p.line(xyz_p2, xyz_p3, j, a, v, v_0, d_s)
	#result = p.line(xyz_p3, xyz_p4, j, a, v, v_0, d_s)	
	#result = p.line(xyz_p4, xyz_s, j, a, v, v_0, d_s)		

def main_circle():
	j = 1.0e-11
	a = 1.0e-6
	v = 0.002
	v_0 = 0
	turn = 0
	joint_0 = [0, 44, -111, 67, 0]

	p = path()

	xyz_0 = p.kinematic.forward(joint_0) 

	xyz_1 = np.array(xyz_0)
	xyz_1[0] += 100

	xyz_2 = np.array(xyz_0)
	xyz_2[0] += 60
	xyz_2[1] += 60
	xyz_2[2] += 60

	# list_tvaj
	xyz_s = [400.0, 0.0 , 180.0, 0.0, 0.0]
	xyz_m = [300.0, 0.0 , 180.0, 0.0, 0.0]
	xyz_e = [349.0, -50.0 , 180.0, 0.0, 0.0]
	result = p.circle(xyz_s, xyz_m, xyz_e, j, a, v, v_0, turn)

	f = open("cmds.txt", "w")
	for r in result:
		f.write(json.dumps(r) + "\n")

	f.close()	
	#result = p.line(xyz_p1, xyz_p2, j, a, v, v_0, d_s)
	#result = p.line(xyz_p2, xyz_p3, j, a, v, v_0, d_s)
	#result = p.line(xyz_p3, xyz_p4, j, a, v, v_0, d_s)	
	#result = p.line(xyz_p4, xyz_s, j, a, v, v_0, d_s)		

	

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

if __name__ == '__main__':

	main_circle()