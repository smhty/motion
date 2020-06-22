import math
import profile
import numpy as np
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

		except Exception as ex:
			print("AAAAAAAAAAAAAAaa: ", ex)
			pass

		return joint

	def inverse_copy(self, xyz):
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
			
			# not valid
			if L > (self.l[1][0] + self.l[2][0]) - self.e or self.l[1][0] > (self.l[2][0] + L) - self.e: # in this case there is no solution
				return joint


			theta_l1_L = math.acos((self.l[1][0] ** 2 + L ** 2 - self.l[2][0] ** 2) / (2 * self.l[1][0] * L))  # l1 angle to L
			theta_L_x = math.atan2(z, x)  # L angle to x axis
			theta_1_0 = theta_L_x + theta_l1_L
			theta_1_1 = theta_L_x - theta_l1_L

			theta_l1_l2 = math.acos((self.l[1][0] ** 2 + self.l[2][0] ** 2 - L ** 2) / (2 * self.l[1][0] * self.l[2][0]))  # l1 angle to l2
			theta_2_0 = theta_l1_l2 - math.pi
			theta_2_1 = -(theta_l1_l2 - math.pi)

			theta_3_0 = a - theta_1_0 - theta_2_0
			theta_3_1 = a - theta_1_1 - theta_2_1

			joint = [[theta_0, theta_1_0, theta_2_0, theta_3_0, b],[theta_0, theta_1_1, theta_2_1, theta_3_1, b]]
			joint = [[math.degrees(j) for j in joint[0]], [math.degrees(j) for j in joint[1]]]

		except Exception as ex:
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
q[n] =n * v + n*(n-1) * a/2 + n*(n-1)*(n+2)*j/6
v[n] =        v +     n * a +      n*(n-1)* j/2 	
a[n] =                        a +       n * j 	

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
	give start and end in xyz plane
	d_s is the sampeling distance, default is 1mm

	start sampeling in xyz and convert it to joint space
	sum all the dsitance in joint space and find the associated profile

	"""
	def line(self, xyz_s, xyz_e, j, a, v, d_s = 1):
		xyz_s = np.array(xyz_s)
		xyz_e = np.array(xyz_e)

		# distance
		d = np.linalg.norm(xyz_e-xyz_s)

		# get the motors
		motor_sample = []
		for i in range(math.floor(d/d_s)):
			xyz = xyz_s + i*d_s*(xyz_e - xyz_s)/d
			joint_xyz = np.array(self.kinematic.inverse(xyz)[0])
			motor_xyz =np.array(self.motor.joint_to_motor(joint_xyz))
			motor_sample.append(motor_xyz)
		xyz = xyz_e	
		joint_xyz = np.array(self.kinematic.inverse(xyz)[0])
		motor_xyz =np.array(self.motor.joint_to_motor(joint_xyz))
		motor_sample.append(motor_xyz)

		# find total distance
		d_motor = [np.linalg.norm(motor_sample[i+1]-motor_sample[i]) for i in range(len(motor_sample) - 1)]
		d_motor = sum(d_motor)

		# find profile for this distance
		result = profile.profile_1(j, a, v, d_motor, v_0)
		list_jerk = [x[1] for x in result["tick_jerk"]]
		list_tick = [x[0] for x in result["tick_jerk"]]
		a, v, q = profile.tick(list_jerk, list_tick, v_0)


	def line_2(self, xyz_s, xyz_e, j, a, v, t_s = 100, v_0 = 0):
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
	j = 2.0e-11
	a = 1.0e-7
	v = 0.01
	v_0 = 0
	t_s = 1000
	joint_s = [0, 30, -30, 0, 0] 
	d = 70
	
	p = path()
	xyz_s = p.kinematic.forward(joint_s) 
	xyz_e = list(xyz_s)
	xyz_e[0] -= d

	# list_tvaj
	result = p.line(xyz_s, xyz_e, j, a, v, t_s, v_0)
	result = p.line(xyz_e, xyz_s, j, a, v, t_s, v_0)
	
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

if __name__ == '__main__':

	main_1()