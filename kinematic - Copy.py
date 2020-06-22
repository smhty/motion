import math
import profile
import numpy as np
import matplotlib.pyplot as plt


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
q[n] =(n+1) * v + n*(n+1) * a/2 + (n-1)*n*(n+1)*j/6
v[n] =        v +     (n+1) * a +      n*(n+1)* j/2 	
a[n] =                        a +       (n+1) * j 	

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
	
	def line(self, xyz_s, xyz_e, j, a, v, t_s = 100, v_0 = 0):
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
		v_m_1 = 0

		# joint p[0] of next segment
		p_0_f = self.kinematic.inverse(xyz_s + q[0]*(xyz_e - xyz_s)/d)[0]
		p_0_f = np.array(p_0_f)

		total_distance = 0
		for i in range(number_of_segments):
			# update p_0 (p start)
			p_s = np.array(p_0_f)
			
			# update p_e and p_0_f
			if i < number_of_segments - 1:
				p_e = np.array(self.kinematic.inverse(xyz_s + q[(i+1)*t_s - 1]*(xyz_e - xyz_s)/d)[0]) 
				p_0_f = np.array(self.kinematic.inverse(xyz_s + q[(i+1)*t_s]*(xyz_e - xyz_s)/d)[0]) 
				t = t_s 
			else:
				p_e = np.array(self.kinematic.inverse(xyz_s + q[- 1]*(xyz_e - xyz_s)/d)[0])  
				p_0_f = np.array(self.kinematic.inverse(xyz_e)[0])
				t = len(q) - i*t_s				

			# q(t_1)
			q_t_1 = np.linalg.norm(p_e - p_s) 
			total_distance += q_t_1 
			
			#v[T-1]
			v_t_1 = np.linalg.norm(p_0_f - p_e)

			inv = np.linalg.inv([[t*(t-1)/2 , t*(t-1)*(t-2)/6], [t, t*(t-1)/2]])
			#inv = np.linalg.inv([[t*(t+1)/2 , t*(t+1)*(t-1)/2], [t , t*(t+1)/2]]) 
			#ans = np.matmul(inv, np.array([[q_t_1 - v_m_1],[v_t_1 - v_m_1]]))
			ans = np.matmul(inv, np.array([[q_t_1 - t* v_m_1],[v_t_1 - v_m_1]]))

			# t, v, a, j
			list_tvaj.append([t, v_m_1, ans[0][0], ans[1][0]])
			print(t,",",ans[0][0],",",ans[1][0])
			# update v_m
			#v_m_1 += (t-1) * ans[0][0] +    t*(t-1)* ans[1][0]/2
			v_m_1 = v_t_1  

		return list_tvaj

def main_1():
	j = 2.0e-10
	a = 1.0e-5
	v = 0.45
	v_0 = 0
	t_s = 200
	joint_s = [0, 0, 0, 0, 0] 
	d = 100
	
	p = path()
	xyz_s = p.kinematic.forward(joint_s) 
	xyz_e = list(xyz_s)
	xyz_e[0] -= d
	
	# list_tvaj
	result = p.line(xyz_s, xyz_e, j, a, v, t_s, v_0)

	
	# plot
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
	
	q_l = [math.ceil(400*x) for x in q_l]
	plt.figure(1)
	
	plt.subplot(211)
	plt.plot(q_l, 'o-')
	plt.title("q_end "+ str(q_l[-1]) )

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
	joint = [0, 0, 0, 0, 0]
	k = kinematic()
	xyz = k.forward(joint)

	print(k.inverse(xyz))

if __name__ == '__main__':

	main_1()