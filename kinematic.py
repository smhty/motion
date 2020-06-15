import math

"""
use kinematic.jpg 
l = [
		[l_00, 0, l_02],
		[l_10, 0, 0],
		[l_20, 0, 0],
		[0, 0, 0],
		[0,  0, l_42]
	]
"""
class kinematic(object):
	"""docstring for kinematic"""
	def __init__(self, l):
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

		except:
			pass

		return joint


if __name__ == '__main__':
	l = [[3.76, 0, 8.11],
		[8, 0 , 0],
		[6, 0 , 0],
		[0, 0 , 0],
		[0, 0, 1.72]
	]

	k = kinematic(l)	

	joint = [45 for i in range(5)]
	xyz = k.forward(joint)
	joint = k.inverse(xyz)

	print(xyz)
	print(joint)