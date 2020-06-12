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
		
	"""
	forward kinematics: joint to xyz
	"""
	def fk(self, joint):
		# joint to radian
		joint = [math.radians(j) for j in joint]

		# first we find x, y, z assuming base rotation is zero (j_0 = 0). Then we rotate everything
		# then we rotate the robot around z axis for j_0
		tmp_d = self.l[0][0] + self.l[1][0] * math.cos(joint[1]) + self.l[2][0] * math.cos(joint[1] + joint[2]) + self.l[4][2] * math.cos(joint[1]+joint[2]+joint[3])
		x = tmp_d * math.cos(joint[0])
		y = tmp_d * math.sin(joint[0])
		z = self.l[0][2] + self.l[1][0] * math.sin(joint[1]) + self.l[2][0] * math.sin(joint[1] + joint[2]) + self.l[4][2] * math.sin(joint[1] + joint[2] + joint[3])
		a = joint[1] + joint[2] + joint[3]
		b = joint[4]
		
		return [x, y, z, a, b]


	"""
	inverse kinematics: xyz to joint
	"""
	def ik(self, xyz):
		joint = []
		# make a copy
		xyz = [x for x in xyz]
		xyz[3] = math.radians(xyz[3])
		xyz[4] = math.radians(xyz[4])		
		
		try:
			# j0
			joint.append(math.atan2(y, x))
		
			# next we assume base is not rotated and everything lives in x-z plane
			xyz[0] = math.sqrt(xyz[0] ** 2 + xyz[1] ** 2)

			# next we update x and z based on base dimensions and hand orientation
			xyz[0] -= (self.l[0][0] + self.l[4][2] * math.cos(xyz[3]))
			xyz[2] -= (self.l[0][2] + self.l[4][2] * math.sin(xyz[3]))

			# at this point x and z are the summation of two vectors one from lower arm and one from upper arm of lengths l1 and l2
			# let L be the length of the overall vector
			# we can calculate the angle between l1 , l2 and L
			L = math.sqrt(xyz[0] ** 2 + xyz[2] ** 2)
		
			# not valid
			if L > (self.l[1][0] + self.[2][0]) or self.l[1][0] > (self.l[2][0] + L) or self.l[2][0] > (self.l[1][0] + L):  # in this case there is no solution
				return None


			teta_l1_L = math.acos((self.l[1][0] ** 2 + L ** 2 - self.l[2][0] ** 2) / (2 * self.l[1][0] * L))  # l1 angle to L
			teta_L_x = math.atan2(xyz[2], xyz[0])  # L angle to x axis
			joint.append(teta_l1_L + teta_L_x)
			# note that the other solution would be to set teta_1 = teta_L_x - teta_l1_L. But for the dynamics of the robot the first solution works better.
			teta_l1_l2 = math.acos((self.l[1][0] ** 2 + self.l[2][0] ** 2 - L ** 2) / (2 * self.l[1][0] * self.l[2][0]))  # l1 angle to l2
			joint.append(teta_l1_l2 - math.pi)
			joint.append(xyz[3] - joint[1] - joint[2])
			joint.append(xyz[4])
			joint = [math.degrees(x) for x in joint]

		return joint