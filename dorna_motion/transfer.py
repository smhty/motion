import numpy as np
import math

class quat(object):
	"""docstring for quat"""
	def __init__(self):
		super(quat, self).__init__()

	def inverse(self, q):
		return[q[0],
				-q[1],
				-q[2],
				-q[3]]

	def mul(self, q, p):
		return [q[0]*p[0]-q[1]*p[1]-q[2]*p[2]-q[3]*p[3],
				q[0]*p[1]+q[1]*p[0]+q[2]*p[3]-q[3]*p[2],
				q[0]*p[2]+q[2]*p[0]-q[1]*p[3]+q[3]*p[1],
				q[0]*p[3]+q[3]*p[0]+q[1]*p[2]-q[2]*p[1]]

	def sum(self, q, p):
		return [q[0]+p[0],
				q[1]+p[1],
				q[2]+p[2],
				q[3]+p[3]]

	def quat_to_point(self, v):
		return [v[1], v[2], v[3]]

	def point_to_quat(self, v):
		return[0,
				v[0],
				v[1],
				v[2]]

	# theta degree around axis u
	def axis_rot_to_quat(self, u, theta):
		cos = np.cos(theta/2)
		sin = np.sin(theta/2)
		L = 1/(math.sqrt(u[0]**2 + u[1]**2 + u[2]**2))
		return [cos,
				sin*u[0],
				sin*u[1],
				sin*u[2]] 

	# rotation: rotate p (pure_quat) around quat q 
	def active_rot(self, pure_quat, q):
		inv = self.inverse(q)
		step_1 = self.mul(q, pure_quat)
		return self.mul(step_1, inv)

	def quat_to_rot_vec(self, q):
		return 0

	def rot_vec_to_quat(self, v):
		return 0


class rot(object):
	"""docstring for rotation_matrix"""
	def __init__(self):
		super(rot, self).__init__()
	
	def x(self, t):
		return np.matrix([[1, 0 , 0], [0, np.cos(t), -np.sin(t)], [0, np.sin(t), np.cos(t)]])	

	def y(self, t):
		return np.matrix([[np.cos(t), 0 , np.sin(t)], [0, 1, 0], [-np.sin(t), 0, np.cos(t)]])	

	def z(self, t):
		return np.matrix([[np.cos(t), -np.sin(t) , 0], [np.sin(t), np.cos(t), 0], [0, 0, 1]])	
