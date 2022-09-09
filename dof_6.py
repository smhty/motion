import numpy as np
import time
import math
import random

"""
Sources:
	http://rasmusan.blog.aau.dk/files/ur5_kinematics.pdf
	https://smartech.gatech.edu/bitstream/handle/1853/50782/ur_kin_tech_report_1.pdf

i | alpha[i-1]        | a[i-1] | d[i] | theta[i]
------------------------------------------------
1 | 0                 | 0      | d[1] | 
2 | alpha[1]= np.pi/2 | a[1]   | 0    | 
3 | 0                 | a[2]   | 0    |
4 | 0                 | a[3]   | d[4] |
5 | alpha[4]= np.pi/2 | 0      | d[5] |
6 | alpha[5]= np.pi/2 | 0      | d[6] |

"""
class dh(object):
	"""docstring for dh"""
	def __init__(self):
		super(dh, self).__init__()

		self.alpha = [0, np.pi/2, 0, 0, np.pi/2, np.pi/2]
		self.a = [0, 1, 2, 1.5, 0, 0]
		self.d = [0, 1, 0, 0, -1, 1, 1]
		self.T_f0_r_base = np.identity(4) # word
		self.T_f_tcp_r6 = np.identity(4) # TCP

	# Ti with respect to i-1 frame 
	def T(self, i, theta):
		ct = np.cos(theta)
		st = np.sin(theta)
		ca = np.cos(self.alpha[i-1])
		sa = np.sin(self.alpha[i-1])

		return np.matrix([
			[ct, -st, 0, self.a[i-1]],
			[st*ca, ct*ca, -sa, -self.d[i]*sa],
			[st*sa, ct*sa, ca, self.d[i]*ca],
			[0, 0, 0, 1]])

	def inv_dh(self, T):
		R = np.matrix([
			[T[0,0], T[1,0], T[2,0]],
			[T[0,1], T[1,1], T[2,1]],
			[T[0,2], T[1,2], T[2,2]]
		])
		S = -np.matmul(R, [[T[0,3]], [T[1,3]], [T[2,3]]])

		return np.matrix([
			[R[0,0], R[0,1], R[0,2], S[0,0]],
			[R[1,0], R[1,1], R[1,2], S[1,0]],
			[R[2,0], R[2,1], R[2,2], S[2,0]],
			[0, 0, 0, 1]	
		])

class dof_6(dh):
	"""docstring for dof_6"""
	def __init__(self):
		super(dof_6, self).__init__()
		self.thr = 0.0000000001

	"""
	joint values are given in radians
	return T_f_tcp_r_base
	"""
	def fw(self, theta):			
		T = self.T_f0_r_base
		for i in range(1, 7):
			T = np.matmul(T, self.T(i, theta[i-1]))
		
		return np.matmul(T, self.T_f_tcp_r6)

	"""
	The robot T_f_tcp_r_base is given
	find all the possible robot orientations 
	"""
	def inv(self, T_f_tcp_r_base):
		# init return
		rtn = []

		T_f6_r0 = np.matmul(T_f_tcp_r_base, self.inv_dh(self.T_f_tcp_r6))
		T_f6_r0 = np.matmul(self.inv_dh(self.T_f0_r_base), T_f6_r0)

		for theta_1 in self.theta_1(T_f6_r0): # 2 x theta_1
			for theta_5 in self.theta_5(T_f6_r0, theta_1): # 2 x theta_5
				theta_6 = self.theta_6(T_f6_r0, theta_1, theta_5) # 1 x theta_6
				for theta_2_3_4 in self.theta_3_2_4(T_f6_r0, theta_1, theta_5, theta_6): # 1 x theta_2, # 2 x theta_3, # 1 x theta_4
					rtn.append([theta_1]+theta_2_3_4+[theta_5, theta_6])
		return rtn

	def theta_1(self, T_f6_r0):
		p5_0 = np.matmul(T_f6_r0, [[0], [0], [-self.d[6]], [1]])
		p5x_0 = p5_0[0,0]
		p5y_0 = p5_0[1,0]
		try:
			alpha = math.asin(-self.d[4]/math.sqrt(p5x_0**2 + p5y_0**2))
			phi_2 = [alpha, math.pi-alpha]

			phi_1 = math.atan2(p5y_0, p5x_0)

			return [phi_1-alpha, phi_1+alpha-math.pi]
		except Exception as ex:
			return []

	def theta_5(self, T_f6_r0, theta_1):
		nom = T_f6_r0[1,3]*np.cos(theta_1)-T_f6_r0[0,3]*np.sin(theta_1)+self.d[4]
		try:
			phi = math.acos(nom/self.d[6])
			return [phi, -phi]
		except Exception as ex:
			return []

	def theta_6(self, T_f6_r0, theta_1, theta_5, theta_6_init=0):
		sgn = 1
		if math.sin(theta_5) == 0:
			return theta_6_init
		elif math.sin(theta_5) < 0:
			sgn = -1
		cos = (-T_f6_r0[0,0]*math.sin(theta_1) + T_f6_r0[1,0]*math.cos(theta_1))/(-sgn)
		sin = (-T_f6_r0[0,1]*math.sin(theta_1) + T_f6_r0[1,1]*math.cos(theta_1))/(-sgn)
		return -math.atan2(sin, cos)

	def theta_3_2_4(self, T_f6_r0, theta_1, theta_5, theta_6):
		rtn = []
		T_f1_r0 = self.T(1, theta_1)
		T_f5_r4 = self.T(5, theta_5)
		T_f6_r5 = self.T(6, theta_6)

		T_f4_r1 = np.matmul(self.inv_dh(T_f1_r0), T_f6_r0)
		T_f4_r1 = np.matmul(T_f4_r1, self.inv_dh(T_f6_r5))
		T_f4_r1 = np.matmul(T_f4_r1, self.inv_dh(T_f5_r4))

		p4x = T_f4_r1[0,3]
		p4z = T_f4_r1[2,3]

		try:
			# theta 3
			p4xz_norm = math.sqrt(p4z**2+(p4x-self.a[1])**2)
			# make sure p4xz_norm is in abs(a[2]-a[3]) and abs(a[2]+a[3])
			if abs(p4xz_norm - abs(self.a[2]-self.a[3])) < self.thr:
				t_3 = math.pi
			elif abs(p4xz_norm - abs(self.a[2]+self.a[3])) < self.thr:
				t_3 = 0
			else:
				t_3 = math.acos((p4xz_norm**2 - self.a[2]**2 - self.a[3]**2)/(2*self.a[2]*self.a[3]))
			
			for theta_3 in [t_3, -t_3]:
				try:
					# theta 2
					phi_3 = math.pi - theta_3
					phi_1 = math.atan2(p4z, (p4x-self.a[1]))
					phi_2 = math.asin(self.a[3]*math.sin(phi_3)/math.sqrt(p4z**2+(p4x-self.a[1])**2))
					theta_2 = phi_1 - phi_2

					# theta_4
					T_f3_r2 = self.T(3, theta_3)
					T_f2_r1 = self.T(2, theta_2)
					T_f4_r3 = np.matmul(self.inv_dh(T_f3_r2), self.inv_dh(T_f2_r1))
					T_f4_r3 = np.matmul(T_f4_r3, T_f4_r1)
					theta_4 = -math.atan2(T_f4_r3[0,1], T_f4_r3[0,0])
					rtn.append([theta_2, theta_3, theta_4])
				except Exception as ex:
					pass
		except Exception as ex:
			pass
		return rtn

	def adjust_degree(self, d):
		d = d%360
		return d - 360* (d > 180)

	def adjust_radian(self, r):
		r = r%(2*math.pi)
		return r - (2*math.pi)*(r>math.pi)

def main_random():
	thr = 0.001
	knmtc = dof_6()
	for i in range(100000):
		flag = True
		theta = [360*random.random()-180, 360*random.random()-180, 0, 0, 360*random.random()-180, 360*random.random()-180]
		_theta = [math.radians(t) for t in theta]
		dist_list = []
		T_f_tcp_r_base = knmtc.fw(_theta)
		theta_all = knmtc.inv(T_f_tcp_r_base)
		for __theta in theta_all:
			__theta = [knmtc.adjust_radian(t) for t in __theta]
			dist = np.linalg.norm(np.array(__theta) - np.array(_theta))
			dist_list.append(dist)
			if dist < thr:
				#print("dist ", dist, __theta, _theta)
				flag = False
				break
		if flag:
			print("theta ", i, theta)
			print("all ", i, [[math.degrees(t) for t in x] for x in theta_all])
			print("dist_list", dist, dist_list)
def main_diagnose():
	theta =  [-46.18777184724834, -124.78004640958426, 0, -140.31962984844893, -0.2814578757161428, 29.903238216683548]
	print(theta)
	_theta = [math.radians(t) for t in theta]
	
	knmtc = dof_6()
	T_f_tcp_r_base = knmtc.fw(_theta)
	T_f6_r0 = np.matmul(T_f_tcp_r_base, knmtc.inv_dh(knmtc.T_f_tcp_r6))
	T_f6_r0 = np.matmul(knmtc.inv_dh(knmtc.T_f0_r_base), T_f6_r0)

	theta_1 = knmtc.theta_1(T_f6_r0)
	print([math.degrees(t) for t in theta_1])
	
	theta_5 = knmtc.theta_5(T_f6_r0, theta_1[1])
	print([math.degrees(t) for t in theta_5])
	
	theta_6 = knmtc.theta_6(T_f6_r0, theta_1[1], theta_5[1])
	print(math.degrees(theta_6))

	theta_3_2_4 = knmtc.theta_3_2_4(T_f6_r0, theta_1[0], theta_5[1], theta_6)
	print([[math.degrees(t) for t in x] for x in theta_3_2_4])

if __name__ == '__main__':
	main_random()
	#main_diagnose()