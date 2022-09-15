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
class DH(object):
	"""docstring for dh"""
	def __init__(self):
		super(DH, self).__init__()

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

class Dof_6(DH):
	"""docstring for dof_6"""
	def __init__(self):
		super(Dof_6, self).__init__()
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
	def inv(self, T_f_tcp_r_base, theta_current, all_sol):
		# init return
		rtn = []

		T_f6_r0 = np.matmul(T_f_tcp_r_base, self.inv_dh(self.T_f_tcp_r6))
		T_f6_r0 = np.matmul(self.inv_dh(self.T_f0_r_base), T_f6_r0)
		if all_sol:
			for theta_1 in self.theta_1(T_f6_r0): # 2 x theta_1
				for theta_5 in self.theta_5(T_f6_r0, theta_1): # 2 x theta_5
					if theta_current:
						theta_6 = self.theta_6(T_f6_r0, theta_1, theta_5, theta_6_init=theta_current[5]) # 1 x theta_6
					else:
						theta_6 = self.theta_6(T_f6_r0, theta_1, theta_5)
					for theta_2_3_4 in self.theta_3_2_4(T_f6_r0, theta_1, theta_5, theta_6): # 1 x theta_2, # 2 x theta_3, # 1 x theta_4
						rtn.append([theta_1]+theta_2_3_4+[theta_5, theta_6])
		
		elif not all_sol and theta_current:
			for theta_1 in self.theta_1(T_f6_r0, t1=theta_current[0]): # 2 x theta_1
				for theta_5 in self.theta_5(T_f6_r0, theta_1, t5=theta_current[4]): # 2 x theta_5
					theta_6 = self.theta_6(T_f6_r0, theta_1, theta_5, theta_6_init=theta_current[5]) # 1 x theta_6
					for theta_2_3_4 in self.theta_3_2_4(T_f6_r0, theta_1, theta_5, theta_6, t3=theta_current[2]): # 1 x theta_2, # 2 x theta_3, # 1 x theta_4
						rtn.append([theta_1]+theta_2_3_4+[theta_5, theta_6])

		return rtn

	def theta_1(self, T_f6_r0, t1=None):
		p5_0 = np.matmul(T_f6_r0, [[0], [0], [-self.d[6]], [1]])
		p5x_0 = p5_0[0,0]
		p5y_0 = p5_0[1,0]
		rtn = []
		try:
			alpha = math.asin(-self.d[4]/math.sqrt(p5x_0**2 + p5y_0**2))
			phi_1 = math.atan2(p5y_0, p5x_0)

			rtn = [phi_1-alpha, phi_1+alpha-math.pi]
			if t1 != None:
				if min(abs(t1-rtn[0]%(2*math.pi)), abs(t1-rtn[0]%(-2*math.pi))) > min(abs(t1-rtn[1]%(2*math.pi)), abs(t1-rtn[1]%(-2*math.pi))):
					rtn.pop(0)
				else:
					rtn.pop(1)
			return rtn
		except Exception as ex:
			return rtn

	def theta_5(self, T_f6_r0, theta_1, t5=None):
		nom = T_f6_r0[1,3]*np.cos(theta_1)-T_f6_r0[0,3]*np.sin(theta_1)+self.d[4]
		rtn = []
		try:
			phi = math.acos(nom/self.d[6])
			rtn = [phi, -phi]
			
			if t5 !=None:
				if t5*rtn[0] >= 0:
					rtn.pop(1)
				else:
					rtn.pop(0)

			return rtn
		except Exception as ex:
			return []

	def theta_6(self, T_f6_r0, theta_1, theta_5, theta_6_init=0):
		sgn = 1
		if abs(math.sin(theta_5)) < self.thr:
			return theta_6_init
		elif math.sin(theta_5) < 0:
			sgn = -1
		cos = (-T_f6_r0[0,0]*math.sin(theta_1) + T_f6_r0[1,0]*math.cos(theta_1))/(-sgn)
		sin = (-T_f6_r0[0,1]*math.sin(theta_1) + T_f6_r0[1,1]*math.cos(theta_1))/(-sgn)
		return -math.atan2(sin, cos)

	def theta_3_2_4(self, T_f6_r0, theta_1, theta_5, theta_6, t3=None):
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
			elif p4xz_norm > min(abs(self.a[2]+self.a[3]), abs(self.a[2]-self.a[3])) and p4xz_norm < max(abs(self.a[2]+self.a[3]), abs(self.a[2]-self.a[3])):
				t_3 = math.acos((p4xz_norm**2 - self.a[2]**2 - self.a[3]**2)/(2*self.a[2]*self.a[3]))
			else:
				return rtn

			t_3_list = [t_3, -t_3]
			
			if t3 !=None:
				if t3*t_3_list[0] >= 0:
					t_3_list.pop(1)
				else:
					t_3_list.pop(0)
			for theta_3 in t_3_list:
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


"""
Euler ZYX mobile
alpha: around z
beta: around mobile y
gamma: around mobile x
"""
class Euler(object):
	"""docstring for Euler"""
	def __init__(self):
		super(Euler, self).__init__()
		

	def rot_to_eul(self, rot):
		# B = math.pi/2
		if max(abs(rot[0, 0]), abs(rot[1, 0])) < 0.000001:
			B = math.pi/2
			A = 0
			G = math.atan2(rot[0, 1], rot[1, 1])
		
		else:
			B = math.atan2(-rot[2, 0], math.sqrt(rot[0, 0]**2 + rot[1, 0]**2))
			A = math.atan2(rot[1, 0], rot[0, 0])
			G = math.atan2(rot[2, 1], rot[2, 2])

		return [A, B, G]


	def eul_to_rot(self, ABG):
		ca = math.cos(ABG[0])
		sa = math.sin(ABG[0])

		cb = math.cos(ABG[1])
		sb = math.sin(ABG[1])

		cg = math.cos(ABG[2])
		sg = math.sin(ABG[2])

		return np.matrix([
			[ca*cb, ca*sb*sg-sa*cg, ca*sb*cg+sa*sg],
			[sa*cb, sa*sb*sg+ca*cg, sa*sb*cg-ca*sg],
			[-sb, cb*sg, cb*cg]])


class C_kntmc(object):
	"""docstring for dorna_c"""
	def __init__(self):
		super(C_kntmc, self).__init__()
		
		# create the 6 degree of freedom robot
		self.dof_6 = Dof_6()

		# create Euler
		self.euler = Euler()

	def joint_to_theta(self, joint):
		theta = list(joint)
		theta[3] +=90 
		return [math.radians(j) for j in theta]


	def theta_to_joint(self, theta):
		joint = [math.degrees(t) for t in theta]
		joint[3] -=90
		return [self.dof_6.adjust_degree(j) for j in joint]

	def fw(self, joint):
		# adjust theta to dof_6
		_theta = self.joint_to_theta(joint)
		
		# fw result
		fw = self.dof_6.fw(_theta)

		# abg
		abg = self.euler.rot_to_eul(fw[0:3, 0:3])
		abg = [math.degrees(r) for r in abg]

		return [fw[0,3], fw[1,3], fw[2,3]] + abg


	def inv(self, xyzabg, joint_current=[], all_sol=True):
		ABG = [math.radians(t) for t in xyzabg[3:]]
		rot = self.euler.eul_to_rot(ABG)

		T_f_tcp_r_base = np.matrix([
			[rot[0,0], rot[0,1], rot[0,2], xyzabg[0]],
			[rot[1,0], rot[1,1], rot[1,2], xyzabg[1]],
			[rot[2,0], rot[2,1], rot[2,2], xyzabg[2]],
			[0, 0, 0, 1]
		])

		# init condition
		theta_current = list(joint_current)
		if theta_current:
			theta_current = self.joint_to_theta(theta_current)
		
		theta_all = self.dof_6.inv(T_f_tcp_r_base, theta_current=theta_current, all_sol=all_sol)
		
		# all the solution
		joint_all = [self.theta_to_joint(theta) for theta in theta_all ]

		# adjust j[5] for infinite rotation
		if joint_current:
			rnd = math.floor(joint_current[5]/360)
			if abs(joint_current[5] - (rnd+1) *360) < abs(joint_current[5] - rnd *360):
				rnd += 1 

			for j in joint_all:
				j[5] += rnd*360
		
		return joint_all

def main_random():
	thr = 0.001
	knmtc = Dof_6()
	for i in range(1):
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
	
	knmtc = Dof_6()
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


def main_dorna_c():
	thr = 0.001
	knmtc = C_kntmc()
	for i in range(100000):
		flag = True
		joint = [360*random.random()-180, 360*random.random()-180, 360*random.random()-180, 360*random.random()-180, 360*random.random()-180, -720*random.random()-360]
		#joint = [110.61048010731821, -170.60639070054873, 0, 0, -56.78812718139007, -125.48839700234943]
		dist_list = []
		xyzabg = knmtc.fw(joint)
		joint_all = knmtc.inv(xyzabg, joint_current=joint, all_sol=False)
		dist = np.linalg.norm(np.array(joint) - np.array(joint_all[0]))
		if dist > 0.001:
			print(i, joint)
			print(joint_all)
			print("######")			
			pass
		#print(i, joint)
		#print(joint_all)
		#print("######")
if __name__ == '__main__':
	#main_random()
	#main_diagnose()
	main_dorna_c()