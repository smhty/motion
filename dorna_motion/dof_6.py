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
def clamp(num, min_value, max_value):
        num = max(min(num, max_value), min_value)
        return num

class DH(object):
	"""docstring for dh"""
	def __init__(self):
		super(DH, self).__init__()

		self.alpha = [0, np.pi/2, 0, 0, np.pi/2, np.pi/2]
		self.a = [0, 1, 1, 1, 0,0,0]
		self.d = [0, 1, 0, 0,-1,1,1]

		self.T_f0_r_base = np.identity(4) # world
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
		self.thr = 0.0001

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
	def d_theta(self,t1,t2):
		dt = t2 - t1
		while dt > math.pi:
			dt = dt - 2 * math.pi
		while dt < -math.pi:
			dt = dt + 2 * math.pi
		return dt

	def angle_space_distance(self,s1 , s2):
		d = 0
		for i in range(len(s1)):
			d = d + self.d_theta(s1[i],s2[i])**2

		return d


	def inv(self, T_f_tcp_r_base, theta_current, all_sol):
		rtn = []
		T_f6_r0 = np.matmul(T_f_tcp_r_base, self.inv_dh(self.T_f_tcp_r6))
		T_f6_r0 = np.matmul(self.inv_dh(self.T_f0_r_base), T_f6_r0)

		for theta_1 in self.theta_1(T_f6_r0): # 2 x theta_1
			for theta_5 in self.theta_5(T_f6_r0, theta_1): # 2 x theta_5
				theta_6 = self.theta_6(T_f6_r0, theta_1, theta_5, theta_current[5])
				for theta_2_3_4 in self.theta_3_2_4(T_f6_r0, theta_1, theta_5, theta_6): # 1 x theta_2, # 2 x theta_3, # 1 x theta_4
					rtn.append([theta_1]+theta_2_3_4+[theta_5, theta_6])
		
		if all_sol:
			return rtn

		if theta_current and len(rtn)>0: 
			best_sol_dist = 1000
			best_sol_indx = 0
			indx = 0
			for sol in rtn:
				dist = self.angle_space_distance(sol,theta_current)
				if dist < best_sol_dist:
					best_sol_dist = dist
					best_sol_indx = indx
				indx  = indx + 1

			if best_sol_dist < 1.0:
				return [rtn[best_sol_indx]]
			else:
				return [theta_current]

		return rtn

	def theta_1(self, T_f6_r0, t1=None):
		p5_0 = np.matmul(T_f6_r0, [[0], [0], [-self.d[6]], [1]])
		p5x_0 = p5_0[0,0]
		p5y_0 = p5_0[1,0]
		rtn = []
		try:
			alpha = math.asin(clamp(-self.d[4]/math.sqrt(p5x_0**2 + p5y_0**2), -1.0 , 1.0))
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
			phi = math.acos(clamp(nom/self.d[6],-1.0,1.0))
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

		if abs(math.sin(theta_5)) < self.thr:
			return theta_6_init

		sgn_t_5 = 1.0
		if math.sin(theta_5)<0:
			sgn_t_5 = -1.0

		cos = -sgn_t_5*(-T_f6_r0[0,0]*math.sin(theta_1) + T_f6_r0[1,0]*math.cos(theta_1))#/(-math.sin(theta_5))
		sin = -sgn_t_5*(-T_f6_r0[0,1]*math.sin(theta_1) + T_f6_r0[1,1]*math.cos(theta_1))#/(-math.sin(theta_5))
		
		res = -math.atan2(sin, cos)

		return res

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
			#if abs(p4xz_norm - abs(self.a[2]-self.a[3])) < self.thr:
			#	t_3 = math.pi
			#elif abs(p4xz_norm - abs(self.a[2]+self.a[3])) < self.thr:
			#	t_3 = 0
			#elif p4xz_norm > min(abs(self.a[2]+self.a[3]), abs(self.a[2]-self.a[3])) and p4xz_norm < max(abs(self.a[2]+self.a[3]), abs(self.a[2]-self.a[3])):
			t_3 = math.acos(clamp((p4xz_norm**2 - self.a[2]**2 - self.a[3]**2)/(2*self.a[2]*self.a[3]),-1.0,1.0))
			#else:
			#	return rtn

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
					phi_2 = math.asin(clamp(self.a[3]*math.sin(phi_3)/math.sqrt(p4z**2+(p4x-self.a[1])**2),-1.0,1.0))
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
		
		self.sum_alpha = 3*np.pi/2
		self.delta = 0

		self.ssa = math.sin(self.sum_alpha)
		self.cd = math.cos(self.delta)
		self.sd = math.sin(self.delta)


	def rot_to_eul(self, rot):

		thr = 0.0001

		sgn = 1.0

		ssa = self.ssa
		cd = self.cd
		sd = self.sd

		sa = -rot[2,2]*ssa
		ca = sgn * math.sqrt(1.0 - sa*sa)

		A = math.atan2(sa,ca)
		B = 0.0
		G = 0.0


		if abs(ca)>thr:

			sb = ssa*rot[2,0]/ca
			cb = ssa*rot[2,1]/ca

			p1 = (rot[1,1]*ca+rot[1,2]*rot[2,1]*rot[2,2]/ca)

			sg = -1.0/rot[2,0]*cd*ssa*p1+rot[1,2]*sd/ca
			cg = -1.0/rot[2,0]*sd*p1 - rot[1,2]*cd*ssa/ca

			B = math.atan2(sb,cb)
			G = math.atan2(sg,cg)

		else:

			sb = - ssa*(rot[1,0]*cd*sa + rot[1,1]*sd)
			cb = ssa*(-rot[1,1]*cd*sa + rot[1,0]*sd)

			B = math.atan2(sb,cb)

		return [A, B, G]


	def eul_to_rot(self, ABG):

		ca = math.cos(ABG[0])
		sa = math.sin(ABG[0])

		cb = math.cos(ABG[1])
		sb = math.sin(ABG[1])

		cg = math.cos(ABG[2])
		sg = math.sin(ABG[2])

		ssa = self.ssa
		cd = self.cd
		sd = self.sd

		return np.matrix([
			[sa*sb*(cd*ssa*sg+cg*sd)+cb*(cg*cd-ssa*sg*sd),
			cg*(-cd*sb+cb*sa*sd)+ssa*sg*(cb*cd*sa+sb*sd),
			ca*(cd*ssa*sg+cg*sd)],
			[cg*ssa*(-cd*sa*sb+cb*sd)+sg*(cb*cd+sa*sb*sd),
			-sb*(cd*sg+cg*ssa*sd)+cb*sa*(-cg*cd*ssa+sg*sd),
			ca*(-cg*cd*ssa+sg*sd)],
			[ca*ssa*sb,
			ca*cb*ssa,
			-ssa*sa]])


class Dorna_c_knmtc(object):
	"""docstring for Dorna_c_knmtc"""
	def __init__(self):
		super(Dorna_c_knmtc, self).__init__()
		
		# create the 6 degree of freedom robot
		self.dof_6 = Dof_6()
		
		self.dof_6.a = [0, 100.0, 300.0, 208.5, 0.000, 0.000]
		self.dof_6.d = [0, 309.7, 0, 0, -133.1, 90.5, 9.707]
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

		return [fw[0,3]/1000, fw[1,3]/1000, fw[2,3]/1000] + abg


	def inv(self, xyzabg, joint_current=[], all_sol=True):
		ABG = [math.radians(t) for t in xyzabg[3:]]
		rot = self.euler.eul_to_rot(ABG)
		xyzabg[0] = 1000*xyzabg[0]
		xyzabg[1] = 1000*xyzabg[1]
		xyzabg[2] = 1000*xyzabg[2]

		T_f_tcp_r_base = np.matrix([
			[rot[0,0], rot[0,1], rot[0,2], xyzabg[0]],
			[rot[1,0], rot[1,1], rot[1,2], xyzabg[1]],
			[rot[2,0], rot[2,1], rot[2,2], xyzabg[2]],
			[0, 0, 0, 1]
		])

		# init condition
		theta_current = None
		if joint_current:
			theta_current = list(joint_current)
			if theta_current:
				theta_current = self.joint_to_theta(theta_current)
		theta_all = self.dof_6.inv(T_f_tcp_r_base, theta_current=theta_current, all_sol=all_sol)
		# all the solution
		joint_all = [self.theta_to_joint(theta) for theta in theta_all ]

		return joint_all

def main_random():
	thr = 0.001
	knmtc = Dof_6()
	knmtc.a = [0, 100.0, 300.0, 208.5, 0, 0]
	knmtc.d = [0, 309.7, 0, 0, -133.1, 90.5, 9.707]
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
				flag = False
				break
		if flag:
			print("theta ", i, theta)
			print("all ", i, [[math.degrees(t) for t in x] for x in theta_all])
			print("dist_list", dist, dist_list)


def main_diagnose():
	theta =  [0.0,0.0,0.0,90.0,0.0,0.0]
	_theta = [math.radians(t) for t in theta]
	
	knmtc = Dof_6()
	knmtc.a = [0, 100.0, 300.0, 208.5, 0.001, 0.001]
	knmtc.d = [0, 309.7, 0, 0, -133.1, 90.5, 9.707]
	T_f_tcp_r_base = knmtc.fw(_theta)

	T_f6_r0 = np.matmul(T_f_tcp_r_base, knmtc.inv_dh(knmtc.T_f_tcp_r6))
	T_f6_r0 = np.matmul(knmtc.inv_dh(knmtc.T_f0_r_base), T_f6_r0)


	theta_1 = knmtc.theta_1(T_f6_r0)
	print("t1:",[math.degrees(t) for t in theta_1])
	
	theta_5 = knmtc.theta_5(T_f6_r0, theta_1[0])
	print("t5:",[math.degrees(t) for t in theta_5])
	
	theta_6 = knmtc.theta_6(T_f6_r0, theta_1[0], theta_5[0])
	print("t6:",math.degrees(theta_6))

	theta_3_2_4 = knmtc.theta_3_2_4(T_f6_r0, theta_1[0], theta_5[0], theta_6)
	print("t3,2,4:",[[math.degrees(t) for t in x] for x in theta_3_2_4])


def main_dorna_c():
	thr = 0.001

	eul = Euler()

	print("hi")
	print(eul.rot_to_eul(np.matrix([[0,-1,0],[1,0,0],[0,0,1]])))
	print(eul.eul_to_rot(eul.rot_to_eul(np.matrix([[0,-1,0],[1,0,0],[0,0,1]]))))

	"""
	knmtc = Dorna_c_knmtc() 
	for i in range(1):
		flag = True
		joint = [-0.001, -0.001 ,-0.001 , -0.001 , -0.001 , -0.001]

		dist_list = []
		xyzabg = knmtc.fw(joint)
		print("xyzabg: ",xyzabg)
		joint_all = knmtc.inv(xyzabg, joint_current=joint, all_sol=False)
		print("final sol:",joint_all)
	"""
if __name__ == '__main__':
	#main_random()
	#main_diagnose()
	main_dorna_c()