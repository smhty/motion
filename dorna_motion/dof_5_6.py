import numpy as np
import time
import math
import random
from .cf import *   #Relative import needed to build package
#from cf import CF  #Import for locally testing the file
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

def d_theta(t1,t2):
	dt = t2 - t1
	while dt > math.pi:
		dt = dt - 2 * math.pi
	while dt < -math.pi:
		dt = dt + 2 * math.pi
	return dt

def angle_space_distance(s1 , s2):
	d = 0
	for i in range(len(s1)):
		d = d + d_theta(s1[i],s2[i])**2
	return d


class DH(object):
	"""docstring for dh"""
	def __init__(self):
		super(DH, self).__init__()

		self.n_dof = 6 # number of degrees of freedom, choose between [5,6]
		self.alpha = [0, 0, np.pi/2, 0, 0, 0, 0] # Rotation of Ci with respect to C(i-1) around the x axis of Ci
		self.delta = [0, 0, 0, 0, 0, np.pi/2, np.pi/2] #Rotation of Ci with respect to C(i-1) around the y axis of Ci
		self.a = [0, 1, 1, 1, 0,0,0] # Translation of Ci with respect to C(i-1) along the x axis of Ci
		self.d = [0, 1, 0, 0,-1,1,1] #Translation of Ci with respect to C(i-1) along the z axis of Ci
		self.base_tr_vector = [1,0,0] # the vector that describes the motion of base on rail

		self.cf_test = CF(ndof = self.n_dof) 

		self.T_f0_r_base = np.identity(4) # world
		self.T_f_tcp_r_last = np.identity(4) # TCP

	# Ti with respect to i-1 frame 
	def T(self, i, theta):
		ct = np.cos(theta)
		st = np.sin(theta)
		ca = np.cos(self.alpha[i])
		sa = np.sin(self.alpha[i])
		cd = np.cos(self.delta[i])
		sd = np.sin(self.delta[i])

		result =  np.matrix([
			[cd*ct, -cd*st, sd, 	self.a[i]*cd + self.d[i]*sd],
			[ct*sa*sd+ca*st, ca*ct-sa*sd*st, -cd*sa,	self.a[i]*(sa*sd) + self.d[i]*(-cd*sa)],
			[-ca*ct*sd + sa*st, ct*sa+ca*sd*st, cd*ca,	self.a[i]*(-ca*sd ) + self.d[i]*(cd*ca)],
			[0, 0, 0, 1]])

		return result

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

class Dof(DH):
	"""docstring for dof"""
	def __init__(self):
		super(Dof, self).__init__()
		self.thr = 0.0001

	"""
	joint values are given in radians
	return T_f_tcp_r_base
	"""
	def fw(self, theta):			
		T = self.T_f0_r_base
		for i in range(1,self.n_dof+1 ):
			T = np.matmul(T, self.T(i, theta[i-1]))
		
		return np.matmul(T, self.T_f_tcp_r_last)

	"""
	The robot T_f_tcp_r_base is given
	find all the possible robot orientations 
	"""



	def inv(self, T_f_tcp_r_base, theta_current, all_sol):
		rtn = []
		T_f_last_r0 = np.matmul(T_f_tcp_r_base, self.inv_dh(self.T_f_tcp_r_last)) 
		T_f_last_r0 = np.matmul(self.inv_dh(self.T_f0_r_base), T_f_last_r0)

		self.cf_test.set_matrix(T_f_last_r0)
		abc = self.cf_test.get_euler()

		for theta_1 in self.theta_1(T_f_last_r0): # 2 x theta_1
			for theta_5 in self.theta_5(T_f_last_r0, theta_1 , abc = abc): # 2 x theta_5
				theta_6 = self.theta_6(T_f_last_r0, theta_1, theta_5, theta_current[5])
				for theta_2_3_4 in self.theta_3_2_4(T_f_last_r0, theta_1, theta_5, theta_6): # 1 x theta_2, # 2 x theta_3, # 1 x theta_4
					rtn.append([theta_1]+theta_2_3_4+[theta_5, theta_6])
		
		if all_sol:
			return rtn

		if theta_current and len(rtn)>0: 
			best_sol_dist = 1000
			best_sol_indx = 0
			indx = 0
			for sol in rtn:
				dist = angle_space_distance(sol,theta_current)
				if dist < best_sol_dist:
					best_sol_dist = dist
					best_sol_indx = indx
				indx  = indx + 1

			if best_sol_dist < 10.0:
				return [rtn[best_sol_indx]]
			else:
				print("bad")
				return [theta_current]

		return rtn

	def theta_1(self, T_f_last_r0, t1=None):
		if(self.n_dof == 6):
			p5_0 = np.matmul(T_f_last_r0, [[0], [0], [-self.d[6]], [1]])
		else:
			p5_0 = np.matmul(T_f_last_r0, [[0], [0], [0], [1]])

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

	def theta_5(self, T_f_last_r0, theta_1, t5=None , abc=[0,0,0]):
		if self.n_dof == 5:
			return [abc[1]]

		nom = T_f_last_r0[1,3]*np.cos(theta_1)-T_f_last_r0[0,3]*np.sin(theta_1)+self.d[4]
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


	def theta_6(self, T_f_last_r0, theta_1, theta_5, theta_6_init=0):

		if(self.n_dof==5):
			return 0

		if abs(math.sin(theta_5)) < self.thr:
			return theta_6_init

		sgn_t_5 = 1.0
		if math.sin(theta_5)<0:
			sgn_t_5 = -1.0

		cos = sgn_t_5*(T_f_last_r0[0,1]*math.sin(theta_1) - T_f_last_r0[1,1]*math.cos(theta_1))#/(-math.sin(theta_5))
		sin = sgn_t_5*(T_f_last_r0[0,0]*math.sin(theta_1) - T_f_last_r0[1,0]*math.cos(theta_1))#/(-math.sin(theta_5))
		
		res = math.atan2(sin, cos)

		return res

	def theta_3_2_4(self, T_f_last_r0, theta_1, theta_5, theta_6, t3=None):
		rtn = []
		T_f1_r0 = self.T(1, theta_1)
		T_f5_r4 = self.T(5, theta_5)


		T_f4_r1 = np.matmul(self.inv_dh(T_f1_r0), T_f_last_r0)
		if(self.n_dof == 6):
			T_f6_r5 = self.T(6, theta_6)
			T_f4_r1 = np.matmul(T_f4_r1, self.inv_dh(T_f6_r5))
		T_f4_r1 = np.matmul(T_f4_r1, self.inv_dh(T_f5_r4))

		p4x = T_f4_r1[0,3]
		p4z = T_f4_r1[2,3]

		try:
			# theta 3
			p4xz_norm = math.sqrt(p4z**2+(p4x-self.a[2])**2)
			t_3 = math.acos(clamp((p4xz_norm**2 - self.a[3]**2 - self.a[4]**2)/(2*self.a[3]*self.a[4]),-1.0,1.0))

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
					phi_1 = math.atan2(p4z, (p4x-self.a[2]))
					phi_2 = math.asin(clamp(self.a[4]*math.sin(phi_3)/math.sqrt(p4z**2+(p4x-self.a[2])**2),-1.0,1.0))
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

class Dorna_c_knmtc(object):
	"""docstring for Dorna_c_knmtc"""
	def __init__(self):
		super(Dorna_c_knmtc, self).__init__()
		
		# create the 6 degree of freedom robot
		self.dof = Dof()
		
		#self.dof.a = [0, 0 , 95.4333, 203.1997, 152.398, 0, 0]#[0, 0 , 100.0, 300.0, 208.5, 0, 0]
		#self.dof.d = [0, 218.529, 0, 0, 1, 13.275, 1]#[0, 309.7, 0, 0,  -133.1, 90.5, 9.707]


	def joint_to_theta(self, joint):
		theta = list(joint)
		return [math.radians(j) for j in theta]


	def theta_to_joint(self, theta):
		joint = [math.degrees(t) for t in theta]
		return [self.dof.adjust_degree(j) for j in joint]

	def fw(self, joint):
		# adjust theta to dof
		_theta = self.joint_to_theta(joint)
		
		# fw result
		fw = self.dof.fw(_theta)
		# abg

		self.dof.cf_test.set_matrix(fw)
		abc = self.dof.cf_test.get_euler()
		abc = [math.degrees(r) for r in abc]

		#give different result: fw, fw can later be changed to pos + abg
		return [fw[0,3]/1000, fw[1,3]/1000, fw[2,3]/1000] + abc


	def inv(self, xyzabc, joint_current=[0,0,0,0,0,0], all_sol=True): #xyzabg
		#print("inv call:",xyzabc)
		ABC = [math.radians(t) for t in xyzabc[3:]]
		
		self.dof.cf_test.set_euler(ABC)
		rot = self.dof.cf_test.local_matrix 

		xyzabc[0] = 1000*xyzabc[0]
		xyzabc[1] = 1000*xyzabc[1]
		xyzabc[2] = 1000*xyzabc[2]

		T_f_tcp_r_base = np.matrix([
			[rot[0,0], rot[0,1], rot[0,2], xyzabc[0]],
			[rot[1,0], rot[1,1], rot[1,2], xyzabc[1]],
			[rot[2,0], rot[2,1], rot[2,2], xyzabc[2]],
			[0, 0, 0, 1]
		])

		# init condition
		theta_current = None
		if joint_current:
			theta_current = list(joint_current)
			if theta_current:
				theta_current = self.joint_to_theta(theta_current)
		theta_all = self.dof.inv(T_f_tcp_r_base, theta_current=theta_current, all_sol=all_sol)
		# all the solution
		joint_all = [self.theta_to_joint(theta) for theta in theta_all ]

		#print("resulting xyzabc:",self.fw(joint_all[0]))

		return joint_all

def main_dorna_c():
	thr = 0.001

	
	knmtc = Dorna_c_knmtc()

	for i in range(1):
		
		flag = True
		
		joint = [5,10,20,30,40,0]#[0,44,-127,83,0,0]
		print("in:",joint)
		dist_list = []
		xyzabg = knmtc.fw(joint)
		print("xyzabg: ",xyzabg)
		joint_all = knmtc.inv(xyzabg, joint_current=joint, all_sol=True)
		print("final sol:",joint_all)



if __name__ == '__main__':
	#main_random()
	#main_diagnose()
	main_dorna_c()