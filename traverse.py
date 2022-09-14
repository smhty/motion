import numpy as np

"""
traverse from point S to E in a linear manner
"""
class linear(object):
	"""docstring for linear"""
	def __init__(self, point_s, point_e):
		super(linear, self).__init__()
		self.d = np.linalg.norm(point_e - point_s)
		self.a = np.array(point_s)
		self.b = (point_e - point_s) / self.d		
		

	def traverse(self, q):
		return self.a + q * self.b


"""
traverse from point S to M and E in a circular manner
"""
class circular(object):
	"""docstring for linear"""
	def __init__(self, point_s, point_m, point_e, turn = 0, dim = 3):

		super(circular, self).__init__()
		
		A = np.array(point_s[0:dim])
		B = np.array(point_e[0:dim])
		C = np.array(point_m[0:dim])
		
		a = A - C
		b = B - C

		# center of circle
		a_a = np.inner(a, a)
		b_b = np.inner(b, b)
		a_b = np.inner(a, b)

		p0 = (b_b * (a_a - a_b) * a) - (a_a * (a_b - b_b) * b)
		p0 = p0 / (2* (a_a * b_b - a_b**2))
		p0 = p0 + C   

		# radius
		r = np.linalg.norm(A-p0)

		# traverse vector
		r1 = A - p0
		r2 = B - p0
		r1_r2 = np.inner(r1, r2)
		r3 = r2 -  (r1* r1_r2/(r**2)) 
		r4 = r3 *(r/np.linalg.norm(r3))
		r5 = C - p0

		# find t_b and t_c
		cos_t_b = min(max(-1, r1_r2 / r**2), 1) 
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

		self.d = d
		self.r = r
		self.p0 = p0
		self.r1 = r1
		self.r2 = r4
		self.a = np.array(point_s[dim:])
		self.b = np.array(point_e[dim:] - point_s[dim:]) / d


	def traverse(self, q):
		t_0 = q / self.r
		xyz = self.p0 + math.cos(t_0)*self.r1 + math.sin(t_0)*self.r2
		ab = self.a + q * self.b 
		return np.append( xyz , ab )	
