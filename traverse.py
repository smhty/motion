import numpy as np
import math

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


"""
curve
"""
class curve(object):
	"""docstring for curve"""
	def __init__(self, point_s, point_m, point_e, corner, hlf=0):
		super(curve, self).__init__()
		BA = point_s - point_m
		BC = point_e - point_m

		LBA = np.linalg.norm(BA)
		LBC = np.linalg.norm(BC)

		BA = BA/LBA
		BC = BC/LBC

		theta = 0.5*(math.pi - math.acos(np.inner(BA, BC)))
		alpha = math.tan(theta)

		y = BA + BC
		y = y/np.linalg.norm(y)

		x = BA - np.inner(BA, y)* y
		x = x/np.linalg.norm(x)

		L = min(corner, 0.5*min(LBA, LBC))

		self.base = point_m # edge position
		self.x = x # curve plane x coordinate
		self.y = y # curve plane y coordinate
		self.l = L # Length we chop from line
		self.alpha = alpha # Tangent value of edge angel min 0 max PI/2
		self.theta = theta # the angle between final line and initial line divided by two
		self.d = 2*L # Integrated length of curve

		# helper variables
		self._L3 = L**3
		self._L4 = L**4
		self._alpha = 1/math.sqrt(1.0 + self.alpha**2)
		self.__alpha = self.alpha / (8.0*math.sqrt(1.0 + self.alpha**2)*self._L3)

		self.hlf = hlf

		self.start = self.traverse(0)
		self.middle = self.traverse(L)
		self.end = self.traverse(2*L)

	def section(self, hlf):
		self.hlf = hlf
		if hlf != 0:
			self.d = self.l

	def traverse(self, q):
		if self.hlf == 2: # traverse second half
			q = self.d+q

		x = (self.l - q) * self._alpha;
		y = self.__alpha * (8.0 *self._L4 - 8.0*self._L3*q + 4.0*self.l*(q**3) - (q**4));

		return self.base + x*self.x + y*self.y


def main():
	import matplotlib.pyplot as plt
	point_s = np.array([0, 0])
	point_m = np.array([0, 3])
	point_e = np.array([0, 6])
	corner = 3.5
	counter = 100

	crv = curve(point_s, point_m, point_e, corner)

	x = []
	y = []
	x.append(point_s[0])
	y.append(point_s[1])

	for i in range(counter):
		point = crv.traverse(crv.d*(i/counter))
		x.append(point[0])
		y.append(point[1])

	x.append(point_e[0])
	y.append(point_e[1])

	x.append(point_m[0])
	y.append(point_m[1])

	x.append(point_s[0])
	y.append(point_s[1])

	plt.figure(1)
	plt.plot(x, y,'r-')
	plt.show()

if __name__ == '__main__':
	main()