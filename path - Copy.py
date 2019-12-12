"""
p1 is the start
p2 is the middle point
p3 is the end point


"""
import numpy as np
import matplotlib.pyplot as plt


def path_2_prm(p_1, p_2, p_3, r):
	# vector
	v_1 = np.subtract(p_1, p_2)
	v_2 = np.subtract(p_3, p_2)

	# the angle forms by p_1, p_2, p_3
	theta = np.arccos(np.inner(v_1, v_2)/(np.linalg.norm(v_1) * np.linalg.norm(v_2)))
	
	# point where v_1 touches the circle
	A = p_2 + r/np.tan(theta/2)*v_1/np.linalg.norm(v_1) 
	
	# point where v_2 touches the circle
	B = p_2 + r/np.tan(theta/2)*v_2/np.linalg.norm(v_2) 	

	# center of the circle
	C = p_2 + (r/np.sin(theta)) *(v_1/np.linalg.norm(v_1) + v_2/np.linalg.norm(v_2))

	"""
	distance
	distance covers by v_1
	distance covers by circle
	distance covers by v_2 
	"""
	d = [np.linalg.norm(v_1 - A),
		r*(np.pi - theta),
		np.linalg.norm(v_2 - B)]
	
	return [v_1, v_2, A, B, C, d] 

if __name__ == '__main__':
	
	p_3 = [10,0]
	p_2 = [5,4]
	p_1 = [22,12]
	r = 3
	points = 100

	[v_1, v_2, A, B, C, d] = path_2_prm(p_1, p_2, p_3, r)

	# make data
	x_1 = []
	y_1 = []

	x_2 = []
	y_2 = []

	x_3 = []
	y_3 = []	

	x_4 = [A[0], B[0], C[0], p_1[0], p_2[0], p_3[0]]
	y_4 = [A[1], B[1], C[1], p_1[1], p_2[1], p_3[1]]

	print(C, (A-C),  - r*v_1/np.linalg.norm(v_1),1/r) 
	for i in range(points):
		l = sum(d)* i/points
		if l <= d[0]:
			p = p_1 + l/d[0]*(A-p_1)
			x_1.append(p[0]) 
			y_1.append(p[1]) 
		elif l <= d[0] + d[1]:
			phi = (l-d[0])/r
			p = C + (A-C)*np.cos(phi) - r*np.sin(phi)*v_1/np.linalg.norm(v_1)
			x_2.append(p[0]) 
			y_2.append(p[1])
		else:
			_d = l - d[0] - d[1]
			p = B + _d/d[2] *(p_3-B)  
			#print(_d/d[2]) 
			x_3.append(p[0]) 
			y_3.append(p[1]) 

	# Create plots with pre-defined labels.
	fig, ax = plt.subplots()
	ax.plot(x_1, y_1, 'bo-', label='in')
	ax.plot(x_2, y_2, 'ro-', label='circle')
	ax.plot(x_3, y_3, 'go-', label='out')
	ax.plot(x_4, y_4, 'ko', label='points')
	

	legend = ax.legend(loc='upper center', shadow=True, fontsize='x-large')

	# Put a nicer background color on the legend.

	plt.show() 	
