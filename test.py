import matplotlib.pyplot as plt
import numpy as np

def p(t, a0 = 0.5, T = 10):
    return 1 + a0*t - (3 + 2*a0*T)*(t**2)/ (T**2) + (2+a0*T)*(t**3)/(T**3)


def main():
	a0 = -0.1
	T = 500
	t1 = np.arange(0, T, 1)

	print()
	plt.figure()
	plt.subplot(111)
	plt.plot(t1, p(t1, a0, T), 'bo')	
	plt.show()



if __name__ == '__main__':
	main()
