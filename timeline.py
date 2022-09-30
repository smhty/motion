import profile
import traverse
import numpy as np
import math

"""
series of commands that form a continuous motion
commands are in abs
input
	point_data: {
		"cmd": "jmove" or "lmove"
		"joint":
		"xyz":
		"t":
		"a":
		"cont":
		"corner":
		"output": None or [0, 0, 1, 1, ]
	}
return
	plot_data:{
		
	}
"""
def points_to_plot(point_data, time_start=0, ticks_per_sec=1000, sample_per_sec=10):
	# adjust cont
	point_data[-1]["cont"] = 0
	point_data[0]["cont"] = 0

	# set the IK and position
	if point_data[0]["cmd"] == "lmove":
		ik = True
		for i in range(len(point_data)):
			point_data[i]["position"] = list(point_data[i]["xyz"])
	else:
		ik = False
		for i in range(len(point_data)):
			point_data[i]["position"] = list(point_data[i]["joint"])
	
	# init plot_data
	plot_data = []

	# init curve
	curve_start = None
	curve_end = None

	# init v_0
	v_0 = 0

	# loop over point_data
	for i in range(1, len(point_data)):
		# travel from A to B
		A = point_data[i-1]["position"]
		B = point_data[i]["position"]

		# form curve_end and v_e
		if point_data[i]["cont"] == 1:
			C = point_data[i+1]["position"]
			curve_end = traverse.curve(np.array(A), np.array(B), np.array(C), point_data[i]["corner"])
			curve_end.section(1)
			v_e = 0.001
		else:
			v_e = 0
			curve_end = None

		# traverse list line
		if curve_start:
			line_A = curve_start.end
		else:
			line_A = np.array(A)
		if curve_end:
			line_B = curve_end.start
		else:
			line_B = np.array(B)

		if all(np.isclose(line_A, line_B)):
			line = None
		else:
			line = traverse.linear(line_A, line_B)			


		# data = [{"time": time_abs, "avq": avq, "position": position}, ..., ]
		# prf = {"tick_jerk": t_j, "v_0": v_0, "v_e": v_e, "d": sum(d_list), "d_e": d_m+v_0*t_0+v_e*t_e, "a_avg": a_avg, "a_avg_e": _a_avg, "t": t, "t_m": t_m}
		prf, data = profile_to_plot([curve_start, line, curve_end], point_data[i]["t"]*ticks_per_sec, point_data[i]["a"], v_0, v_e, point_data[i]["j"], time_start, ticks_per_sec, sample_per_sec) 
		plot_data +=data

		# time_start
		time_start += prf["t_m"]/ticks_per_sec

		# adjust point_data
		point_data[i]["t_m"] = prf["t_m"]/ticks_per_sec
		point_data[i]["time_end"] = time_start
		point_data[i]["a_e"] = prf["a_avg_e"]

		#v_0
		v_0 = prf["v_e"]

		# from curve start
		if curve_end:
			curve_start = traverse.curve(np.array(A), np.array(B), np.array(C), point_data[i]["corner"])
			curve_start.section(2)
		else:
			curve_start = None

	# adjust first point
	point_data[0]["time_end"] = 0

	"""
	gen_cmd
	first run the motion and then run the IO
	give each command an ID starting from 1000
	"""
	cmds = []
	for point in point_data:
		# init cmd
		cmd = {k: point[k] for k in ["cmd", "t", "a", "cont", "corner"]}
		cmd["rel"] = 0
		cmd["id"] = 1000 + point["time_end"]
		if cmd["cmd"] == "jmove":
			for i in range(len(point["joint"])):
				cmd["j"+str(i)] = point["joint"][i]
		else:
			for i in range(len(point["xyz"])):
				cmd[["xyzabcde"][i]] = point["xyz"][i]
		cmds.append(cmd)
		
		# check if IO is available
		if "output" in point and point["output"]:
			cmd = {"out"+str(i):point["output"][i] for i in range(len(point["output"]))}
			cmd["id"] = 2000 + point["time_end"]
			cmd["cmd"] = "output"
			cmd["queue"] = 0
			cmds.append(cmd)

	return point_data, plot_data, cmds

"""
traverse object
	.d -> distance
	.traverse(q) -> return the value along the path after traversing distance q (from 0 to max distance) 

traverse list
trv_list: [start_curve, middle part, end_curve]
	trv_list[0] constant vel v_0 on start curve, 
	trv_list[1] middle part with accel, and decel
	trv_list[2] constant vel v_e on end curve 

t: total time
a_avg: input and output accel
v_0: initial v
v_e:ending v
j: jerk ratio from 0 to 1

time_start (seconds): is the time that tick 0 appears

we form a list, by sampling from the profile
each element has a time, vel, joint and xyz value
the sample points are the beginning of each jerk interval and also every frq_sample
"""
def profile_to_plot(trv_list, t, a_avg, v_0, v_e, j, time_start=0, ticks_per_sec=1000, sample_per_sec=10):
	# init data
	data = []

	# form d_list
	d_list = []
	for trv in trv_list:
		if trv == None:
			d_list.append(0)
		else:
			d_list.append(trv.d)

	# profile
	# {"tick_jerk": t_j, "v_0": v_0, "v_e": v_e, "d": sum(d_list), "d_e": d_m+v_0*t_0+v_e*t_e, "a_avg": a_avg, "a_avg_e": _a_avg, "t": t, "t_m": t_m}
	prf = profile.profile(d_list, t, a_avg, v_0, v_e, j)

	# initial avq
	q_init = 0 # 
	v_init = v_0
	a_init = 0

	for t_j in prf["tick_jerk"]:
		# length is not zero
		if t_j[0] == 0:
			continue
		
		# time_end
		time_end = time_start + (t_j[0]-1)/ticks_per_sec

		"""
		time_abs 
		"""
		time_abs_list = [k/sample_per_sec for k in range(math.ceil(time_start*sample_per_sec), math.floor(time_end*sample_per_sec)+1)]
		#print("time_end: ", time_end)
		if len(time_abs_list) == 0  or time_abs_list[-1] != time_end:
			time_abs_list.append(time_end)
		"""
		"""
		#print("time_abs_list: ", time_abs_list)
		# sample points
		# k/sample_per_sec >= time_start and k/sample_per_sec <= time_end
		#for k in range(math.ceil(time_start*sample_per_sec), math.floor(time_end*sample_per_sec)+1):
		for time_abs in time_abs_list:
			# time
			#time_abs = k/sample_per_sec	
			# tick
			tick =  math.floor((time_abs-time_start) * ticks_per_sec) 
			
			# avq
			avq = profile.tavq(tick+1, t_j[1], a0 = a_init, v0 = v_init, q0 = q_init) # accel, vel, q
				
			# position
			if avq[2] > d_list[0] + d_list[1]: 
				trv = trv_list[2]
			elif d_list[0] <= avq[2] <= d_list[0] + d_list[1]:
				trv = trv_list[1]
			else:
				trv = trv_list[0]

			# check if traverse exists
			if trv:
				position = trv.traverse(avq[2])
				data.append({"time": time_abs, "avq": avq, "position": position})

		# update final avq
		avq = profile.tavq(t_j[0], t_j[1], a0 = a_init, v0 = v_init, q0 = q_init)
		q_init = avq[2]
		v_init = avq[1]
		a_init = avq[0]

		# update time_start
		time_start = time_end + (1)/ticks_per_sec

	return prf, data


def main():
	import matplotlib.pyplot as plt
	from mpl_toolkits import mplot3d
	import time
	point_data = [
		{"cmd": "jmove", "joint":[0, 0], "xyz": [0], "t": 10, "a": 5e-9, "j": 0, "cont": 1, "corner": 100},
		{"cmd": "jmove", "joint":[1500, 1500], "xyz": [0], "t": 18, "a": 5e-9, "j": 0, "cont": 1, "corner": 100},
		{"cmd": "jmove", "joint":[3000, 0], "xyz": [0], "t": 16, "a": 5e-9, "j": 0, "cont": 1, "corner": 100},
		{"cmd": "jmove", "joint":[4500, 1500], "xyz": [0], "t": 14, "a": 5e-9, "j": 0, "cont": 1, "corner": 100},
		{"cmd": "jmove", "joint":[6000, 0], "xyz": [0], "t": 18, "a": 5e-9, "j": 0, "cont": 1, "corner": 100},
		{"cmd": "jmove", "joint":[7500, 1500], "xyz": [0], "t": 10, "a": 5e-9, "j": 0, "cont": 1, "corner": 100},
		{"cmd": "jmove", "joint":[9000, 0], "xyz": [0], "t": 12, "a": 5e-9, "j": 0, "cont": 1, "corner": 100},
	]

	point_data, plot_data, cmds = points_to_plot(point_data, time_start=0, ticks_per_sec=100000, sample_per_sec=5)	

	x = [d["time"] for d in plot_data]
	y = [d["avq"][1] for d in plot_data]

	plt.figure(1)
	plt.plot(x, y,'r-')
	plt.show()

if __name__ == '__main__':
	main()