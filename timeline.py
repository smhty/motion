import profile
"""
series of commands are given to the robot and break it to groups
each group is a set of points that follow a continuous motion
"""
class Parse(object):
	"""docstring for waypoint"""
	def __init__(self, arg):
		super(waypoint, self).__init__()
		self.current = []
		self.grps = []
	

	"""
	convert commands to groups
	"""
	def cmds_to_groups(self, cmds):
		grps = []
		while len(cmds) != 0
			# init group
			grp = []
			while len(cmds) != 0:
				if len(grp) == 0:
					grp.append(cmds.pop(0))

					# circle
					if grp[0]["cmd"] == "cmove":
						break

					# cont == 0 or no cont
					if "cont" not in grp[0] or grp[0]["cont"] == 0:
						break

				else:
					# different cmd
					if cmds[0]["cmd"] != grp[-1]["cmd"]:
						break

					# add to group
					grp.append(cmds.pop(0))

					# cont == 0
					if "cont" in grp[-1] and grp[-1]["cont"] == 0:
						break
			
			# check for the size
			if len(grp):
				grps.append(grp)

		self.grps = grps


	def group_to_plot()

	def plot_data(self, group_index=0):


"""
Each group is a series of continuous commands
"""
class Group(object):
	"""docstring for Group"""
	def __init__(self, cmds):
		super(Group, self).__init__()
		self.cmds = cmds

	def plot_data_create(cmds):


	def plot_data_update

"""
traverse object
	.distance -> distance
	.q(d) -> return the value along the path after traversing distance d (from 0 to max distance) 

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
def profile_plot_data(trv_list, t, a_avg, v_0, v_e, j, time_start=0, ticks_per_sec=1000, sample_per_sec=100):
	# init data
	data = []

	# profile
	# {"tick_jerk": t_j, "v_0": v_0, "v_e": v_e, "d": sum(d_list), "d_e": d_m+v_0*t_0+v_e*t_e, "a_avg": a_avg, "a_avg_e": _a_avg, "t": t, "t_m": t_m}
	prf = profile.profile([trv.d for trv in trv_list], t, a_avg, v_0, v_e, j)
	
	# initial avq
	q_1 = 0 # 
	v_1 = v_0
	a_1 = 0

	for t_j in prf["tick_jerk"]:
		# length is not zero
		if t_j[0] == 0:
			continue
		
		# time_end
		time_end = time_start + (t_j[0]-1)/ticks_per_sec

		# sample points
		# k/sample_sec >= time_start and k/sample_sec <= time_end
		for k in range(math.ceil(time_start/sample_per_sec), math.floor(time_end/sample_per_sec)+1)
			# time to tick
			tick =  math.floor(k*ticks_per_sec) # tick
			
			# avq
			avq = profile.tavq(tick, t_j[1], a0 = a_1, v0 = v_1, q0 = q_1) # accel, vel, q
			
			# position 
			position = trv.q(avq[2]) # position
			
			data.append({"tick": tick, "avq": avq, "position": position})

		# update final avq
		avq = profile.tavq(t_j[0]-1, t_j[1], a0 = a_1, v0 = v_1, q0 = q_1)
		q_1 = avq[2]
		v_1 = avq[1]
		a_1 = avq[0]

		# update time_start
		time_start = time_end + (1)/ticks_per_sec

	return data