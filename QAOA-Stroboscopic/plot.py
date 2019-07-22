import matplotlib.pyplot as plt

data_arr = [[1, 6, 9.202439046895883e-09, 3.295932529668139e-09, 42.0, 2.0],
[2, 6, 9.73979369378597e-11, 1.963814960070796e-10, 130.5, 64.60456640207408],
[3, 6, 1.9005354096046645e-09, 1.7104121857591144e-09, 232.0, 8.0],
[4, 6, 3.1413711909511477e-09, 2.6066976266763884e-09, 352.5, 148.13422967025548],
[1, 8, 0.18816457320336366, 0.18816457320335903, 26.0, 6.0],
[2, 8, 5.085792897929764e-13, 4.869084597684386e-13, 126.0, 0.0],
[3, 8, 2.9280466939951566e-12, 2.836930187255757e-12, 172.0, 52.0], 
[4, 8, 1.6343805891994023e-09, 3.8910018992207236e-09, 405.0, 46.636895265444075],
[1, 10, 0.1153054270799852, 0.1153054270798704, 28.0, 8.0],
[2, 10, 1.2340128918708615e-13, 1.1618483952702263e-13, 111.0, 27.0],
[3, 10, 0.07107603048471955, 0.07107603048459732, 228.0, 108.0], 
[4, 10, 0.08924869687660188, 0.08924869662691905, 446.25,86.88174434252572],
[1, 12, 3.4861002973229915e-14, 2.942091015256665e-14, 22.0, 2.0]]

X_accuracy = [[] for _ in range(5)]
Y_accuracy = [[] for _ in range(5)]
E_accuracy = [[] for _ in range(5)]

X_time = [[] for _ in range(5)]
Y_time = [[] for _ in range(5)]
E_time = [[] for _ in range(5)]

for data in data_arr: 
	index = int((data[1] / 2)-3)
	X_accuracy[index].append(data[0])
	Y_accuracy[index].append(data[2])
	E_accuracy[index].append(data[3])

	X_time[index].append(data[0])
	Y_time[index].append(data[4])
	E_time[index].append(data[5])

for index in range(5):
	plt.errorbar(X_accuracy[index], Y_accuracy[index], yerr = E_accuracy[index], capsize=2, elinewidth=1, fmt='.--', label='{t} vertices'.format(t=(index+3)*2))

plt.xlabel("p")
plt.ylabel("1-r")
plt.legend(loc='best')
plt.show()

import matplotlib.pyplot as plt
for index in range(5):
	plt.errorbar(X_time[index], Y_time[index], yerr = E_time[index], capsize=2, elinewidth=1, fmt='.--', label='{t} vertices'.format(t=(index+3)*2))

plt.xlabel("p")
plt.ylabel("optimization time")
plt.legend(loc='best')
plt.yscale("log")
plt.show()

