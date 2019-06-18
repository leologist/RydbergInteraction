import itertools
import numpy as np

def create_plaq(x, y):
	return [[x, y, 'a'], [x+1, y, 'p'], [x, y+1, 'p'], [x-1, y, 'p'], [x, y-1, 'p']]

Q = np.zeros((29, 29))

index = 0;
finallist = []
for i, j in [(0, 0), (-1, -1), (1, -1), (-2, -2), (0, -2), (2, -2), (-3, -3), (-1, -3), (1, -3), (3, -3)]:
	for itm1 in create_plaq(i, j):
		for itm2 in create_plaq(i, j):
			if(itm1!=itm2):
				x = finallist.index(itm1) if itm1 in finallist else None
				y = finallist.index(itm2) if itm2 in finallist else None
				if(x is None):
					finallist.append(itm1)
					x = index
					index+=1

					
				if(y is None):
					finallist.append(itm2)
					y = index
					index+=1
				if(itm1[2] == 'a' or itm2[2] == 'a'):
					Q[x][y] += 2.0
					Q[y][x] += 2.0
				else:
					Q[x][y] += 1.0
					Q[y][x] += 1.0
print(Q)

tosort = np.zeros((2**29, 2))

index = 0

for i in itertools.product([-1, 1], repeat=29):
	tosort[index][0] = np.dot(i,np.dot(Q, i))
	tosort[index][1] = abs(i[4] + i[7] + i[10] + i[15])
	index+=1
	if index % 1000000 == 0:
		print(index)
		break

np.save("tosort.npy", tosort)