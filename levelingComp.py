
BS = 1.25

IS = [1.68, 1.72, 1.82]

listIS = []

listIS.append(BS)
listIS.extend(IS)

rise = []
fall = []

RL_const = 10.00

RL_LIST = []

for i,x in enumerate(listIS):
	try:
		asd = listIS[i: i+2]
		ans = asd[0] - asd[1]
		if (ans > 0):
			rise.append(ans)
		else:
			fall.append(abs(ans))
	except IndexError:
		print 'done'


print "RISE: ", rise

print "FALL: ", fall
ans = RL_const
for f,xx in enumerate(fall):
	hh = ans - xx
	ans = hh
	print hh
