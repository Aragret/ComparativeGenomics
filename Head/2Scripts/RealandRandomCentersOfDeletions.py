import time
import timeit
import os


import random

import statistics


# MajorArc: 5781 - 16569


first = open('breaks.csv','r')
head = first.readline().strip()


# reps = open("breaks_comp.csv", "w")
# reps.write("Group;Middle;")
# reps.write("%s\n" % head)


arr = []

while True:
	read1 = first.readline().strip()
	if read1 == "":
		break

	#reps.write("Real;")

	row1 = read1.split(';')

	# starts = ';'.join(row1[:13])
	# endura = ','.join(row1[13:])

	middle = int(row1[1]) + int(row1[3].replace('"',''))//2

	# reps.write("%s;%s;%s\n" % (middle, starts, endura ))

	arr.append(middle)


sd_real = statistics.stdev(arr)










#reps = open("MitoBreakDB_121219_rand.csv", "w")
#reps.write("Middle;")
# reps.write("Deletions;5_breakpoint;3_breakpoint;Deletion_length\n")




big = 0
same = 0


#megaaar = []

for i in range(10000):

	first = open('breaks.csv','r')
	head = first.readline().strip()

	arr = []

	while True:
		read1 = first.readline().strip()
		if read1 == "":
			break


		#reps.write("Rand;")

		row1 = read1.split(';')


		s = 0
		f = 40000

		while s < 5781 or f > 16569:
			middle = random.randint(5781,16569)

			fulllen = int(row1[3].replace('"',''))
			half = fulllen//2

			if fulllen % 2 == 1:
				s = middle - half
				f = middle + half
			else:
				s = middle - half
				f = middle + half + 1


				#n = ranstart + ranplus//2
				#middle = ranstart + ranplus//2


		


		arr.append(middle)



		#reps.write("%s;%s;%s;%s;%s\n" % (middle, row1[0], s, f, row1[3] ))
		#reps.write("%s;%s;%s;%s;%s\n" % (middle, row1[0], row1[1], row1[2], row1[3] ))


	sd_curr = statistics.stdev(arr)
	print(sd_curr)

	if sd_curr > sd_real:
		big+=1

	if sd_curr == sd_real:
		same+=1



print("")
print(sd_real)
print("")
print(same)
print("")
print(big)