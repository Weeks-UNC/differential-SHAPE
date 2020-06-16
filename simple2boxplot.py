#!/usr/bin/python
# 06/27/12 Gregg Rice, gmr@unc.edu,
#          UNC Chemistry - Weeks Lab
#
# This script converts the simple normalization method
# to the boxplot normalization. Method described below
#
# 1.find the scaling factor by ranking  all of the shape values
# 2. take 1.5*abs(Q1-Q3) as a cutoff
# 3. remove either the top 10% of the RNA or the positions above this cutoff, whichever is smaller
# 3.1 if RNA is less than 100 bp, use 5%
# 4. Average the next 10% from the origional length of the RNA --> this is the scaling factor
#
# end of method

import sys,math,csv



def calcQuartile(x,q,qtype=1):
	#source: http://adorio-research.org/wordpress/?p=125
	# x = array, q = quartile (in % as a decimal)
	y=x
	n = len(y)
	abcd = [(0,   0, 1, 0), # inverse empirical distrib.function., R type 1
	  (0.5, 0, 1, 0), # similar to type 1, averaged, R type 2
	  (0.5, 0, 0, 0), # nearest order statistic,(SAS) R type 3
	  (0,   0, 0, 1), # California linear interpolation, R type 4
	  (0.5, 0, 0, 1), # hydrologists method, R type 5
	  (0,   1, 0, 1), # mean-based estimate(Weibull method), (SPSS,Minitab), type 6 
	  (1,  -1, 0, 1), # mode-based method,(S, S-Plus), R type 7
	  (1.0/3, 1.0/3, 0, 1), # median-unbiased ,  R type 8
	  (3/8.0, 0.25, 0, 1)   # normal-unbiased, R type 9.
	 ]
	a, b, c, d = abcd[qtype-1]
	g, j = math.modf( a + (n+b) * q -1)
	if j<0:
		return x[0]
	elif j>=n:
		return x[n-1]
	j = int(math.floor(j))
	if g == 0:
		return x[j]
	else:
		return y[j] + (y[j+1]- y[j])* (c + d * g)

def findBoxplotFactor(array):
	x,o,a = [],[],0
	for i in array:
		if i<-500: i = 0.00
		x.append(i)
	x.sort()
	tenPct = len(x)/10
	fivePct = len(x)/20
	#calculate the interquartile range *1.5
	qLimit = 1.5*abs(calcQuartile(x,0.25)-calcQuartile(x,0.75))
	tenLimit = x[len(x)-1 - tenPct]
	fiveLimit = x[len(x)-1 - fivePct]
	#choose the cutoff that eliminates the fewest points
	#print qLimit, tenLimit
	limit = max(qLimit,tenLimit)
	if len(x)<100:
		limit = max(qLimit,fiveLimit)
		print '5limit'
	#make new list without the outliers
	for i in range(len(x)):
		if x[i]<limit:
			o.append(x[i])
	#avg next ten percent
	for i in range(-tenPct,0):
		a = o[i] + a
	normFactor = a/tenPct
	#print "Normalization Factor:",normFactor
	return normFactor

def normalizeData(array):
	normFactor = findBoxplotFactor(shapeArray)
	for i in range(len(array)):
		if array[i][1] > -500:
			array[i][1] = array[i][1]/normFactor
	return array

def renormalizeArray(array,bpfactor=False):
	if not bpfactor:
		bpfactor = findBoxplotFactor(array)
	for i in range(len(array)):
		if array[i] > -500:
			array[i] = array[i]/bpfactor
	return array,bpfactor

if __name__ == '__main__':
	if len(sys.argv) < 3:
		print 'Usage:<input> <output>'
		sys.exit()
	
	input  = sys.argv[1]
	output = sys.argv[2]
	
	w = open(output,"w")
	
	data = []
	shapeArray = []
	seqArr = []
	errArr = []
	
	for i in open(input).readlines():
		line = i.rstrip().split()
		data.append([line[0],float(line[1])])
		shapeArray.append(float(line[1]))
		try:
			errArr.append(float(line[2]))
			seqArr.append(float(line[3]))
		except:
			pass
	
	#for line in csv.reader(open(input).readlines()[0:],delimiter = "\t"):
	#	data.append([line[0],float(line[1])])
	#	shapeArray.append(float(line[1]))

	norm = normalizeData(data)
	
	for i in norm:
		line = "%s\t%s\n" % (i[0],round(i[1],3))
		w.write(line)
