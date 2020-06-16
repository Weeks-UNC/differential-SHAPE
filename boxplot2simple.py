#!/usr/bin/python
# 08/29/11 Gregg Rice
#
# This script converts the simple normalization method
# to the simple norm normalization. Method described below
#
# 1. rank order the data points
# 2. ignore the top 2% of data points
# 3. average the next 8% to find the normalization factor
#
# end of method

import sys,math,csv

if len(sys.argv) == 1:
	print 'Usage: <boxPlot.inp> <simpleNorm.out>'
	sys.exit()

input  = sys.argv[1]
output = sys.argv[2]

w = open(output,"w")

data = []

for line in csv.reader(open(input).readlines()[0:],delimiter = "\t"):
	data.append([line[0],float(line[1])])

def calcQuartile(x,q,qtype=7):
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
		if i[1]>-500:
			x.append(i[1])
	x.sort()
	tenPct = len(x)/10
	#calculate the interquartile range *1.5
	qLimit = 1.5*abs(calcQuartile(x,0.25)-calcQuartile(x,0.75))
	tenLimit = x[len(x)-1 - tenPct]
	#choose the cutoff that eliminates the fewest points
	limit = max(qLimit,tenLimit)
	#make new list without the outliers
	for i in range(len(x)):
		if x[i]<limit:
			o.append(x[i])
	#avg next ten percent
	for i in range(-tenPct,0):
		a = o[i] + a
	normFactor = a/tenPct
	print "Normalization Factor:",normFactor
	return normFactor

def findSimpleNormFactor(array):
	x,o,a = [],[],0
	for i in array:
		if i[1]>-500:
			x.append(i[1])
		else: x.append(0.0)
	x.sort()
	#print len(x)
	twoPct = int(round(len(x)*.02))
	eightPct = int(round(len(x)*.08))
	for i in range(0,len(x)-twoPct):
		o.append(x[i])
	for i in range(-eightPct,0):
		a = o[i] + a
	normFactor = a/eightPct
	print "Normalization Factor:", normFactor
	return normFactor

def normalizeData(array):
	normFactor = findSimpleNormFactor(array)
	for i in range(len(array)):
		if array[i][1] > -500:
			array[i][1] = array[i][1]/normFactor
	return array

norm = normalizeData(data)

for i in norm:
	line = "%s\t%s\n" % (i[0],round(i[1],3))
	w.write(line)
