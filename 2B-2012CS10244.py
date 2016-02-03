import math
# trainingFile = list(open("training.txt"))

def get_islandLocs(cpg_testFile):
	cpg_testList = list(open(cpg_testFile,'r'))
	# print cpg_testList[2].replace('\n','')
	islandLocs = []#[[int(s),int(e)] for ]
	for island in cpg_testList:
		island = island.replace('\n','')
		line = island.split()
		# print line
		s = int(line[0])
		e = int(line[1])
		islandLocs.append([s,e])
	return islandLocs
# print get_islandLocs("cpg_island_training.txt")

def IsInCpG_island(charNumber,islandLoc):
	res = False
	for locs in islandLoc:
		if(charNumber>=locs[0] and charNumber<=locs[1]):
			res = True
			break
	return res

def trainHMM(cpg_testFile,trainingFile):
	trainingList = list(open(trainingFile,'r'))
	islandLocs = get_islandLocs(cpg_testFile)
	
	transMatrix = {"A+":{},"A-":{},"C+":{},"C-":{},"G+":{},"G-":{},"T+":{},"T-":{}} # 8*8
	for key in transMatrix:
		transMatrix[key] = {"A+":0,"A-":0,"C+":0,"C-":0,"G+":0,"G-":0,"T+":0,"T-":0}

	emiMatrix = {"A+":{},"A-":{},"C+":{},"C-":{},"G+":{},"G-":{},"T+":{},"T-":{}} # 8*4
	for key in emiMatrix:
		emiMatrix[key] = {'A':0,'C':0,'G':0,'T':0}

	startStateMatrix ={"A+":0,"A-":0,"C+":0,"C-":0,"G+":0,"G-":0,"T+":0,"T-":0} #8*1

	baseNumber = 0;
	lastBase = ''
	lastBaseInCpG = ''
	for line in trainingList:
		line = line.replace('\n','')
		for i in range(len(line)):
			baseNumber += 1
			inCpG = IsInCpG_island(baseNumber,islandLocs) 
			if(baseNumber==1):
				lastBase = line[i]
				if(inCpG):
					lastBaseInCpG = '+'
				else:
					lastBaseInCpG = '-'
					startStateMatrix[lastBase+'-'] += 1
			else:
				curBase = line[i]
				if(inCpG):
					transMatrix[lastBase+lastBaseInCpG][curBase+'+'] += 1
					lastBaseInCpG ='+'
				else:
					transMatrix[lastBase+lastBaseInCpG][curBase+'-'] += 1
					lastBaseInCpG ='-'
					startStateMatrix[curBase+'-'] += 1
				lastBase = curBase
	
	for key in emiMatrix:
		for value in emiMatrix[key]:
			if(key[0]==value):
				emiMatrix[key][value] = 1000

	sumMatrix = {"A+":0,"A-":0,"C+":0,"C-":0,"G+":0,"G-":0,"T+":0,"T-":0}
	for key in transMatrix:
		for subkey in transMatrix[key]:
			sumMatrix[key] += transMatrix[key][subkey]
			if(transMatrix[key][subkey]==0):
				transMatrix[key][subkey]=1
	for key in transMatrix:
		# print "sumMatrix: ",key,sumMatrix[key]
		for subkey in transMatrix[key]:
			transMatrix[key][subkey] = (transMatrix[key][subkey]/float(sumMatrix[key]))*1000

	# for key in transMatrix:
	# 	for value in transMatrix[key]:
	# 		transMatrix[key][value] = (float(transMatrix[key][value]+1)/float(baseNumber-1))*1000
	
	totalNonCpG = sum([startStateMatrix[key] for key in startStateMatrix])
	# print "totalNonCpG: ",totalNonCpG
	for key in startStateMatrix:
		startStateMatrix[key] = (float(startStateMatrix[key])/float(totalNonCpG))*1000
		# print key,startStateMatrix[key]

	return [transMatrix,emiMatrix,startStateMatrix]

def printMatrix(matrix):
	for key in matrix:
		for value in matrix[key]:
			print matrix[key][value]
		print "\n"


def viterbiAlgo(trans_prob,emi_prob,start_prob,observationsFile):
	observationsList = list(open(observationsFile,'r'))
	states=["A+","A-","C+","C-","G+","G-","T+","T-"]
	observations = (''.join(observationsList)).replace('\n','')
	T = {}
	for s in states:
		T[s] = {"prob":start_prob[s],"v_path":[s],"v_prob":start_prob[s]}
		# print "T: ",s,T[s]["v_prob"]

	for i in range(len(observations)):
		U = {}
		for ns in states:
			total = 0
			argmax = []
			valmax = 0
			for s in states:
				prob = T[s]["prob"]
				v_path = T[s]["v_path"]
				v_prob = T[s]["v_prob"]
				p = emi_prob[s][observations[i]]*trans_prob[s][ns]
				# print "p: ",emi_prob[s][observations[i]],trans_prob[s][ns]
				prob  *= p
				v_prob *= p
				total += prob
				if(v_prob>valmax):
					argmax = v_path + [ns]
					valmax = v_prob

			U[ns] = {"prob":total,"v_path":argmax,"v_prob":valmax}
		normalizingConstant = 0
		v_normalizingConstant = 0
		for key in U:
			v_normalizingConstant += U[key]["v_prob"]
			# print key,U[key]["v_prob"]

		for key in U:
			U[key]["v_prob"] = U[key]["v_prob"]/float(v_normalizingConstant)
			# U[key]["v_prob"] = U[key]["v_prob"]/float(1)
		T = U
	total = 0
	argmax = []
	valmax = 0
	for s in states:
		prob = T[s]["prob"]
		v_path = T[s]["v_path"]
		v_prob = T[s]["v_prob"]
		total += prob
		if(v_prob > valmax):
			argmax = v_path
			valmax = v_prob
	predFile = open("results.txt",'w')
	for r in argmax:
		predFile.write(r+" ")
	return [argmax,math.log(valmax)]

def getCpGIslands(hiddenStates):
	tillPos = False
	start = 0
	end = 0
	CpG_islands = []
	for i in range(1,len(hiddenStates)):
		if(hiddenStates[i][1]!=hiddenStates[i-1][1]):
			if(not tillPos):
				start = i + 1
				tillPos = True
			else:
				end = i + 1
				tillPos = False
				CpG_islands.append([start,end])
	if(hiddenStates[len(hiddenStates)-1][1]=='+'):
		end = len(hiddenStates)
	CpG_islands.append([start,end])
	return CpG_islands

# trainedData = trainHMM("cpg_island_training.txt","training.txt")
# res = viterbiAlgo(trainedData[0],trainedData[1],trainedData[2],"testing.txt")
# CpG_islands = getCpGIslands(res[0])

trainData = raw_input("Training File: ")
CpGTrain = raw_input("CpG Island File: ")
testData = raw_input("Testing File: ")
trainedData = trainHMM(CpGTrain,trainData)
res = viterbiAlgo(trainedData[0],trainedData[1],trainedData[2],testData)
CpG_islands = getCpGIslands(res[0])
print "Predicted CpG Islands are: "
print CpG_islands
print "Log-Probability: " ,res[1]


