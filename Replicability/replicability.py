#! /usr/bin/env python

from unittest import TestCase, makeSuite, main

from Model.Plot.performance import Performance
from numpy import where
from scipy.io import loadmat, savemat


class Replicability:
  def __init__(self, projectName="Mousefunc"):
    self.__projectname = projectName
  def __del__(self):
    del self.__projectname

  def loadPerformanceMatlab(self, fileName):
    performance = loadmat(fileName, squeeze_me=False, mat_dtype=True)
    performance = performance['performance']
    return performance

  def loadLabelsMatlab(self, fileName):
    labels = loadmat(fileName, squeeze_me=False, mat_dtype=True)
    labels = labels['labels']
    return labels

  def loadCVMatlab(self, fileName):
    cv = loadmat(fileName, squeeze_me=False, mat_dtype=True)
    cv = cv['CV']
    return cv

  def selectReproducibilitySet(self, label, score, proteinid, cv):
    cvall = where(cv == -1)[0]
    labelall = label[cvall]
    scoreall = score[cvall]
    proteinall = proteinid[cvall]

    indeces = where(labelall == 0.0)[0]
    protein = proteinall[indeces]
    score = scoreall[indeces]

    return (protein, score)

  def avgAUROCProtein(self, sets, proteinid, top=10):
    numset = len(sets)

    i = 0
    setid = arange(numset)
    meanroc = []
    while i < numset:
      # extract label-and score sets
      labelid = []
      for id in setid:
        if i != id:
          labelid.append(id)
				
      # create average scores
      avgscores = []
      k = 0
      for s in sets[0]:
        score = 0.0
        for id in labelid:
          score += sets[id][k]
        score = score/float(len(labelid))
        avgscores.append(score)
        k += 1
      x = sort(avgscores)

      # create label
      k = 0
      label = zeros(len(proteinid))
      for score in x[:top]:
        label[where(avgscores==score)] = 1
        index = where(avgscores==score)[0]

      roc = self.AUROCGillis(label, sets[i])
      meanroc.append(roc)
      i += 1
    roc_mean = reduce(lambda x, y: x + y / float(len(meanroc)), meanroc, 0)
    roc_std = std(meanroc)
    return (roc_mean, roc_std)

  def defineSample(self, numPred=10, simdir='./simulation_3_27_2015/'):
    while 1:
      ppi = randint(0, 5) 
      co = randint(0, 5)
      sm = 0
      if ppi+co < numPred:
        sm = numPred - (ppi + co)	
      else:
        co = numPred - ppi

      numAlg = randint(0, 6)
      if numAlg == 0:
        numAlg = 1
        filename1 = simdir+"simrep_" + str(ppi) + str(co) + str(sm) + str(numAlg) + ".mat"
        if not exists(filename1):
          return (ppi, co, sm, numAlg)

  def sample(self, ppi, co, sm, qAlg, goid, proteinid, top):
    projectsppi = ["HumanPPIBioGRID_8_20_2015", "HumanPPIHIPPIE_8_20_2015", \
                   "HumanPPIIntAct_8_20_2015", "HumanPPII2D_8_20_2015", "HumanPPIGM_8_20_2015"]

    projectssm = ["HumanSMKEGG_8_20_2015", "HumanSMReact_8_20_2015", \
                  "HumanSMCarta_8_20_2015", "HumanSMPfam_8_20_2015", "HumanSMInterPro_8_20_2015"]

    projectsco = ["HumanCoF20A_3_24_2015", "HumanCoF20B_3_24_2015", \
                  "HumanCoFGPL57020A_5_22_2014", "HumanCoFGPL57020B_5_22_2014"]
 
    algorithms = ["NeighborVoting", "GeneMANIA", "LogisticRegression", \
                  "SGD", "PassiveAggressive", "RWR"]

    predictors = []

    i = 0
    while i < ppi:
      perffile = "performance_" + algorithms[randint(0, qAlg)] + "_" + \
                  projectsppi[randint(0, 5)] + ".mat"
      perf = self.loadPerformanceMatlab(perffile)
      predictors.append(perf)

      i += 1

    i = 0
    while i < co:
      perffile = "performance_" + algorithms[randint(0, qAlg)] + "_" + \
					projectsco[randint(0, 4)] + ".mat"
      perf = self.loadPerformanceMatlab(perffile)
      predictors.append(perf)

      i += 1
	
    i = 0
    while i < sm:
      perffile = "performance_" + algorithms[randint(0, qAlg)] + "_" + \
                 projectssm[randint(0, 5)] + ".mat"
      perf = self.loadPerformanceMatlab(perffile)
      predictors.append(perf)

      i += 1
		
    # aggregation
    perf = predictors[0] + predictors[1] + predictors[2] + predictors[3] + predictors[4] + \
           predictors[5] + predictors[6] + predictors[7] + predictors[8] + predictors[9]
    perf = perf / float(10.0)

    labelsfilename = "performance_Label_Human_4_24_2014.mat"
    labels = self.loadLabelsMatlab(labelsfilename)
    print labels.shape

    cvfilename = "performance_CV_Human_4_24_2014.mat"
    cv = self.loadCVMatlab(cvfilename)
    print cv.shape

    aggregation = []
    i = 0
    for go in goid:
      print "_________ GOID= %d; i= %d _________" % (int(go[0]), i)

      label = labels[:,i]
      scores = perf[:, i]

      cv0 = where(cv == 0)[0]
      cv1 = where(cv == 1)[0]
      cv2 = where(cv == 2)[0]

      labelcv0 = label[cv0]
      labelcv1 = label[cv1]
      labelcv2 = label[cv2]	

      scorescv0 = scores[cv0]
      scorescv1 = scores[cv1]
      scorescv2 = scores[cv2]

      meanroc = []
      meanfmax = []

      roc = self.AUROCGillis(labelcv0, scorescv0)
      meanroc.append(roc)

      roc = self.AUROCGillis(labelcv1, scorescv1)
      meanroc.append(roc)

      roc = self.AUROCGillis(labelcv2, scorescv2)
      meanroc.append(roc)

      roc_mean = reduce(lambda x, y: x + y / float(len(meanroc)), meanroc, 0)

      aggregation.append(roc_mean)

      i += 1
		
    # Reproducibility	
    reproducibility = []
    i = 0
    for go in goid:
      sets = []

      (protein1, score1) = self.selectReproducibilitySet(labels[:,i], predictors[0][:,i], \
									proteinid, cv)
      (protein2, score2) = self.selectReproducibilitySet(labels[:,i], predictors[1][:,i], \
										proteinid, cv)
      (protein3, score3) = self.selectReproducibilitySet(labels[:,i], predictors[2][:,i], \
										proteinid, cv)
      (protein4, score4) = self.selectReproducibilitySet(labels[:,i], predictors[3][:,i], \
										proteinid, cv)
      (protein5, score5) = self.selectReproducibilitySet(labels[:,i], predictors[4][:,i], \
										proteinid, cv)
      (protein6, score6) = self.selectReproducibilitySet(labels[:,i], predictors[5][:,i], \
										proteinid, cv)
      (protein7, score7) = self.selectReproducibilitySet(labels[:,i], predictors[6][:,i], \
										proteinid, cv)
      (protein8, score8) = self.selectReproducibilitySet(labels[:,i], predictors[7][:,i], \
										proteinid, cv)
      (protein9, score9) = self.selectReproducibilitySet(labels[:,i], predictors[8][:,i], \
										proteinid, cv)
      (protein10, score10) = self.selectReproducibilitySet(labels[:,i], predictors[9][:,i], \
										proteinid, cv)
		
      sets.append(score1)
      sets.append(score2)
      sets.append(score3)
      sets.append(score4)
      sets.append(score5)
      sets.append(score6)
      sets.append(score7)
      sets.append(score8)
      sets.append(score9)
      sets.append(score10)

      (rocmean, rocstd) = self.avgAUROCProtein(sets, protein1, top=top)
      reproducibility.append(rocmean)
	
      i += 1

    return (reproducibility, aggregation)

  def simulation(self, N, numPred=10, top=10, simdir='./simulation_3_27_2015/'):
    go = Ontology()
    goid = go.loadGOIDMatlab(goDIR+"GOIDSlim_Human_4_24_2014.mat")
    proteinidfilename = "performance_ProteinID_Human_4_24_2014.mat"
    proteinid = go.loadProteinIDMatlab(proteinidfilename)
    del go

    i = 0
    while i < N:
      # select ratio of data resources
      (ppi, co, sm, numAlg) = self.defineSample()

      filename1 = simdir+"simrep_" + str(ppi) + str(co) + str(sm) + str(numAlg) + ".mat"
      filename2 = simdir+"simagg_" + str(ppi) + str(co) + str(sm) + str(numAlg) + ".mat"

      reproducibilitysampling = lil_matrix((goid.shape[0], numPred))
      aggregationsampling = lil_matrix((goid.shape[0], numPred))

      j = 0
      while j < numPred:
        (reproducibility, aggregation) = self.sampleRun(ppi, co, sm, numAlg, \
							goid, proteinid, top)

        reproducibility = asarray(reproducibility)
        reproducibility = reproducibility.reshape(len(reproducibility), 1)

        aggregation = asarray(aggregation)
        aggregation = aggregation.reshape(len(aggregation), 1)

        reproducibilitysampling[:, j] = reproducibility
        aggregationsampling[:, j] = aggregation

        j += 1

      savemat(filename1, {'sampling':reproducibilitysampling})
      savemat(filename2, {'sampling':aggregationsampling})
      i += 1

