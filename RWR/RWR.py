#! /usr/bin/env python

from unittest import TestCase, makeSuite, main

from math import floor, sqrt, fabs
from numpy import zeros, count_nonzero, where, array_equal, transpose, asarray, float64, fill_diagonal, dot, copy
from numpy.random import permutation, seed, shuffle
from numpy.linalg import inv
from scipy.sparse import lil_matrix
from scipy.io import loadmat, savemat

biogridDIR = "/home/wverleye/work/scripts/Data/BioGRID/"

class RWR:
	def __init__(self, projectname, dbName="funcpred", goidTable="SGD_Mousefunc_GOID", \
			proteinTable="SGD_Mousefunc_Protein"):
		seed(66)
		self.__projectname = projectname
		self.__database = Database(dbName=dbName)
		self.__goidtable = goidTable
		self.__proteintable = proteinTable
		self.__epsilon = 0.0
		self.__theta = 0.0
	def __del__(self):
		del self.__database
		del self.__goidtable
		del self.__proteintable
		del self.__epsilon
		del self.__theta

	def loadGOData(self, proteinID, GOID, function, network):
		"""
			An interface function to load GO data
		"""
	        go = Ontology()
                self.__proteins = go.loadProteinIDMatlab(proteinID)
                self.__numproteins = len(self.__proteins)
                self.__goid = go.loadGOIDMatlab(GOID)
                self.__go2gene = go.loadGO2GeneMatlab(function)
		net = Network()
		self.__network = net.loadNetworkMatlab(network)
		del net
                del go

	def selectAnnotatedProteinsMousefunc(self, goid):
                annotations = self.__go2gene[:, goid]
                return annotations	

	def model(self, epsilon, theta):
		"""
			epsilon: convergence parameter
			theta: restart probability
		"""

		self.__epsilon = epsilon
		self.__theta = theta

	def createTransitionMatrix(self):
		(x, y) = self.__network.shape

		D = lil_matrix(self.__network.shape)

		nodedegree = self.__network.sum(axis=1)
		nodedegree = asarray(nodedegree)
		nodedegree = nodedegree.ravel()

		D.setdiag(nodedegree)
		Q = dot(inv(D.todense()), self.__network)

		filename = "transition_" + self.__projectname + ".mat"
		savemat(filename, {'Q':Q})

	def fit(self, network, label):
		# Load transition matrix
		filename = "transition_" + self.__projectname + ".mat"
		Q = loadmat(filename, squeeze_me=False, mat_dtype=True)
		Q = Q['Q']

		print Q.shape

		# Initialize p_t
		(x, y) = network.shape
		p_t1 = zeros(x)
		Vm = len([x for x in label if x == 1])
		indeces = where(label==1)[0]
		p_t1[indeces] = 1/float(Vm)
		p_t0 = copy(p_t1)

		i = 0
		d1 = 0.0
		while 1:
			p_t2 = (1 - self.__theta) * dot(Q.T, p_t1) + self.__theta * p_t0 
			d2 = sum(p_t2-p_t1)

			i += 1
			if i == 100:
				return p_t2
			elif d2 == 0.0:
				return p_t2
			elif (d2+d1) == 0.0:
				return p_t2
			d1 = sum(p_t2-p_t1)
			p_t1 = p_t2

	def runMousefunc(self, nFold=3):
		"""
			CV: -1 => total model (no cv)
			CV: nFold => mean metric over cv
		"""
		self.__database.createGOIDView(self.__goidtable, double=["AUROC", "AUPR", "Fmax"], drop=True)
		self.__database.createProteinView(self.__proteintable, \
						double=["ProteinID", "Label", "Score"], drop=True)
		
		# Get labels
		test = 0
		pp = permutation(self.__numproteins)
		resultid = 0
		for goid in self.__goid:
			print "____________ GOID= %d ____________" % goid
			# Get label for GOID
			goidindex = where(self.__goid==goid)
			goidindex = int(goidindex[0])
			annotations = self.selectAnnotatedProteinsMousefunc(goidindex)

			print "0s=", len([x for x in annotations if x == 0])
			print "1s=", len([x for x in annotations if x == 1])
			print "-1s=", len([x for x in annotations if x == -1])

			annotation = []
                        for value in annotations:
                                annotation.append(value)

			annotation = asarray(annotation).astype(float64)
                        annotation = annotation.ravel()

			self.model(0.01, 0.8)
			score = self.fit(self.__network, annotation)

			scores = []
			for value in score:
				scores.append(1-value)

			per = Performance(annotations, scores)
			roc = per.AUROCGillis()
                        print "AUROC= ", roc
                        pr = per.AUPRGillis()
                        print "AUPR= ", pr
			fmax = per.Fmax()
                        print "Fmax= ", fmax

			self.__database.insertProteinView(self.__proteintable, resultid, goid[0], -1, \
						self.__proteins, annotations, scores)
			self.__database.insertGOIDView(self.__goidtable, resultid, goid[0], -1, [roc, pr, fmax])
			resultid += 1

			del per

			labelIx = range(self.__numproteins)
			offset = 0
			fold = 0
			meanroc = []
			meanpr = []
			meanfmax = []

			while fold < nFold:
				print "____________ Fold= %d ____________" % fold
				lastelem = min(self.__numproteins, offset+floor(self.__numproteins/nFold))
				ix = []
				for index in pp[offset+1:lastelem]:
					ix.append(labelIx[index])
				
				offset = lastelem
	
				labeltmp = []
				for value in annotations:
					labeltmp.append(float(value))
	
				for index in ix:
					labeltmp[index] = 0

				labeltmp = asarray(labeltmp).astype(float64)
				labeltmp = labeltmp.ravel()

				print "0s=", len([x for x in labeltmp if x == 0])
				print "1s=", len([x for x in labeltmp if x == 1])
				print "-1s=", len([x for x in labeltmp if x == -1])

				self.model(0.01, 0.8)
				scores = self.fit(self.__network, labeltmp)

				score = []
				annotation = []
				proteins = []
				for index in ix:
					score.append(float(1-scores[index]))
					annotation.append(annotations[index])
					proteins.append(self.__proteins[index])

				per = Performance(annotation, score)
				roc = per.AUROCGillis()
                	        print "AUROC= ", roc
				meanroc.append(roc)
                        	pr = per.AUPRGillis()
	                        print "AUPR= ", pr			
				meanpr.append(pr)
				fmax = per.Fmax()
	                        print "Fmax= ", fmax
				meanfmax.append(fmax)

				self.__database.insertGOIDView(self.__goidtable, resultid, goid[0], fold,\
									[roc, pr, fmax])
				self.__database.insertProteinView(self.__proteintable, resultid, goid[0],\
								fold, proteins, annotation, score)

				del proteins
				del annotation
				del score
				del per
				fold += 1
				resultid += 1

			roc_mean = reduce(lambda x, y: x + y / float(len(meanroc)), meanroc, 0)
			print "Mean AUROC= ", roc_mean
			pr_mean = reduce(lambda x, y: x + y / float(len(meanpr)), meanpr, 0)
			print "Mean AUPR= ", pr_mean
			fmax_mean = reduce(lambda x, y: x + y / float(len(meanfmax)), meanfmax, 0)
			print "Mean Fmax= ", fmax_mean
			
			self.__database.insertGOIDView(self.__goidtable, resultid, goid[0], nFold, \
								[roc_mean, pr_mean, fmax_mean])
			resultid += 1

			test += 1


class TestRWR(TestCase):
	def setUp(self):
		projectname = "HumanPPIBioGRID_8_20_2015"
		goidtable = "RWR_" + projectname + "_GOID"
                proteintable = "RWR_" + projectname + "_Protein"
		self.rwrhuman = RWR(projectname, dbName="funcpred", goidTable=goidtable, \
					proteinTable=proteintable)

	def tearDown(self):
		del self.rwrhuman
	def testARun(self):
		projectname = "Human_8_20_2015"
		proteinid = goDIR + "entrezGeneID_" + projectname + ".mat"
                goidname = goDIR + "GOID_" + projectname + ".mat"
                go2genename = goDIR + "Gene2GO_" + projectname + ".mat"
                networkname = biogridDIR + "ppiNetwork_" + projectname + ".mat"
		self.rwrhuman.loadGOData(proteinid, goidname, go2genename, networkname)
		self.rwrhuman.createTransitionMatrix()
		self.rwrhuman.runMousefunc()


def suite():
        suite = makeSuite(TestRWR, 'test')
        return suite

if __name__ == '__main__':
        main(defaultTest='suite', argv=['-d'])
