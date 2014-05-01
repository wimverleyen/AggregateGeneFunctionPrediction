#! /usr/bin/env python

from unittest import TestCase, makeSuite, main

from math import floor, sqrt, fabs
from numpy import count_nonzero, where, array_equal, transpose, asarray, float64
from numpy.random import permutation, seed, shuffle

from sklearn.linear_model import SGDClassifier

class StochasticGradientDescent:
	def __init__(self, dbName="funcpred", goidTable="SGD_Mousefunc_GOID", \
			proteinTable="SGD_Mousefunc_Protein"):
		seed(66)
		self.__database = Database(dbName=dbName)
		self.__goidtable = goidTable
		self.__proteintable = proteinTable

		self.__proteins = []
                self.__numproteins = 0

		self.__goid = []
	
		self.__network = None
	def __del__(self):
		del self.__database
		del self.__goidtable
		del self.__proteins
		del self.__numproteins
		del self.__goid
		del self.__network
	def loadMousefunc(self, network):
                mousefunc = MouseFunc()
                self.__proteins = mousefunc.proteins()
                self.__numproteins = mousefunc.numProteins()
                self.__goid = mousefunc.loadGOIDMat()
                self.__network = mousefunc.loadNetwork(network)
                self.__go2gene = mousefunc.getGO2Gene()
                del mousefunc
	def loadGOData(self, proteinID, GOID, function, network):
		"""
			A new interface function to load GO data
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
	def convertScore(self, score):
                tmp = []
                for s in score:
                        tmp.append(fabs(s))
                mag = max(tmp)
                del tmp

                scores = []
                for s in score:
                        value = s/float(mag)
                        scores.append(1 - ((value/float(2))+0.5))
                return scores
        def selectAnnotatedProteinsMousefunc(self, goid):
                annotations = self.__go2gene[:, goid]
                return annotations	
	def selectAnnotatedProteins(self, goid):
		annotations = self.__go2gene[goid, :]
		annotations = transpose(annotations)
		print annotations.shape
		return annotations
	def run(self, nFold=3, iter=10, verbose=1, loss='modified_huber', penalty='l2', shuffle=True):
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

			model = SGDClassifier(loss=loss, class_weight='auto', penalty=penalty, \
						n_iter=iter, shuffle=shuffle, verbose=verbose)
			model.fit(self.__network, annotation)
			scores = model.decision_function(self.__network)
                        scores = self.convertScore(scores)
			
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

				print "0s=", len([x for x in labeltmp if x == 0])
				print "1s=", len([x for x in labeltmp if x == 1])
				print "-1s=", len([x for x in labeltmp if x == -1])

				model = SGDClassifier(loss=loss, class_weight='auto', penalty=penalty, \
							n_iter=iter, shuffle=shuffle, verbose=verbose)
				model.fit(self.__network, labeltmp)
				scores = model.decision_function(self.__network)
	                        scores = self.convertScore(scores)

				score = []
				annotation = []
				proteins = []
				for index in ix:
					score.append(float(scores[index]))
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
			#if test == 5:
			#	break
	
class TestStochasticGradientDescent(TestCase):
	def setUp(self):
		projectname = "YeastPPI_9_11_2013"
		goidtable = "SGD_" + projectname + "_GOID"
                proteintable = "SGD_" + projectname + "_Protein"
		self.sgdyeast = StochasticGradientDescent(dbName="funcpred", goidTable=goidtable, \
					proteinTable=proteintable)
	def tearDown(self):
		del self.sgdyeast
	def testARun(self):
		projectname = "Yeast_9_11_2013"
		proteinid = goDIR + "sgdID_" + projectname + ".mat"
                goidname = goDIR + "GOID_" + projectname + ".mat"
                go2genename = goDIR + "Gene2GO_" + projectname + ".mat"
                networkname = ppiDIR + "ppiNetworkRandom_" + projectname + ".mat"
		self.sgdyeast.loadGOData(proteinid, goidname, go2genename, networkname)
		self.sgdyeast.run(iter=10, verbose=0, loss='hinge', penalty='l2')

def suite():
        suite = makeSuite(TestStochasticGradientDescent, 'test')
        return suite

if __name__ == '__main__':
        main(defaultTest='suite', argv=['-d'])
