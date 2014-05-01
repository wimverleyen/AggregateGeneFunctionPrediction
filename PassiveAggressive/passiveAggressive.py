#! /usr/bin/env python

from unittest import TestCase, makeSuite, main

from math import floor, sqrt, fabs
from numpy import count_nonzero, where, array_equal, transpose, asarray, float64
from numpy.random import permutation, seed

from sklearn.linear_model import PassiveAggressiveClassifier

class PassiveAggressive:
	def __init__(self, dbName="funcpred", goidTable="PassiveAggressive_Mousefunc_GOID", \
			proteinTable="PassiveAggressive_Mousefunc_Protein"):
		seed(66)
		self.__database = Database(dbName=dbName)
		self.__goidtable = goidTable
		self.__proteintable = proteinTable

		self.__proteins = []
                self.__numproteins = 0
                self.__goid = []
                self.__network = None
                self.__go2gene = None
	def __del__(self):
		del self.__database
		del self.__goidtable
		del self.__proteintable
		del self.__proteins
		del self.__numproteins
		del self.__goid
		del self.__network
		del self.__go2gene
	def loadGOData(self, proteinID, GOID, function, network):
                """
                        A new interface function to load GO data
                """
                go = Ontology()
                self.__proteins = go.loadProteinIDMatlab(proteinID)
                self.__numproteins = len(self.__proteins)
                self.__goid = go.loadGOIDMatlab(GOID)
                self.__go2gene = go.loadGO2GeneMatlab(function)
                del go
                net = Network()
                self.__network = net.loadNetworkMatlab(network)
                del net
        def selectAnnotatedProteins(self, goid):
                annotations = self.__go2gene[goid, :]
                annotations = transpose(annotations)
                print annotations.shape
                return annotations
	def selectAnnotatedProteinsMousefunc(self, goid):
                annotations = self.__go2gene[:, goid]
                return annotations
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
	def run(self, nFold=3, iter=10, verbose=1):
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
                        print goidindex
			annotations = self.selectAnnotatedProteinsMousefunc(goidindex)

			print "0s=", len([x for x in annotations if x == 0])
			print "1s=", len([x for x in annotations if x == 1])
			print "-1s=", len([x for x in annotations if x == -1])
						
			annotation = []
                        for value in annotations:
                                annotation.append(value)

			annotation = asarray(annotation).astype(float64)
                        annotation = annotation.ravel()

			model = PassiveAggressiveClassifier(loss='hinge', n_iter=iter, verbose=verbose)
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
				
				labeltmp = asarray(labeltmp).astype(float64)
                                labeltmp = labeltmp.ravel()
                                print labeltmp.shape			
	
				for index in ix:
					labeltmp[index] = 0

				print "0s=", len([x for x in labeltmp if x == 0])
				print "1s=", len([x for x in labeltmp if x == 1])
				print "-1s=", len([x for x in labeltmp if x == -1])

				model = PassiveAggressiveClassifier(loss='hinge', \
							n_iter=iter, verbose=verbose)
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
			#print sum(meanroc)/float(len(meanroc))
			pr_mean = reduce(lambda x, y: x + y / float(len(meanpr)), meanpr, 0)
			print "Mean AUPR= ", pr_mean
			#print sum(meanpr)/float(len(meanpr))
			fmax_mean = reduce(lambda x, y: x + y / float(len(meanfmax)), meanfmax, 0)
			print "Mean Fmax= ", fmax_mean

			self.__database.insertGOIDView(self.__goidtable, resultid, goid[0], nFold, \
						[roc_mean, pr_mean, fmax_mean])
			resultid += 1

			test += 1
			#if test == 5:
			#	break

class TestPassiveAggressive(TestCase):
	def setUp(self):
		projectname = "YEastPPI_9_11_2013"
                goidtable = "PassiveAggressive_" + projectname + "_GOID"
                proteintable = "PassiveAggressive_" + projectname + "_Protein"
                self.pahuman = PassiveAggressive(dbName="funcpred", goidTable=goidtable, \
                                        proteinTable=proteintable)
	def tearDown(self):
		del self.payeast
	def testARun(self):
		projectname = "Yeast_9_11_2013"
                proteinid = goDIR + "sgdID_" + projectname + ".mat"
                goidname = goDIR + "GOID_" + projectname + ".mat"
                go2genename = goDIR + "Gene2GO_" + projectname + ".mat"
                networkname = ppiDIR + "ppiNetwork_" + projectname + ".mat"
                self.pahuman.loadGOData(proteinid, goidname, go2genename, networkname)
                self.pahuman.run(iter=10, verbose=0)

def suite():
        suite = makeSuite(TestPassiveAggressive, 'test')
        return suite

if __name__ == '__main__':
        main(defaultTest='suite', argv=['-d'])
