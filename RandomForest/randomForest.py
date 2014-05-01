#! /usr/bin/env python

from unittest import TestCase, makeSuite, main

from math import floor
from numpy.random import permutation, seed
from numpy import count_nonzero, delete, where, transpose, asarray, float64
from sklearn.ensemble import RandomForestClassifier

class RandomForest:
	def __init__(self, dbName="funcpred", goidTable="RandomForest_Mousefunc_GOID",\
				proteinTable="RandomForest_Mousefunc_Protein"):
		seed(66)
		self.__database = Database(dbName=dbName)
		self.__goidtable = goidTable
		self.__proteintable = proteinTable

		self.__leafnodes = {}
		
		self.__proteins = []
                self.__numproteins = 0
                self.__goid = []
                self.__go2gene = None
	def __del__(self):
		del self.__database
		del self.__goidtable
		del self.__proteintable
		del self.__leafnodes
		del self.__proteins
		del self.__numproteins
		del self.__goid
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
		net = Network()
		self.__network = net.loadNetworkMatlab(network)
		del net
                del go
	def selectAnnotatedProteins(self, goid):
		annotations = self.__go2gene[goid, :]
		annotations = transpose(annotations)
		print annotations.shape
		return annotations
	def selectAnnotatedProteinsMousefunc(self, goid):
		annotations = self.__go2gene[:, goid]
		print annotations.shape
		return annotations
	def runMousefunc(self, nFold=3, nTree=250, criterion="gini", density=0.1):
		"""
			CV: -1 => total model (no cv)
			CV: nFold => mean metric over cv
		"""
		self.__database.createGOIDView(self.__goidtable, double=["AUROC", "AUPR", "Fmax"], drop=True)
		self.__database.createProteinView(self.__proteintable, double=["ProteinID", "Label", "Score"], \
							drop=True)

		# Get labels
		test = 0
		pp = permutation(self.__numproteins)
		resultid = 0
		for goid in self.__goid:
			print "____________ GOID %d ____________" % goid
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
                        print annotation.shape
	
			#model = RandomForestClassifier(n_estimators=100, n_jobs=-1)
			model = RandomForestClassifier(n_estimators=nTree, criterion=criterion, compute_importances=True,\
							min_density=density, bootstrap=True, max_features="auto", \
							n_jobs=8, verbose=1)
			model.fit(self.__network, annotation)
			
			scores = model.predict_proba(self.__network)
			per = Performance(annotations, scores[:,0])
			roc = per.AUROCGillis()
                        print "AUROC= ", roc
                        pr = per.AUPRGillis()
                        print "AUPR= ", pr
			fmax = per.Fmax()
                        print "Fmax= ", fmax

			self.__database.insertProteinView(self.__proteintable, resultid, goid[0], -1, self.__proteins, \
								annotations, scores[:,0])
			self.__database.insertGOIDView(self.__goidtable, resultid, goid[0], -1, [roc, pr, fmax])
			resultid += 1

			del per
			#break

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
				
				model = RandomForestClassifier(n_estimators=nTree, criterion=criterion, compute_importances=True,\
                                                                min_density=density, bootstrap=True, max_features="auto", \
                                                                n_jobs=8, verbose=1)
				model.fit(self.__network, labeltmp)
				scores = model.predict_proba(self.__network)
				
				scores = scores[:,0]

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

				self.__database.insertGOIDView(self.__goidtable, resultid, goid[0], fold, [roc, pr, fmax])
				self.__database.insertProteinView(self.__proteintable, resultid, goid[0], fold, \
									proteins, annotation, score)

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
			
class TestRandomForest(TestCase):
	def setUp(self):
		projectname = "YeastPPI_9_11_2013"
                goidtable = "RandomForest_" + projectname + "_GOID"
                proteintable = "RandomForest_" + projectname + "_Protein"
		self.rfyeast = RandomForest(dbName="funcpred", goidTable=goidtable, \
					proteinTable=proteintable)
	def tearDown(self):
		del self.rfyeast
	def testARun(self):
		projectname = "Yeast_9_11_2013"
		proteinid = goDIR + "sgdID_" + projectname + ".mat"
		goid = goDIR + "GOID_" + projectname + ".mat"
                function = goDIR + "Gene2GO_" + projectname + ".mat"
                networkname = coDIR + "ppiNetwork_" + projectname + ".mat"
		self.rfyeast.loadGOData(proteinid, goid, function, networkname)
		self.rfyeast.runMousefunc(nTree=300, criterion="gini")

def suite():
        suite = makeSuite(TestRandomForest, 'test')
        return suite

if __name__ == '__main__':
        main(defaultTest='suite', argv=['-d'])
