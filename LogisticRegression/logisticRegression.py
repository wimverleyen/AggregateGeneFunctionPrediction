#! /usr/bin/env python

from unittest import TestCase, makeSuite, main

from math import floor
from numpy import count_nonzero, where, array_equal, transpose
from numpy.random import permutation, seed

from joblib import Parallel, delayed

from sklearn.linear_model import LogisticRegression

class LogRegression:
	def __init__(self, dbName="funcpred", goidTable="LogisticRegression_Mousefunc_GOID",\
				proteinTable="LogisticRegression_Mousefunc_Protein"):
		seed(66)
		self.__database = Database(dbName=dbName)
		self.__goidtable = goidTable
		self.__proteintable = proteinTable

		self.__proteins = []
                self.__numproteins = 0

		self.__goid = []

		self.__go2gene = None
		self.__network = None
	def __del__(self):
		del self.__database
		del self.__goidtable
		del self.__proteintable
		del self.__proteins
		del self.__numproteins
		del self.__goid
		del self.__go2gene
		del self.__network
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
	def selectAnnotatedProteinsMousefunc(self, goid):
		annotations = self.__go2gene[:, goid]
		return annotations
	def run(self, nFold=3, dual=True, penalty='l2', C=1.0, tol=0.0001):
		"""
			CV: -1 => total model (no cv)
			CV: nFold => mean metric over cv
		"""
		self.__database.createGOIDView(self.__goidtable, double=["AUROC", "AUPR", "Fmax"], drop=True)
		self.__database.createProteinView(self.__proteintable, double=["ProteinID", "Label", "Score"], \
							drop=True)
		
		print "network shape= ", self.__network.shape

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

			model = LogisticRegression(penalty=penalty, dual=dual, class_weight="auto", C=C, tol=tol)
			model.fit(self.__network, annotation)
			scores = model.predict_proba(self.__network)
			
			per = Performance(annotations, scores[:,0])
			roc = per.AUROCGillis()
                        print "AUROC= ", roc
                        pr = per.AUPRGillis()
                        print "AUPR= ", pr
			fmax = per.Fmax()
                        print "Fmax= ", fmax

			self.__database.insertProteinView(self.__proteintable, resultid, goid[0], -1,\
								self.__proteins, annotations, scores[:,0])
			self.__database.insertGOIDView(self.__goidtable, resultid, goid[0], -1, \
							[roc, pr, fmax])
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

				model = LogisticRegression(penalty=penalty, dual=dual, class_weight="auto", C=C, tol=tol)
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

				self.__database.insertGOIDView(self.__goidtable, resultid, goid[0], fold,\
								 [roc, pr, fmax])
				self.__database.insertProteinView(self.__proteintable, resultid, goid[0], fold, \
									proteins, annotation, score)

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
	def runMousefuncTestSet(self, dual=True, penalty='l2', C=1.0, tol=0.0001):
		"""
			CV: -1 => total model (no cv)
			CV: nFold => mean metric over cv
		"""
		self.__database.createGOIDView(self.__goidtable, double=["AUROC", "AUPR", "Fmax", "P20R"], drop=True)
		self.__database.createProteinView(self.__proteintable, double=["ProteinID", "Label", "Score"], \
							drop=True)
		
		print "network shape= ", self.__network.shape

		mf = MouseFunc()
		testset = mf.getTestSet()
		print testset.shape
	
		proteintestset = mf.proteinsTestSet()
                proteins = self.__proteins.transpose()
                proteins = proteins[0]
	
		del mf

		i = 0
                indecestestset = []
                for protein in proteins:
                        if protein in proteintestset:
                                indecestestset.append(i)
                        i += 1

		print len(indecestestset)

		test = 0
		resultid = 0
		for goid in self.__goid:
			print "____________ GOID= %d ____________" % goid
			# Get label for GOID
			goidindex = where(self.__goid==goid)
                        goidindex = int(goidindex[0])
                        print goidindex
			annotations = self.selectAnnotatedProteinsMousefunc(goidindex)

			print "0s= ", len([x for x in annotations if x == 0])
			print "1s= ", len([x for x in annotations if x == 1])
			print "-1s= ", len([x for x in annotations if x == -1])

			annotation = []
			for value in annotations:
				annotation.append(value)

			model = LogisticRegression(penalty=penalty, dual=dual, class_weight="auto", C=C, tol=tol)
			model.fit(self.__network, annotation)
			scores = model.predict_proba(self.__network)
			scores = scores[:,0]

			per = Performance(annotation, scores)
			roc = per.AUROCGillis()
                        print "AUROC= ", roc
                        pr = per.AUPRGillis()
                        print "AUPR= ", pr
			fmax = per.Fmax()
                        print "Fmax= ", fmax
			p20r = per.precision(recall=0.2)
                        print "P20R= ", p20r

			self.__database.insertProteinView(self.__proteintable, resultid, goid[0], -1,\
								self.__proteins, annotation, scores)
			self.__database.insertGOIDView(self.__goidtable, resultid, goid[0], -1, \
							[roc, pr, fmax, p20r])
			resultid += 1

			del per

			annotationheldout = testset[goidindex,:]
                        print annotationheldout.shape

			annotationtest = []
                        k = 0
                        for annotation in annotations:
                                annotationtest.append(annotation)
                                if annotationheldout[k] == 1:
                                        annotationtest[k] = 1
                                if annotationheldout[k] == -1:
                                        annotationtest[k] = 0

                                k += 1

                        print "0s=", len([x for x in annotationtest if x == 0])
                        print "1s=", len([x for x in annotationtest if x == 1])
                        print "-1s=", len([x for x in annotationtest if x == -1])

                        annotationtestset = []
                        scoretestset = []
                        proteintestset = []
                        k = 0
                        for protein in annotationtest:
                                if k in indecestestset:
                                        annotationtestset.append(annotationtest[k])
                                        scoretestset.append(scores[k])
                                        proteintestset.append(proteins[k])
                                k += 1

			print "0s=", len([x for x in annotationtestset if x == 0])
                        print "1s=", len([x for x in annotationtestset if x == 1])
                        print "-1s=", len([x for x in annotationtestset if x == -1])

			per = Performance(annotationtestset, scoretestset)
			roc = per.AUROCGillis()
                        print "AUROC= ", roc
                        pr = per.AUPRGillis()
                        print "AUPR= ", pr
			fmax = per.Fmax()
                        print "Fmax= ", fmax
			p20r = per.precision(recall=0.2)
                        print "P20R= ", p20r

			self.__database.insertProteinView(self.__proteintable, resultid, goid[0], 3,\
								proteintestset, annotationtestset, scoretestset)
			self.__database.insertGOIDView(self.__goidtable, resultid, goid[0], 3, \
							[roc, pr, fmax, p20r])
			resultid += 1

			del per

			test += 1
			#if test == 1:
			#	break
	
class TestLogRegression(TestCase):
	def setUp(self):
		projectname = "YeastPPI_9_11_2013"
                goidtable = "LogisticRegression_" + projectname + "_GOID"
                proteintable = "LogisticRegression_" + projectname + "_Protein"
		self.lryeast = LogRegression(dbName="funcpred", goidTable=goidtable, \
					proteinTable=proteintable)
	def tearDown(self):
		del self.lryeast
	def testARun(self):
		projectname = "Yeast_9_11_2013"
		proteinid = goDIR + "sgdID_" + projectname + ".mat"
                goidname = goDIR + "GOID_" + projectname + ".mat"
                go2genename = goDIR + "Gene2GO_" + projectname + ".mat"
                networkname = ppiDIR + "ppiNetwork_" + projectname + ".mat"
		self.lryeast.loadGOData(proteinid, goidname, go2genename, networkname)
		self.lryeast.run()

def suite():
        suite = makeSuite(TestLogRegression, 'test')
        return suite

if __name__ == '__main__':
        main(defaultTest='suite', argv=['-d'])
