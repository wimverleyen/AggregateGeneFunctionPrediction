#! /usr/bin/env python

from unittest import TestCase, makeSuite, main

from os import system
from math import floor
from numpy import where, tile, transpose, dot, array, int32, equal, asarray, ones_like, subtract, reshape, isnan
from numpy.random import shuffle, random_sample
from numpy.random import permutation, seed
from scipy.sparse import issparse, lil_matrix

biogridDIR = "/home/wverleye/work/scripts/Data/BioGRID/"


class NeighborVoting:
	def __init__(self, dbName, goidTable="NeighborVoting_Mousefunc_GOID", proteinTable="NeighborVoting_Mousefunc_Protein"):
		seed(66)
		self.__database = Database(dbName=dbName)
		self.__goidtable = goidTable
		self.__proteintable = proteinTable

		self.__proteins = []
		self.__numproteins = 0
		self.__goid = []

		self.__go2gene = []
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
                annotations = self.__go2gene[:,goid]
                return annotations
		
	def _fit(self, annotation):
                sumin = dot(self.__network, annotation)
                sumall = tile(sum(self.__network), (1, 1))
                sumall = sumall.transpose()
                scores = sumin/sumall

                del sumin
                del sumall

                scores = subtract(ones_like(scores), scores)
                return scores

	def runMousefunc(self, nFold=3):
		"""
			CV: -1 => total model (no cv)
			CV: nFold => mean metric over cv
		"""
		self.__database.createGOIDView(self.__goidtable, \
				double=["AUROC", "AUPR", "Fmax"], drop=True)
		self.__database.createProteinView(self.__proteintable, \
				double=["ProteinID", "Label", "Score"], drop=True)

		# Get labels
		test = 0
		pp = permutation(self.__numproteins)
		resultid = 0
		for goid in self.__goid:
			print "____________ GOID %d ____________" % goid
			goidindex = where(self.__goid==goid)
                        goidindex = int(goidindex[0])
			annotations = self.selectAnnotatedProteinsMousefunc(goidindex)
			annotations = asarray(annotations).ravel()
			annotations = reshape(annotations, (annotations.shape[0], 1))

			print "0s=", len([x for x in annotations if x == 0])
			print "1s=", len([x for x in annotations if x == 1])
			print "-1s=", len([x for x in annotations if x == -1])

			scores = self._fit(annotations)

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
			del per

			resultid += 1

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

				labeltmp = asarray(labeltmp).ravel()
				labeltmp = reshape(labeltmp, (labeltmp.shape[0], 1))
				scores = self._fit(labeltmp)

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
			self.__database.insertGOIDView(self.__goidtable, resultid, goid[0], nFold, [roc_mean, pr_mean, fmax_mean])
			resultid += 1
			test += 1
	

class TestNeighborVoting(TestCase):
	def setUp(self):
		projectname = "HumanPPIBioGRID_8_20_2015"
		goidtable = "NeighborVoting_" + projectname + "_GOID"
		proteintable = "NeighborVoting_" + projectname + "_Protein"
		self.nvhuman = NeighborVoting(dbName="funcpred", goidTable=goidtable, \
				proteinTable=proteintable)

	def tearDown(self):
		del self.nvhuman

	def testARun(self):
		projectname = "Human_8_20_2015"
		proteinid = goDIR + "entrezGeneID_" + projectname + ".mat"
		goid = goDIR + "GOID_" + projectname + ".mat"
                function = goDIR + "Gene2GO_" + projectname + ".mat"
                networkname = biogridDIR + "ppiNetwork_" + projectname + ".mat"
		self.nvhuman.loadGOData(proteinid, goid, function, networkname)
		self.nvhuman.runMousefunc()


def suite():
        suite = makeSuite(TestNeighborVoting, 'test')
        return suite

if __name__ == '__main__':
        main(defaultTest='suite', argv=['-d'])
