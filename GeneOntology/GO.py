#/usr/bin/env python

from unittest import TestCase, makeSuite, main

from scipy.io import loadmat, savemat

class Ontology:
	def __init__(self):
		pass
	def __del__(self):
		pass
	def loadGOIDMatlab(self, fileName):
		goid = loadmat(fileName, squeeze_me=False, mat_dtype=True)
		goid = goid['GOID']
		return goid
	def loadProteinIDMatlab(self, fileName):
		proteinid = loadmat(fileName, squeeze_me=False, mat_dtype=True)
		proteinid = proteinid['proteinid']
		return proteinid
	def loadSGDIDMat(self, fileName):
                sgdid = loadmat(fileName, squeeze_me=False, mat_dtype=True)
                sgdid = sgdid["sgdID"]	
		return sgdid
	def loadEntrezGeneIDMat(self, fileName):
		entrezGene = loadmat(fileName, squeeze_me=False, mat_dtype=True)
		entrezGene = entrezGene['proteinid']
		#entrezGene = entrezGene[0]
		return entrezGene
	def loadGO2GeneMatlab(self, fileName):
		go2geneMatrix = loadmat(fileName, squeeze_me=False, mat_dtype=True)
		go2geneMatrix = go2geneMatrix['go2geneMatrix']
		if issparse(go2geneMatrix):
			return go2geneMatrix.todense()
		else:
			return go2geneMatrix
		
class TestOntology(TestCase):
	def setUp(self):
		self.goyeast = Ontology()
	def tearDown(self):
		del self.goyeast
	def testAYeast(self):
		pass

def suite():
	suite = makeSuite(TestOntology, 'test')
	return suite

if __name__ == '__main__':
	main(defaultTest='suite', argv=['-d'])
