#! /usr/bin/env python

from unittest import TestCase, makeSuite, main

from random import uniform
from numpy import nan, std, mean, arange, sort, argsort
from scipy.stats import rankdata

from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_curve
from sklearn.metrics import auc, auc_score, f1_score

class Performance:
	def __init__(self, labels, scores, projectName=""):
		self.__labels = labels
		self.__scores = scores
		self.__projectname = projectName
		self.__f1 = []
		
	def __del__(self):
		del self.__labels
		del self.__scores
	def AUPRGillis(self):
		ranks = rankdata(self.__scores)
		Npos = len([label for label in self.__labels if label > 0])
		Nneg = len([label for label in self.__labels if label <= 0])
		ranksidx = argsort(ranks)

		ranksum = 0.0
		count = 1
		for index in ranksidx:
			if self.__labels[index] == 1:
				ranksum += count/float(ranks[index])
				count += 1
		if Npos > 0:
			auc = ranksum/float(Npos)
		else:
			auc = nan
		return auc
	def AUROC(self):
		try:
			self.__fpr, self.__tpr, thresholds = roc_curve(self.__labels, self.__scores)
		except Exception as e:
			print "roc_curve exception"
			print e
			return nan
		try:
			self.__rocarea = auc(self.__fpr, self.__tpr)
		except Exception as e:
			print "roc auc exception"
			print e
			return nan
		return self.__rocarea
	def AUROCGillis(self):
		print len(self.__labels)
		print len(self.__scores)
		ranks = rankdata(self.__scores)
		Npos = len([label for label in self.__labels if label > 0])
		Nneg = len([label for label in self.__labels if label <= 0])

		ranksum = 0.0
		index = 0
		for rank in ranks:
			if self.__labels[index] == 1:
				ranksum += rank
			index += 1
		value = ranksum - ((Npos * (Npos + 1))/float(2))
		if Npos > 0:
			value = value/float(Npos * Nneg)
		else:
			value = 0.5
		auc = 1 - value
		return auc
	def AUROCScore(self):
		try:
			self.__rocarea = auc_score(self.__labels, self.__scores)
		except Exception as e:
			print "roc_curve exception"
			print e
			return nan
		return self.__rocarea
	def Fmax(self):
		ranks = rankdata(self.__scores)
		ranksidx = argsort(ranks)
		Npos = len([label for label in self.__labels if label > 0])

		rankedscores = zeros(len(self.__labels))
		rankedlabels = zeros(len(self.__labels))
		i = 0
		for index in ranksidx:
			rankedscores[i] = self.__scores[index]
			rankedlabels[i] = self.__labels[index]
			i += 1

		tp = 0
		fp = 0
		index = 0
		for score in rankedscores:
			if rankedlabels[index] == 1:
				tp += 1
			else:
				fp += 1
			prec = tp/float(tp+fp)
			if Npos > 0:
				rec = tp/float(Npos)
			else:
				rec = 0.0
			if float(prec+rec) > 0.0:
				f = (2*prec*rec)/float(prec+rec)
			else:
				f = 0.0
			self.__f1.append(f)
			index += 1
			if tp == Npos:
				break
		return max(self.__f1)
	def F1(self):
		self.__prec, self.__recall, thresholds = precision_recall_curve(self.__labels, self.__scores)
		index = 0
		for value in thresholds:
			f1 = self.__prec[index] * self.__recall[index]/float(self.__prec[index] + \
				self.__recall[index])
			print f1
			self.__f1.append(float)
			index += 1
		return self.__f1
	
class TestPerformance(TestCase):
        def setUp(self):
                self.performance = Performance([], [], projectName="YeastPPI_9_11_2013")
        def tearDown(self):
                del self.performance
        def testAPerformance(self):
		pass
	
def suite():
        suite = makeSuite(TestPerformance, 'test')
        return suite

if __name__ == '__main__':
        main(defaultTest='suite', argv=['-d'])
