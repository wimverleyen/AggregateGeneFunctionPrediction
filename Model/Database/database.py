#! /usr/bin/env python

from unittest import TestCase, makeSuite, main
import os

from numpy import mean, std
from scipy.stats import rankdata
import storm
from storm.databases.mysql import MySQL
from storm.locals import *

class Database:
	def __init__(self, dbName):
		"""
			Connect 2 database
		"""
		self.__uri = "mysql://:3306/" + dbName
		self.__database = create_database(self.__uri)
		self.__store = Store(self.__database)
	def __del__(self):
		del self.__uri
		self.__store.close()
	def dropTable(self, tableName):
		dropquery = "DROP TABLE IF EXISTS " + tableName
		print "dropquery= ", dropquery
		try:
			self.__store.execute(dropquery)
			self.__store.commit()
			self.__store.flush()
		except Exception as e:
			print "drop table exception"
			print e
	def createGOIDView(self, tableName, double, drop=False):
		"""
			Performance table
		"""
		if drop == True:
			try:
				self.dropTable(tableName)
			except Exception as e:
				print "drop query: exception caught"
				print e
		createquery = "CREATE TABLE IF NOT EXISTS " + tableName + " ( result_id INTEGER PRIMARY KEY, GOID INTEGER, CV INTEGER, "
		index = 0
		for column in double:
			createquery += column + " DOUBLE PRECISION"
			if index < len(double)-1:
				createquery += ", "
			index += 1
		createquery += " )"
		print "createquery= ", createquery
		try:
			self.__store.execute(createquery)	
			self.__store.flush()
		except Exception as e:
			print "create query: exception caught"
			print e
	def insertGOIDView(self, tableName, resultid, goid, cv, double):
		insertquery = "INSERT INTO " + tableName + " VALUES ( " + str(resultid) + ", "
		insertquery += str(goid) + ", "
		insertquery += str(cv) + ", "

		index = 0
		for value in double:
			insertquery += str(value)
			if index < len(double)-1:
				insertquery += ", "
			index += 1
		insertquery += " )"
		try:
			self.__store.execute(insertquery)	
			self.__store.flush()
		except Exception as e:
			print "insert query: exception caught"
			print insertquery
			print e
	def createProteinView(self, tableName, double, drop=False):
		"""
			Performance table
		"""
		if drop == True:
			try:
				self.dropTable(tableName)
			except Exception as e:
				print "drop query: exception caught"
				print e
		createquery = "CREATE TABLE IF NOT EXISTS " + tableName + " ( result_id INTEGER AUTO_INCREMENT PRIMARY KEY, GOID INTEGER, CV INTEGER, "
		index = 0
		for column in double:
			createquery += column + " DOUBLE PRECISION"
			if index < len(double)-1:
				createquery += ", "
			index += 1
		createquery += " )"
		print "createquery= ", createquery
		try:
			self.__store.execute(createquery)	
			self.__store.flush()
		except Exception as e:
			print "create query: exception caught"
			print e
			print createquery
	def insertProteinView(self, tableName, resultid, goid, cv, protein, label, score):
		for index in range(len(label)):
			double = []
			double.append(float(protein[index]))
			double.append(float(label[index]))
			double.append(float(score[index]))
			insertquery = "INSERT INTO " + tableName + " (GOID, CV, ProteinID, Label, Score) VALUES ( " #+ str(resultid) + ", "
			insertquery += str(goid) + ", "
			insertquery += str(cv) + ", "

			index = 0
			for value in double:
				insertquery += str(value)
				if index < len(double)-1:
					insertquery += ", "
				index += 1
			insertquery += " )"
			#print "insertequery= ", insertquery
			try:
				self.__store.execute(insertquery)	
				self.__store.flush()
			except Exception as e:
				print "insert query: exception caught"
				print e
				print insertquery
			del double
	
class TestDatabase(TestCase):
        def setUp(self):
                self.database = Database("funcpred")
        def tearDown(self):
                del self.database
        def testAGOIDView(self):
		pass

def suite():
        suite = makeSuite(TestDatabase, 'test')
        return suite

if __name__ == '__main__':
        main(defaultTest='suite', argv=['-d'])

