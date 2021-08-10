import numpy as np

import matplotlib.pyplot as plt

from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.ensemble import AdaBoostClassifier
h = .02  # step size in the mesh

from sklearn.metrics import roc_curve, auc
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, classification_report, confusion_matrix

import luigi
from ..time_counter import TimeTaskMixin
from joblib import dump, load
import os

class ClassifyEvaluateTask(luigi.Task, TimeTaskMixin):
	folder = luigi.Parameter()
	prefix = luigi.Parameter()
	mode = luigi.Parameter()
	model = luigi.Parameter()

	def run(self):
		print("Running visualization")
		if(self.mode=='train'):
			yield PlotResults(self.folder)
		self.output().open("w").close()
		
	def requires(self):
		if(self.mode=='train'):
			return ClassifyCrossValidateTask(self.folder)
		if(self.mode=='test'):
			return ApplyModelForTestTask(self.folder, self.prefix, self.model)

	def output(self):
		return luigi.LocalTarget( os.path.join(self.folder+"log", "classification_evaluation.txt") )

class ApplyModelForTestTask(luigi.Task, TimeTaskMixin):
	folder = luigi.Parameter()
	prefix = luigi.Parameter()
	model = luigi.Parameter()

	def run(self):
		print("Running application of models for testing")
		pairs=[]
		f=open(self.folder+self.prefix+"dataset.txt","r")
		for line in f:
			if(line!=""):
				l=line.replace("\n","").split("\t")
				pairs.append([l[0], l[1]])
		f.close()
		
		clf = load(self.model) 
		X = np.load(self.folder+"x.npy", allow_pickle=True)
		predictions = clf.predict(X)
		
		c=0
		f=open(self.folder+"predictions.tsv","w")
		for p in predictions:
			f.write("%s\t%s\t%.3f\n" %(pairs[c][0], pairs[c][1], p) )
			c+=1
		f.close()
		
		self.output().open("w").close()
		
	def requires(self):
		yield ProcessDataForClassification(self.folder)

	def output(self):
		return luigi.LocalTarget( os.path.join(self.folder+"log", "model_application.txt") )

class PlotResults(luigi.Task, TimeTaskMixin):
	folder = luigi.Parameter()

	def run(self):
		print("Running visualization")

		predictions = np.load(self.folder+"predictions.npy", allow_pickle=True)
		y_test = np.load(self.folder+"y_test.npy", allow_pickle=True)

		X = np.load(self.folder+"x.npy", allow_pickle=True)
		aux_y = np.load(self.folder+"y.npy", allow_pickle=True)

		features=['sequence','go-cc','go-bp', 'go-mf','pfam','kegg','predrep']
		for i in range(X.shape[1]+1):
			if(i==X.shape[1]):
				y=y_test
				x=predictions
			else:
				y=aux_y
				x=X[:,i]

			false_positive_rate, true_positive_rate, thresholds = roc_curve(y, x)
			roc_auc = auc(false_positive_rate, true_positive_rate)
			f = plt.figure()
			plt.title('ROC Curve - Combination ('+features[i]+')')
			plt.plot(false_positive_rate, true_positive_rate, 'b', label='AUC = %0.2f'% roc_auc)
			plt.legend(loc='lower right')
			plt.plot([0,1],[0,1],'r--')
			plt.xlim([-0.1,1.2])
			plt.ylim([-0.1,1.2])
			plt.ylabel('True Positive Rate')
			plt.xlabel('False Positive Rate')
			f.savefig(self.folder+features[i]+'-roc_auc.png')

		f = plt.figure()
		plt.title('ROC Curve')
		ax = plt.subplot(111)
		for i in range(X.shape[1]+1):
			if(i==X.shape[1]):
				y=y_test
				x=predictions
			else:
				y=aux_y
				x=X[:,i]

			false_positive_rate, true_positive_rate, thresholds = roc_curve(y, x)
			roc_auc = auc(false_positive_rate, true_positive_rate)
			plt.plot(false_positive_rate, true_positive_rate, 'b', label=features[i]+' - AUC = %0.2f'% roc_auc)
		chartBox = ax.get_position()
		ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*0.6, chartBox.height])
		ax.legend(loc='upper center', bbox_to_anchor=(1.45, 0.8), shadow=True, ncol=1)
		plt.xlim([-0.1,1.2])
		plt.ylim([-0.1,1.2])
		plt.ylabel('True Positive Rate')
		plt.xlabel('False Positive Rate')
		f.savefig(self.folder+'all-roc_auc.png')
		
		self.output().open("w").close()

	def output(self):
		return luigi.LocalTarget( os.path.join(self.folder+"log", "plot_results.txt") )

class ClassifyCrossValidateTask(luigi.Task, TimeTaskMixin):
	folder = luigi.Parameter()

	def run(self):
		print("Running classification and evaluation")
		yield DoClassificationLogTask(self.folder)
		yield DoCrossValidationTask(self.folder)
		self.output().open("w").close()
		
	def requires(self):
		yield ProcessDataForClassification(self.folder)

	def output(self):
		return luigi.LocalTarget( os.path.join(self.folder+"log", "classification_evaluation.txt") )

class DoClassificationLogTask(luigi.Task, TimeTaskMixin):
	folder = luigi.Parameter()

	def run(self):

		X_train = np.load(self.folder+"x_train.npy", allow_pickle=True)
		X_test = np.load(self.folder+"x_test.npy", allow_pickle=True)
		y_train = np.load(self.folder+"y_train.npy", allow_pickle=True)
		y_test = np.load(self.folder+"y_test.npy", allow_pickle=True)
		
		# consider changing to stacking ensemble https://machinelearningmastery.com/stacking-ensemble-machine-learning-with-python/
		name = "PredRep"
		clf = AdaBoostClassifier()
		clf.fit(X_train, y_train)
		dump(clf, self.folder+'model_trained.joblib')
		predictions = clf.predict(X_test)

		threshold=0.8

		np.save(self.folder+"predictions", predictions)

		f=open(self.folder+"evaluation_log-"+str(threshold)+".txt","w")
		f.close()

		tn=0.0
		fn=0.0
		tp=0.0
		fp=0.0
		
		p=0.0
		n=0.0

		a=0
		for pr in predictions:
			if(y_test[a]==1):
				p+=1

				if(float(pr)>=threshold):
					tp+=1
				else:
					fn+=1
			elif(y_test[a]==0):
				n+=1
				
				if(float(pr)<threshold):
					tn+=1
				else:
					fp+=1
			a+=1

		with open(self.folder+"evaluation_log-"+str(threshold)+".txt","a") as fg:
			fg.write("Results :\n")
			fg.write("	Positives: "+str(p)+"\n")
			fg.write("	Negatives: "+str(n)+"\n")
		
			fg.write("	Confusion Matrix:\n")
			fg.write("		tp ("+str(tp)+") | fn ("+str(fn)+")\n")
			fg.write("		fp ("+str(fp)+") | tn ("+str(tn)+")\n")

			fg.write("	Saving ROC Curve plot...\n")
			fg.write("	Metrics:\n")
			fg.write('		Accuracy:'+str(accuracy_score(y_test, predictions))+"\n")
			fg.write('		F1 score:'+str(f1_score(y_test, predictions))+"\n")
			fg.write('		Recall:'+str(recall_score(y_test, predictions))+"\n")
			fg.write('		Precision:'+str(precision_score(y_test, predictions))+"\n")
			fg.write('		\nClassification report:\n'+str(classification_report(y_test,predictions)))
			fg.write('		\nConfusion matrix:\n'+str(confusion_matrix(y_test, predictions))+"\n")

		self.output().open("w").close()

	def output(self):
		return luigi.LocalTarget( os.path.join(self.folder+"log", "do_classification_log.txt") )

class DoCrossValidationTask(luigi.Task, TimeTaskMixin):
	folder = luigi.Parameter()

	def run(self):
		X = np.load(self.folder+"x.npy", allow_pickle=True)
		y = np.load(self.folder+"y.npy", allow_pickle=True)

		clf = AdaBoostClassifier()

		y_pred = cross_val_predict(clf, X, y, cv=10)
		
		np.save(self.folder+"final_score", y_pred)

		f1=cross_val_score(clf, X, y, scoring='f1', cv=10)
		precision=cross_val_score(clf, X, y, scoring='precision', cv=10)
		recall=cross_val_score(clf, X, y, scoring='recall', cv=10) 
		accuracy=cross_val_score(clf, X, y, scoring='accuracy', cv=10)

		f=open(self.folder+"result_cross-validation.txt","w")
		f.close()

		with open(self.folder+"result_cross-validation.txt", "a") as g:
			g.write("f1;precision;recall;accuracy\n")

		for i in range(len(f1)):
			with open(self.folder+"result_cross-validation.txt", "a") as g:
				g.write(str(f1[i])+";"+str(precision[i])+";"+str(recall[i])+";"+str(accuracy[i])+"\n")
		with open(self.folder+"result_cross-validation.txt", "a") as g:
			g.write("\n")

		self.output().open("w").close()

	def output(self):
		return luigi.LocalTarget( os.path.join(self.folder+"log", "do_cross_validation.txt") )

class ProcessDataForClassification(luigi.Task, TimeTaskMixin):
	folder = luigi.Parameter()

	def run(self):
		X=[]
		y=[]
		f=open(self.folder+"dataset_ppi.txt","r")
		for line in f:
			if(line!=""):
				l=line.replace("\n","").split(" ")
				if(l[-1]=="+1"):
					y.append(1)
				else:
					y.append(0)

				x_temp=[]
				for feature in l[:-1]:
					if(feature!=""):
						x_temp.append(float(feature))
				X.append(x_temp)
		f.close()

		X=np.array(X)
		y=np.array(y)

		#X = StandardScaler().fit_transform(X)

		X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.4, random_state=42)
		np.save(self.folder+"x_train", X_train)
		np.save(self.folder+"x_test", X_test)
		np.save(self.folder+"y_train", y_train)
		np.save(self.folder+"y_test", y_test)
		np.save(self.folder+"x", X)
		np.save(self.folder+"y", y)

		self.output().open("w").close()

	def output(self):
		return luigi.LocalTarget( os.path.join(self.folder+"log", "process_data.txt") )
