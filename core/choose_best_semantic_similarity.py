from pygosemsim import graph
from pygosemsim import annotation
from pygosemsim import similarity
from pygosemsim import term_set
import functools

import numpy as np
import os
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, auc

class Evaluation_semantic_similarities_go:
	def __init__(self, folder, prefix):
		self.main_folder=folder
		self.prefix=prefix
		features=['go-bp', 'go-cc', 'go-mf']
		sem_sim_methods=['resnik','lin','jiang','pekar','wang']
		combo=[]
		for f in features:
			for m in sem_sim_methods:
				combo.append(f+"_"+m)
		self.combo=combo

	def generate_features(self):
		folder=self.main_folder
		prefix=self.prefix

		G = graph.from_resource("go")
		similarity.precalc_lower_bounds(G)
		annot = annotation.from_resource("goa_human")

		branches=['biological_process', 'cellular_component', 'molecular_function' ]
		
		features=['go-bp', 'go-cc', 'go-mf']
		sem_sim_methods=['resnik','lin','jiang','pekar','wang']
		combo=[]
		for f in features:
			for m in sem_sim_methods:
				combo.append(f+"_"+m)

		f=open("dataset_report_go_semsim.tsv", "w")
		f.close()

		c=0
		ds=['false','positive']
		for d in ds:
			class_=0
			if(d=='positive'):
				class_=1

			f=open(folder+prefix+d+".txt","r")
			for line in f:
				l=line.split("\t")
				p1=l[0]
				p2=l[1]

				metrics=[]
				if(p1 in annot.keys() and p2 in annot.keys()):
					fp1 = annot[p1]["annotation"].keys()
					fp1_={}
					for b in branches:
						fp1_[b]=[]

					for go in fp1:
						fp1_[G.node[go]['namespace']].append(go)

					fp2 = annot[p2]["annotation"].keys()
					fp2_={}
					for b in branches:
						fp2_[b]=[]

					for go in fp2:
						fp2_[G.node[go]['namespace']].append(go)

					for b in branches:
						for m in sem_sim_methods:
							sf = functools.partial(term_set.sim_func, G, eval('similarity.'+m))
							metrics.append(str(term_set.sim_bma(fp1_[b], fp2_[b], sf)))
				else:
					for b in branches:
						for m in sem_sim_methods:
							metrics.append("0.0")

				with open(folder+"dataset_report_go_semsim.tsv", "a") as fg:
					fg.write(("\t".join(metrics))+"\t"+str(class_)+"\n")

				c+=1
				print(c)
		f.close()

	def evaluation(self):
		folder=self.main_folder
		combo=self.combo

		X=[]
		y=[]
		f=open(folder+"dataset_report_go_semsim.tsv","r")
		for line in f:
			l=line.replace("\n","").split("\t")
			x_temp=[]
			for x in l:
				if(x=='None'):
					x=0.0
				x_temp.append(float(x))
			X.append(x_temp[:-1])

			y.append(x_temp[-1])
		f.close()

		X=np.array(X)
		y=np.array(y)

		os.system("rm -rf "+folder+"charts_gosemsim_analysis")
		os.system("mkdir "+folder+"charts_gosemsim_analysis")
		X = StandardScaler().fit_transform(X)

		f = plt.figure()
		plt.title('ROC Curve')
		for i in range(X.shape[1]):
			self.plot_auc_curve(y, X[:,i], combo[i], folder+"charts_gosemsim_analysis/")
			false_positive_rate, true_positive_rate, thresholds = roc_curve(y, X[:,i])
			roc_auc = auc(false_positive_rate, true_positive_rate)

			plt.plot(false_positive_rate, true_positive_rate, 'b', label=combo[i]+' - AUC = %0.2f'% roc_auc)

		plt.xlim([-0.1,1.2])
		plt.ylim([-0.1,1.2])
		plt.ylabel('True Positive Rate')
		plt.xlabel('False Positive Rate')

		f.savefig(folder+'all-roc_auc.png')

	def plot_auc_curve(self, ytest, preds, id_, folder):
		false_positive_rate, true_positive_rate, thresholds = roc_curve(ytest, preds)
		roc_auc = auc(false_positive_rate, true_positive_rate)
		
		f = plt.figure()
		plt.title('ROC Curve - Combination ('+id_+')')
		plt.plot(false_positive_rate, true_positive_rate, 'b', label='AUC = %0.2f'% roc_auc)
		plt.legend(loc='lower right')
		plt.plot([0,1],[0,1],'r--')
		plt.xlim([-0.1,1.2])
		plt.ylim([-0.1,1.2])
		plt.ylabel('True Positive Rate')
		plt.xlabel('False Positive Rate')

		f.savefig(folder+id_.lower()+'-roc_auc.png')

	def run(self):
		#print("Calculating metrics and generating features...")
		#self.generate_features()
		print("Doing evaluation for isolated methods ...")
		self.evaluation()

import sys
a=Evaluation_semantic_similarities_go(sys.argv[1],sys.argv[2])
a.run()