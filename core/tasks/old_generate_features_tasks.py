#from sprint_tasks import SprintTask

import numpy as np

from pygosemsim import graph
from pygosemsim import similarity

import luigi
import os
from ..time_counter import TimeTaskMixin

from core.tasks.calc_features_tasks import GetCalculatedFeatureGoTask, GetCalculatedFeaturePfamTask, GetCalculatedFeaturePathwayTask, GetCalculatedFeatureSequenceTask, InitializeVariablesTask

class GenerateFeaturesTask(luigi.Task, TimeTaskMixin):
	folder = luigi.Parameter()
	prefix = luigi.Parameter()

	def run(self):
		print("Running protein annotation")
		yield GetFeatureVectorTask(self.folder, self.prefix)
		self.output().open("w").close()
		
	def requires(self):
		yield InitializeVariablesTask(self.folder, self.prefix)

	def output(self):
		return luigi.LocalTarget(self.folder+"log/generate_features.txt")

class GetFeatureVectorTask(luigi.Task, TimeTaskMixin):
	folder = luigi.Parameter()
	prefix = luigi.Parameter()

	def run(self):
		print("Running get feature vector")
		self.output().open("w").close()
		
	def requires(self):
		# colocar no require de retrieveInfoPair
		#self.initialize_protein_features(folder)
		#self.initiatize_score_sequence(folder)

		G = graph.from_resource("go")
		similarity.precalc_lower_bounds(G)

		#u = Utils_search()
		pairs={}
		count=0
		data = ['false', 'positive']
		for d in data:
			class_="-1"
			if(d=='positive'):
				class_="+1"
			f=open(self.folder+self.prefix+d+".txt","r")
			for line in f:
				l=line.split("\t")
				p1=l[0]
				p2=l[1]
				
				pairs[p1+","+p2]=class_
				#features = self.unite_features(p1, p2, d, folder)
				#u.write_pair_features(features, folder, class_)

			f.close()

		for p in pairs.keys():
			pr=p.split(",")
			yield GetPairFeatureWrite(pr, pairs[p], self.folder, G)

	def output(self):
		return luigi.LocalTarget(self.folder+"log/get_feature_vector.txt")

class GetPairFeatureWrite(luigi.Task, TimeTaskMixin):
	pair = luigi.Parameter()
	class_ = luigi.Parameter()
	folder = luigi.Parameter()
	graph_go = luigi.Parameter()

	def run(self):
		print("Running calc features")
		
		with open(self.folder+"dataset_ppi.txt", "a") as g:
			g.write(" "+self.class_+"\n")

		infos=[ ['go_cc','go_bp','go_mf'], ['pfam'], ['ko'] ]
		for information in infos:
			os.system("rm "+self.folder+self.pair[0]+"_"+self.pair[1]+"-"+(("-").join(information))+"-get_info_pair.npy")

		self.output().open("w").close()
		
	def requires(self):
		yield GetCalculatedFeaturePathwayTask(self.pair, self.folder)
		yield GetCalculatedFeaturePfamTask(self.pair, self.folder)
		yield GetCalculatedFeatureGoTask(self.pair, self.graph_go, self.folder)
		yield GetCalculatedFeatureSequenceTask(self.pair, self.folder)
		
	def output(self):
		return luigi.LocalTarget(self.folder+"log/unite_features_write.txt")