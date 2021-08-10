#from sprint_tasks import SprintTask

import numpy as np

from core.tasks.pygosemsim.graph import *
from core.tasks.pygosemsim.similarity import *

import luigi
import os
from ..time_counter import TimeTaskMixin
from ..utils import Utils_search

from core.tasks.calc_features_tasks import GetCalculatedFeature, InitializeVariablesTask

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
		return luigi.LocalTarget( os.path.join(self.folder+"log", "generate_features.txt") )

class GetFeatureVectorTask(luigi.Task, TimeTaskMixin):
	folder = luigi.Parameter()
	prefix = luigi.Parameter()

	def run(self):
		print("Running get feature vector")

		G = from_resource("go")
		precalc_lower_bounds(G)

		#u = Utils_search()
		pairs={}
		count=0
		data = ['positive','false']
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

		u=GetCalculatedFeature()
		graph_go = G
		for p in pairs.keys():
			class_=pairs[p]
			pair=p.split(",")

			metrics=u.calc_write_go(pair, graph_go, self.folder)
			metrics.append(u.calc_write_pfam(pair, self.folder))
			metrics.append(u.calc_write_pathway(pair, self.folder))
			metrics.append(u.calc_write_sequence(pair, self.folder))

			with open(self.folder+"dataset_ppi.txt", "a") as g:
				g.write((" ".join(metrics))+" "+class_+"\n")

			#infos=[ ['go_cc','go_bp','go_mf'], ['pfam'], ['ko'] ]
			#for information in infos:
			#	os.system("rm "+self.folder+pair[0]+"_"+pair[1]+"-"+(("-").join(information))+"-get_info_pair.npy")

		self.output().open("w").close()

	def output(self):
		return luigi.LocalTarget( os.path.join(self.folder+"log", "get_feature_vector.txt") )
