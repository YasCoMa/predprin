import luigi

import numpy as np
import os

from core.tasks.pygosemsim.similarity import *
from core.tasks.pygosemsim.term_set import *
import functools

from core.tasks.sprint_tasks import GenerateFileScoresTask
from ..time_counter import TimeTaskMixin

class GetCalculatedFeature():

	def get_info_pair(self, pair, folder, information):
		info={}
		for i in information:
			info[i] = np.load(folder+i+".npy", allow_pickle=True)

		features=[]
		p=pair.split(",")
		ps=[p[0], p[1]]
		for p in ps:
			feature=[]
			for i in information:
				try:
					feature.append(eval("info['"+i+"'].item().get('"+p+"')"))
				except:
					if(i=='ko'):
						feature.append("None")
					else:
						feature.append(['None'])

			features.append(feature)

		return features

	def calc_write_go(self, pair, graph_go, folder):
		information=['go_cc','go_bp','go_mf']
		features=self.get_info_pair(pair[0]+","+pair[1], folder, information)
		
		components=['cc','bp','mf']
		c=0
		for branch in components:
			l1=features[0][c]
			l2=features[1][c]

			m=0
			if(l1!=None and l1!='None' and l2!=None and l2!='None'):
				sf = functools.partial(sim_func, graph_go, pekar)
				m = sim_bma(l1, l2, sf)
				if(m==None):
					m=0.0
				
			with open(folder+"dataset_ppi.txt", "a") as g:
				g.write(" "+str(m))
			c+=1

	def calc_write_pfam(self, pair, folder):
		information=['pfam']
		features=self.get_info_pair(pair[0]+","+pair[1], folder, information)
		
		pfam_ids1=features[0][0]
		pfam_ids2=features[1][0]

		pfam_interaction = np.load(folder+"pfam_interaction.npy", allow_pickle=True )
		pfam_=0.0
		if(pfam_ids1!=None and pfam_ids1!='None' and pfam_ids2!=None and pfam_ids2!='None'):
			for pfam1 in pfam_ids1:
				for pfam2 in pfam_ids2:
					if( (pfam1+"-"+pfam2) in pfam_interaction.item().keys()):
						if(pfam_interaction.item().get(pfam1+"-"+pfam2)=="1"):
							pfam_=1.0
							break

					if( (pfam2+"-"+pfam1) in pfam_interaction.item().keys()):
						if(pfam_interaction.item().get(pfam2+"-"+pfam1)=="1"):
							pfam_=1.0
							break

		with open(folder+"dataset_ppi.txt", "a") as g:
			g.write(" "+str(pfam_))

	def calc_write_pathway(self, pair, folder):
		information=['ko']
		features=self.get_info_pair(pair[0]+","+pair[1], folder, information)

		kos = np.load(folder+"kos.npy", allow_pickle=True)
		ko1=features[0][0]
		ko2=features[1][0]
		maps_ko1=[]
		maps_ko2=[]
		try:
			maps_ko1=kos.item().get(ko1)
			if(maps_ko1==None):
				maps_ko1=[]
		except:
			pass

		try:
			maps_ko2=kos.item().get(ko2)
			if(maps_ko2==None):
				maps_ko2=[]
		except:
			pass

		in_common=list(np.intersect1d(maps_ko1,maps_ko2))
		kegg=0.0
		if(len(in_common)>1):
			kegg=1.0

		with open(folder+"dataset_ppi.txt", "a") as g:
			g.write(" "+str(kegg))

	def calc_write_sequence(self, pair, folder):
		seq_dict = np.load(folder+"seq_dict.npy", allow_pickle=True)

		score=0.0
		key = pair[0]+"-"+pair[1]
		if( key in seq_dict.item().keys()):
			score=seq_dict.item().get(key)
		
		key = pair[1]+"-"+pair[0]
		if( key in seq_dict.item().keys()):
			score=seq_dict.item().get(key)
		if(float(score) > 20.0):
			score=1.0

		with open(folder+"dataset_ppi.txt", "a") as g:
			g.write(" "+str(score))
			
class InitializeVariablesTask(luigi.Task, TimeTaskMixin):
	folder = luigi.Parameter()
	prefix = luigi.Parameter()

	def run(self):
		print("Running initilizing variables")
		
		seq_dict ={}
		f=open(self.folder+"scores.txt","r")
		for l in f:
			l=l.replace("\n","").split(" ")
			seq_dict[l[0]+"-"+l[1]]=float(l[2])
		f.close()
		np.save(self.folder+"seq_dict", seq_dict)

		pfam_interaction={}
		f=open("iddi_1-0.txt","r")
		for line in f:
			l=line.split("\t")
			pfam_interaction[l[0]+"-"+l[1]]=l[26]
		f.close()
		np.save(self.folder+"pfam_interaction", pfam_interaction)

		kos={}
		f=open("map_kegg.txt","r")
		for line in f:
			l=str(line).replace("\n","").split("\t")

			if(l[0].find("path:map")!=-1):
				ko_=l[1].split(":")[1]
				pathway=l[0].split(":")[1]
				if(not ko_ in kos.keys()):
					kos[ko_]=[]
				kos[ko_].append(pathway)

		f.close()
		np.save(self.folder+"kos", kos)

		go_cc={}
		go_bp={}
		go_mf={}
		pfam={}
		ko={}
		
		proteins=[]
		f=open(self.folder+self.prefix+"dataset.txt","r")
		for line in f:
			l=line.split("\t")
			p1=l[0]
			p2=l[1]
			
			if(not p1 in proteins):
				proteins.append(p1)
			
			if(not p2 in proteins):
				proteins.append(p2)

		f.close()

		for p in proteins:
			f=open("annotation_data/"+p+".tsv","r")
			for line in f:
				l=line.replace("\n","").split("\t")
				
				go_cc[l[0]]=l[1].split(" ")
				go_bp[l[0]]=l[2].split(" ")
				go_mf[l[0]]=l[3].split(" ")
				pfam[l[0]]=l[5].split(" ")
				ko[l[0]]=l[4]

				break
			f.close()

		np.save(self.folder+"go_cc", go_cc)
		np.save(self.folder+"go_bp", go_bp)
		np.save(self.folder+"go_mf", go_mf)
		np.save(self.folder+"pfam", pfam)
		np.save(self.folder+"ko", ko)

		self.output().open("w").close()

	def requires(self):
		return GenerateFileScoresTask(self.folder, self.prefix)

	def output(self):
		return luigi.LocalTarget(self.folder+"log/initilizing_variables.txt")
