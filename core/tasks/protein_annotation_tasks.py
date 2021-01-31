from ..utils import Utils_search, PreProcessing

from rdflib import Graph
import luigi
import datetime, os

from core.tasks.pygosemsim.download import *
from core.tasks.pygosemsim.graph import *

from ..time_counter import TimeTaskMixin

class ProteinAnnotationTask(luigi.Task, TimeTaskMixin):
	folder = luigi.Parameter()
	prefix = luigi.Parameter()

	def run(self):
		print("Running protein annotation")
		self.output().open("w").close()
		
	def requires(self):
		return GetInfoPairsTask(self.folder, self.prefix)

	def output(self):
		return luigi.LocalTarget(self.folder+"log/protein_annotation.txt")

class GetInfoPairsTask(luigi.Task, TimeTaskMixin):
	folder = luigi.Parameter()
	prefix = luigi.Parameter()

	def run(self):
		print("Running protein annotation")

		folder=self.folder
		prefix=self.prefix

		try:
			G = from_resource("go")
		except:
			obo("go")
			G = from_resource("go")
		
		u = Utils_search()
		list_annotated=u.get_annotated_proteins(folder)

		new_proteins=[]
		f=open(folder+prefix+"dataset.txt","r")
		for line in f:
			l=line.split("\t")
			p1=l[0]
			p2=l[1]
			
			if(not p1 in list_annotated and not p1 in new_proteins):
				new_proteins.append(p1)
			
			if(not p2 in list_annotated and not p2 in new_proteins):
				new_proteins.append(p2)

		f.close()

		folder = self.folder
		prefix = self.prefix
		graph_go = G
		for p in new_proteins:
			protein = p

			features=[]
			go_cc_ids=[]
			go_mf_ids=[]
			go_bp_ids=[]
			ko_ids=""
			pfam_ids=[]

			try:
				g=Graph()
				g.parse("rdf_data/"+protein+".rdf", format="xml")
				
				results=g.query("""
				SELECT distinct ?o ?c
				WHERE{
						?s <http://purl.uniprot.org/core/classifiedWith> ?o .
						?s <http://www.w3.org/2000/01/rdf-schema#seeAlso> ?c .
				}
				""")
				namespace=""
				for row in results.result:
					info=row[0]

					if(info.find("obo/GO_")!=-1):
						go=info.replace("http://purl.obolibrary.org/obo/","")
						go_term=go.replace("_", ":")
						
						#namespace = str(self.get_namespace_from_go(go_term.replace(":","_")))
						if(go_term in graph_go.nodes):
							namespace = graph_go.nodes[go_term]['namespace']
							#print(go_term, namespace)
							if(namespace=="cellular_component"):
								if(not go_term in go_cc_ids):
									go_cc_ids.append(go_term)
							if(namespace=="molecular_function"):
								if(not go_term in go_mf_ids):
									go_mf_ids.append(go_term)
							if(namespace=="biological_process"):
								if(not go_term in go_bp_ids):
									go_bp_ids.append(go_term)	
						#go_ids.append(str(go_term[0]))

					info2=row[1]
					if(info2.find("pfam")!=-1 and not(info2.find("/supfam")!=-1)):
						pfam=info2
						pfam=pfam.replace("http://purl.uniprot.org/pfam/","")
						if(not pfam in pfam_ids):
							pfam_ids.append(pfam)

					if(info2.find("ko")!=-1):
						ko=info2
						ko_ids=ko.replace("http://purl.uniprot.org/ko/","")

				if(len(go_cc_ids)==0):
					go_cc_ids.append("None")
				if(len(go_mf_ids)==0):
					go_mf_ids.append("None")
				if(len(go_bp_ids)==0):
					go_bp_ids.append("None")
				if(len(pfam_ids)==0):
					pfam_ids.append("None")
				if(ko_ids==""):
					ko_ids="None"

			except:
				go_cc_ids.append("None")
				go_mf_ids.append("None")
				go_bp_ids.append("None")
				pfam_ids.append("None")
				ko_ids="None"

			features.append(go_cc_ids)
			features.append(go_mf_ids)
			features.append(go_bp_ids)
			features.append(ko_ids)
			features.append(pfam_ids)

			txt=protein+"\t"
			b=0
			for feature in features:
				if(not(str(type(feature))=="<class 'str'>")):
					a=0
					for fea in feature:
						txt+=str(fea)
						if(a!=len(feature)-1):
							txt+=" "
						a+=1

					if(b!=len(features)-1):
						txt+="\t"
					b+=1
				else:
					txt+=feature+"\t"
			
			with open("annotation_data/"+protein+".tsv", "w") as g:
				g.write(txt+"\n")
	
		self.output().open("w").close()

	def requires(self):
		return DownloadRDFData_proteinTask(self.folder, self.prefix)

	def output(self):
		return luigi.LocalTarget(self.folder+"log/get_info_pairs.txt")

class DownloadRDFData_proteinTask(luigi.Task, TimeTaskMixin):
	folder = luigi.Parameter()
	prefix = luigi.Parameter()

	def run(self):
		print("Running Download rdf data")
		folder=self.folder
		prefix=self.prefix

		if( not os.path.isdir("rdf_data") ):
			os.system("mkdir rdf_data")
		if( not os.path.isdir("sequence_data") ):
			os.system("mkdir sequence_data")

		proteins=[]
		data = ['positive', 'false']
		for d in data:
			f=open(folder+prefix+d+".txt","r")
			for line in f:
				l=line.replace("\n","").split("\t")
				p1=l[0]
				p2=l[1]

				if( not(p1 in proteins) and not(os.path.isfile("rdf_data/"+p1+".rdf")) and not(os.path.isfile("sequence_data/"+p1+".fasta")) ):
					proteins.append(p1)
					
				if( not(p2 in proteins) and not(os.path.isfile("rdf_data/"+p2+".rdf")) and not(os.path.isfile("sequence_data/"+p2+".fasta")) ):
					proteins.append(p2)
					
			f.close()

		u=PreProcessing()
		for p in proteins:
			u.download_protein_info(folder, p)

		self.output().open("w").close()

	def output(self):
		return luigi.LocalTarget(self.folder+"log/download_rdf_data.txt")

