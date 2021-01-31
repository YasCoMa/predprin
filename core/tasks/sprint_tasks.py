import luigi
#from luigi.contrib.external_program import ExternalProgramTask

from ..time_counter import TimeTaskMixin

import urllib.request, os, subprocess

# prepare files
# initialize proteins with sequences (maybe not)
# get sequences -> download sequence
# concat proteins sequences
# predict interactions (external task)
# separate files
# generate file scores

# generate file scores -> separate files -> predict interactions (external task) -> concat sequences

class GenerateFileScoresTask(luigi.Task, TimeTaskMixin):
	folder = luigi.Parameter()
	prefix = luigi.Parameter()

	def run(self):
		print("Running Generate file scores")

		pairs={"pos": [], "neg": []}
		ds=['positive','false']
		for d in ds:
			cl='neg'
			if(d=='positive'):
				cl='pos'
			f=open(self.folder+"train_"+d+".txt","r")
			for line in f:
				l=line.replace("\n","").split(" ")
				pairs[cl].append(l[0]+","+l[1])
			f.close()

		f=open(self.folder+"scores.txt","w")
		f.close()
		
		ds=['pos','neg']
		for d in ds:
			c=0
			f=open(self.folder+"result.txt."+d)
			for line in f:
				ps=pairs[d][c].split(",")

				l=line.replace("\n","").split(" ")
				with open(self.folder+"scores.txt","a") as fg:
					fg.write(ps[0]+" "+ps[1]+" "+str(l[0])+"\n")
					fg.write(ps[1]+" "+ps[0]+" "+str(l[0])+"\n")
				c+=1
			f.close()

		self.output().open("w").close()

	def requires(self):
		return RunSprintPredictionTask(self.folder, self.prefix)
		
	def output(self):
		return luigi.LocalTarget(self.folder+"log/generate_file_score.txt")

class RunSprintPredictionTask(luigi.Task, TimeTaskMixin):
	folder = luigi.Parameter()
	prefix = luigi.Parameter()

	def run(self):
		arguments = [os.getcwd()+"/core/sprint/bin/predict_interactions", "-p "+self.folder+"new_protein.seq", "-h "+os.getcwd()+"/core/sprint/HSP/pre_computed_HSP", "-tr "+self.folder+"train_positive.txt", "-pos "+self.folder+"train_positive.txt", "-neg "+self.folder+"train_false.txt", "-o "+self.output().path ]
		subprocess.call(" ".join(arguments), shell=True)
		#self.output().open("w").close()
	
	def requires(self):
		return ConcatSequencesDsToDatabaseTask(self.folder, self.prefix)

	def output(self):
		return luigi.LocalTarget(self.folder+"result.txt")
"""
class RunSprintPredictionTask(luigi.contrib.external_program.ExternalProgramTask, TimeTaskMixin):
	folder = luigi.Parameter()
	prefix = luigi.Parameter()

	def run(self):
		return ["./core/sprint/bin/predict_interactions", "-p "+self.folder+"new_protein.seq", "-h core/sprint/HSP/pre_computed_HSP", "-tr "+self.folder+"train_positive.txt", "-e", "-o "+self.output().path ]
	
	def requires(self):
		return ConcatSequencesDsToDatabaseTask(self.folder, self.prefix)

	def output(self):
		return luigi.LocalTarget(self.folder+"result_interactome.txt")
"""

class ConcatSequencesDsToDatabaseTask(luigi.Task, TimeTaskMixin):
	folder = luigi.Parameter()
	prefix = luigi.Parameter()

	def run(self): # Remove " and protein entries without sequence occurrences from all_proteins.seq before compute HSPs
		f=open(self.folder+"new_protein.seq","w")
		f.close()

		proteins=[]
		f=open("core/sprint/uniprot_seq_all.fasta","r")
		for line in f:
			l=line.replace("\n","")
			if(l.find(">")!=-1):
				proteins.append(l)

			with open(self.folder+"new_protein.seq","a") as fg:
				fg.write(l.replace('"',"")+"\n")
		f.close()

		f=open(self.folder+"protein.seq","r")
		for line in f:
				l=line.replace("\n","")        
				if(l.find(">")!=-1):
					if(not (l in proteins) ):
						with open(self.folder+"new_protein.seq","a") as fg:
							fg.write(l+"\n")
						flag=True
					else:
						flag=False

				elif(l!="" and flag):
					with open(self.folder+"new_protein.seq","a") as fg:
						fg.write(l.replace('"',"")+"\n")
		f.close()

		self.output().open("w").close()

	def requires(self):
		return GetSequencesTask(self.folder, self.prefix)

	def output(self):
		return luigi.LocalTarget(self.folder+"log/concat_sequences.txt")

class GetSequencesTask(luigi.Task, TimeTaskMixin):
	folder = luigi.Parameter()
	prefix = luigi.Parameter()

	def run(self):
		print("Getting sequences...") # use this way of writing tasks when you have to get all the results before return to the task which is waiting

		f=open(self.folder+"protein.seq", "w")
		f.close()
		
		proteins=[]
		ds=['positive', 'false']
		for d in ds:
			f=open(self.folder+self.prefix+d+".txt")
			for line in f:
				l=line.replace("\n","").split("\t")
				p1=l[0]
				p2=l[1]
				if( not(p1 in proteins) ):
					proteins.append(p1)

				if( not(p2 in proteins) ):	
					proteins.append(p2)

			f.close()

		for p in proteins:
			protein = p
			link=protein+".fasta not found"
			try:
				f=open("sequence_data/"+protein+".fasta","r")
				for line in f:
					with open(self.folder+"protein.seq", "a") as myfile:
						myfile.write(line)
				f.close()
			except:
				print("Not working: "+link)

		self.output().open("w").close()

	def requires(self):
		return PrepareFilesTask(self.folder, self.prefix)

	def output(self):
		return luigi.LocalTarget(self.folder+"log/get_sequences.txt")

"""
class DownloadProteinSequenceTask(luigi.Task, TimeTaskMixin):
	folder = luigi.Parameter()
	protein = luigi.Parameter()

	def run(self):
		link="http://www.uniprot.org/uniprot/"+self.protein+".fasta"
		try:
			f = urllib.request.urlopen(link)
			file = f.read()
			f=open(self.folder+"tempseq.txt","w")
			f.writelines(str(file))
			f.close()

			f=open(self.folder+"tempseq.txt","r")
			for line in f:
				with open(self.folder+"protein.seq", "a") as myfile:
					myfile.write(">"+self.protein+"\n")
				l=str(line).replace("b'","").replace("'","").split("\\n")
				c=0
				for l_ in l:
					if(l_!="" and l_.find(">")==-1):
						with open(self.folder+"protein.seq", "a") as myfile:
							myfile.write(l_)
					c+=1
					if(c==len(l)):
						with open(self.folder+"protein.seq", "a") as myfile:
							myfile.write("\n")
			f.close()

			os.remove(self.folder+"tempseq.txt")
		except:
			print("Not working: "+link)

		self.output().open("w").close()

	def output(self):
		return luigi.LocalTarget(self.folder+"log/download_save_sequence.txt")
"""

class PrepareFilesTask(luigi.Task, TimeTaskMixin):
	folder = luigi.Parameter()
	prefix = luigi.Parameter()

	def run(self):
		ds=['positive','false']
		for d in ds:
			f=open(self.folder+"train_"+d+".txt","w")
			f.close()

			f=open(self.folder+self.prefix+d+".txt","r")
			for line in f:
				l=line.replace("\n","").split("\t")
				p1=l[0]
				p2=l[1]

				with open(self.folder+"train_"+d+".txt","a") as fg:
					fg.write(p1+" "+p2+"\n")
			f.close()

		self.output().open("w").close()

	def output(self):
		return luigi.LocalTarget(self.folder+"log/prepare_files.txt")