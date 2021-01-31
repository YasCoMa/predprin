from core.tasks.protein_annotation_tasks import ProteinAnnotationTask
from core.tasks.generate_features_tasks import GenerateFeaturesTask
from core.tasks.classify_evaluate_tasks import ClassifyEvaluateTask
from core.time_counter import TimeTaskMixin

from core.utils import Communication
from core.utils import PreProcessing

import luigi
import os

class RunPPIExperiment(luigi.Task, TimeTaskMixin):
	parameters_file = luigi.Parameter()
	mode = luigi.Parameter()
	model = luigi.Parameter()

	def run(self): 
		print("Running experiment on datasets")
		
		self.output().open("w").close()

	def requires(self):
		p=PreProcessing()
		valid, data, error_message =  p.validate_parameter_file(self.parameters_file, self.mode, self.model)
		if(valid):
			for d in data["datasets"]:
				p.cleaning(d["folder"])
		
			print("Running Experiment: "+data["name"])
			os.environ['email_owner']=data['email']
			
			f=open(d["folder"]+"init.txt","w")
			f.close()

			for d in data["datasets"]:
				yield PredRep_ppi_Workflow(d["folder"], d["prefix"], self.mode, self.model)
		else:
			print("ERROR >>> "+error_message)

	def output(self):
		return luigi.LocalTarget("run_experiment.txt")

class PredRep_ppi_Workflow(luigi.Task, TimeTaskMixin):
	folder = luigi.Parameter()
	prefix = luigi.Parameter()
	mode = luigi.Parameter()
	model = luigi.Parameter()

	def run(self): # identifier for the positive and false part
		# only input: two files of interactions, one positive and other false, with the pair and class, separated by tab
		print("Running PredRep workflow")
		
		self.output().open("w").close()

	def requires(self):
		return PredRep_ppi_step3(self.folder, self.prefix, self.mode, self.model)

	def output(self):
		return luigi.LocalTarget(self.folder+"log/ppi_workflow.txt")


class PredRep_ppi_step3(luigi.Task, TimeTaskMixin):
	folder = luigi.Parameter()
	prefix = luigi.Parameter()
	mode = luigi.Parameter()
	model = luigi.Parameter()

	def run(self):
		print("Running classification and evaluation")

		yield ClassifyEvaluateTask(self.folder, self.prefix, self.mode, self.model)
		
		self.output().open("w").close()

	def requires(self):
		return PredRep_ppi_step2(self.folder, self.prefix)

	def output(self):
		return luigi.LocalTarget(self.folder+"log/pred_rep_step3.txt")

class PredRep_ppi_step2(luigi.Task, TimeTaskMixin):
	folder = luigi.Parameter()
	prefix = luigi.Parameter()

	def run(self):
		print("Running Generating features")
		
		yield GenerateFeaturesTask(self.folder, self.prefix)
		
		self.output().open("w").close()

	def requires(self):
		return PredRep_ppi_step1(self.folder, self.prefix)

	def output(self):
		return luigi.LocalTarget(self.folder+"log/pred_rep_step2.txt")

class PredRep_ppi_step1(luigi.Task, TimeTaskMixin):
	folder = luigi.Parameter()
	prefix = luigi.Parameter()

	def run(self):
		print("Running Protein annotation")
		
		yield ProteinAnnotationTask(self.folder, self.prefix)
		
		self.output().open("w").close()

	def output(self):
		return luigi.LocalTarget(self.folder+"log/pred_rep_step1.txt")

# Semantic description
# https://link.springer.com/article/10.1186/s13321-016-0168-9
# https://www.sciencedirect.com/science/article/pii/S1570826811000229?via%3Dihub
# https://www.w3.org/TR/2013/REC-prov-o-20130430/

# command to run luigi: luigid (using the file client.cfg)
# PYTHONPATH='.' luigi --module main RunPPIExperiment --parameters-file params_server.json --workers 3 --worker-keep-alive
#  python3 -m luigi --module main RunPPIExperiment --parameters-file params_server.json --workers 3

import sys, os
#if __name__ == '__main__':
#	os.system("rm run_experiment.txt")
#	p=PreProcessing()
#	for i in range(1,4):
#		#p.cleaning("/home/yasmmin/Dropbox/lncc/tese/predrep_as_workflow/aux_ds_"+str(i)+"/")
#		p.cleaning("/home/yasminn/Documentos/tese/predrep_as_workflow/data_test_workflow/aux_ds_"+str(i)+"/")

PredRep_ppi_Workflow.event_handler(luigi.Event.SUCCESS)(Communication().send_success_email)
#	luigi.build( [ RunPPIExperiment(sys.argv[1]) ] ) # parameters file
