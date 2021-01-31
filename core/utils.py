import luigi
import datetime
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import os
import json
import urllib.request

class Utils_search:

	def get_annotated_proteins(self, folder):
		list_annotated=[]
		if(not os.path.isdir("annotation_data")):
			os.system("mkdir annotation_data")
			
		for f in os.listdir("annotation_data"):
			list_annotated.append(f.replace(".tsv",""))
				
		return list_annotated

class Communication:

	def send_success_email(self, task: luigi.Task) -> None:
	    sender_email=""
	    sender_server=""
	    sender_port=""
	    sender_password=""
	    
	    if(sender_email!=""):
		    fromaddr = sender_email
		    toaddr = os.environ['email_owner']
		    msg = MIMEMultipart()
		    msg['From'] = fromaddr
		    msg['To'] = toaddr
		    name_task_class = task.task_family
		    msg['Subject'] = name_task_class+" task ran successfully on PredRep workflow"
		    body = "<p>\n"
		    body += " ".join((
			    name_task_class,
			    "task finished successfully at",
			    datetime.datetime.now().strftime("%H:%M:%S on %Y-%m-%d"),
			    ))
		    body += "\n</p>\n"
		    # Then start mining task attributes, pick some off, and dump them to a table
		    body += "<table align=\"left\" border=1>\n"
		    body += "\n".join([
			    "<tr><td>{}</td><td>{}</td></tr>".format(key, value)
			    for key, value in (
				    ("Task ID", task.task_id),
				    ("Task Class Name", name_task_class),
				    ("Parameter(s)", task.get_params()),
				    ("Requirement(s)", task.requires()),
				    ("Input(s)", task.input()),
				    ("Output(s)", task.output()),
				    )
			    ])
		    body += "\n</table>\n"
		    msg.attach(MIMEText(body, 'html'))

		    try:
			    server = smtplib.SMTP(sender_server, sender_port)
			    server.starttls()
			    server.login(fromaddr, sender_password)
			    text = msg.as_string()
			    server.sendmail(fromaddr, toaddr, text)
			    server.quit()

			    valid=True
		    except:
			    valid=False

		return None

class PreProcessing:
	def cleaning(self, folder):
		print("Cleaning")
		# Cleaning files to execute whole workflow
		#if(os.path.isdir(folder+"log")):
		#	os.system("rm -rf "+folder+"log")
		#os.system("rm "+folder+"log/initilizing_variables.txt")
		#os.system("rm "+folder+"log/generate_features.txt")
		#os.system("rm "+folder+"log/get_feature_vector.txt")
		#os.system("rm "+folder+"log/pred_rep_step2.txt")

		#os.system("rm "+folder+"log/classification_evaluation.txt")
		#os.system("rm "+folder+"log/plot_results.txt")
		#os.system("rm "+folder+"log/classification_evaluation.txt")
		#os.system("rm "+folder+"log/do_classification_log.txt")
		#os.system("rm "+folder+"log/do_cross_validation.txt")
		#os.system("rm "+folder+"log/process_data.txt")
		#os.system("rm "+folder+"log/pred_rep_step3.txt")

		#os.system("rm "+folder+"log/ppi_workflow.txt")

		#os.system("rm run_experiment.txt")

		# Step 1
		#os.system("rm rdf_data/*.rdf")
		#os.system("rm "+folder+"new_info_proteins_*")
		
		# Step 2
		#os.system("rm "+folder+"train_*")
		#os.system("rm "+folder+"protein.seq")
		#os.system("rm "+folder+"new_protein.seq")
		#os.system("rm "+folder+"result.txt*")
		#os.system("rm "+folder+"dataset_ppi.txt")

		# Step 3
		#os.system("rm "+folder+"*.png")
		#os.system("rm "+folder+"*.npy")
		#os.system("rm "+folder+"evaluation_log-0.8.txt")
		#os.system("rm "+folder+"result_cross-validation.txt")

	def validate_parameter_file(self, file, mode, model):
		valid=False
		data=None
		error_message=""
		try:
			if(not mode in ['train', 'test']):
				error_message+="\nThe mode parameter should be either train or test"
		
			if(mode=='test'):
				if(model==None or model==''):
					error_message+="\nThe model file was not given as input"
				else:
					if(not os.path.isfile(model)):
						error_message+="\nThe model file was not found"
					else:
						if(not model.endswith(".joblib")):
							error_message+="\nThe model has not the correct extension (joblib)"
			
			with open(file) as json_file:
				data = json.load(json_file)

			if(not("name" in data.keys()) or not("description" in data.keys()) or not("owner" in data.keys()) or not("email" in data.keys()) or not("datasets" in data.keys()) or not(type(data["datasets"])==list) ):
				error_message+="\nThe input does not have all or some of the following fields: datasets, name, description, owner and email"
			elif (data["name"]=="" or data["owner"]=="" or data["email"]=="" or data["name"]==None or data["owner"]==None or data["email"]==None):
				error_message+="\nYou did not fill all the fields: name, description, owner and email\n"
			elif( not type(data["name"])==str or not type(data["owner"])==str or not type(data["email"])==str ):
				error_message+= "\nThe following fields must be string: name, description, owner and email"
			elif(len(data["datasets"])==0):
				error_message+="\nThere is no datasets in the list\n"
			else:
				for d in data["datasets"]:
					if ( not os.path.isdir(d["folder"]) ):
						error_message+="\nThe dataset folder"+d["folder"]+" does not exist\n"
					else:
						types=["false","positive"]
						for t in types:
							keyfile=d["folder"]+d["prefix"]+t+".txt"
							if ( not os.path.isfile(keyfile) ):
								error_message+="\nThe dataset file "+d["folder"]+" with "+t+" pairs with uniprot identifiers separated by tab does not exist (using .txt extension) \n"
			
			if(error_message==""):
				valid=True
			else:
				data=None
		except:
			valid=False
			data=None
			error_message="Problem on loading json file"

		return valid, data, error_message

	def download_protein_info(self, folder, p):
		
		try:
			link = "https://www.uniprot.org/uniprot/"+p+".rdf"
			f = urllib.request.urlopen(link)
			file = f.read()
			f=open("rdf_data/"+p+".rdf","w")
			f.writelines(str(file).replace("b'","").replace("'",'"').replace('\\"','"').replace("\\n",'\n').replace('>"','>'))
			f.close()
			
			new_link=""
			id_=p
			f=open("rdf_data/"+p+".rdf","r")
			for line in f:
			    l=line.replace("\n","")
			    if(l.find("replacedBy")!=-1):
			        new_link=l.split("=")[1].replace('"',"").replace("/>","")
			        break
			f.close()
			
			if(new_link!=""):
			    id_=new_link.split("/")[-1]
			    f = urllib.request.urlopen(new_link+".rdf")
			    file = f.read()
			    f=open("rdf_data/"+p+".rdf","w")
			    f.writelines(str(file).replace("b'","").replace("'",'"').replace('\\"','"').replace("\\n",'\n').replace('>"','>'))
			    f.close()

			link = "https://www.uniprot.org/uniprot/"+id_+".fasta"
			f = urllib.request.urlopen(link)
			file = f.read()
			f=open(folder+"tempseq.txt","w")
			f.writelines(str(file))
			f.close()

			f=open("sequence_data/"+p+".fasta", "w")
			f.close()

			f=open(folder+"tempseq.txt","r")
			for line in f:
				with open("sequence_data/"+p+".fasta", "a") as myfile:
					myfile.write(">"+p+"\n")
				l=str(line).replace("b'","").replace("'","").split("\\n")
				c=0
				for l_ in l:
					if(l_!="" and l_.find(">")==-1):
						with open("sequence_data/"+p+".fasta", "a") as myfile:
							myfile.write(l_)
					c+=1
					if(c==len(l)):
						with open("sequence_data/"+p+".fasta", "a") as myfile:
							myfile.write("\n")
			f.close()
		except:
			pass
