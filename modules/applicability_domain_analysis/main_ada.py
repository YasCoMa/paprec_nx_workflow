import os
import json
import shutil
import pandas as pd

import sys
workflow_path = os.environ.get('paprecPath')
if(workflow_path[-1]=='/'):
    workflow_path = workflow_path[:-1]

sys.path.insert(0, os.path.join( workflow_path, 'modules', "datasets") )
from treat_dataset import Dataset

sys.path.insert(0, os.path.join( workflow_path, 'modules', "applicability_domain_analysis") )
from perform_ada_analysis import ADAnalysis

import argparse

class ADAManager:

	def __init__(self, args):
		self.dataDir = args.data_directory

		self.ds = Dataset( self.dataDir )

		self.sim_metrics = ["cityblock", "euclidean", "rogerstanimoto", "correlation", "cosine"]

		with open( args.parameter_file, 'r' ) as g:
			self.config = json.load(g)

		tempFolder = os.path.join( self.dataDir, 'tmp' )
		if( not os.path.isdir(tempFolder) ):
			os.makedirs( tempFolder )

		tempFolderStep = os.path.join( tempFolder, 'ada' )
		self.tmpDir = tempFolderStep

		if( not os.path.isdir(tempFolderStep) ):
			os.mkdir( tempFolderStep )
			os.mkdir( os.path.join(tempFolderStep, 'tasks') )
			os.mkdir( os.path.join(tempFolderStep, 'ready') )

	def _mark_as_ready(self, prefix, task_id):
		outtask = open( os.path.join(self.tmpDir, 'ready', f"{prefix}_{task_id}.ready" ), 'w' )
		outtask.close()
		outtask = open( os.path.join(self.tmpDir, 'ada.ready' ), 'w' )
		outtask.close()

	def split_screening_tasks(self, mode):
		exps = self.config["experiment_combinations"]
		
		valid_queue = set()
		for e in exps:
			if( (e['task'] == 'ada') and (mode=='ada') ):
				target = e["target"]
				method = e["extraction_method"]
				proportion = e["perc_testset"]

				sim_metrics = self.sim_metrics
				if( 'similarity_metrics' in e ):
					mets = set( e['similarity_metrics'] ).intersection(sim_metrics)
					if( len(mets) > 0 ):
						sim_metrics = list(mets)
				sim_metrics = ';'.join(sim_metrics)

				all_all = False
				if( 'is_all_against_all' in e ):
					all_all = e['is_all_against_all']

				if( all_all ):
					exec_dss = self.config["datasets"]
					exec_dss = list( filter( lambda x: x in self.ds.datasets, exec_dss ))
					exec_dss = list( filter( lambda x: self.ds.datasets[x]["target"] == target, exec_dss ))
					for base in exec_dss:
						for test in exec_dss:
							_ide = f"{method},{base},{test},{proportion},{sim_metrics}"
							valid_queue.add(_ide)

				else:
					for p in e['dataset_pairs']:
						base = p['base']
						test = p['test']
						if( (base in self.ds.datasets) and (test in self.ds.datasets) ):
							if( (self.ds.datasets[base]["target"] == target) and (self.ds.datasets[test]["target"] == target) ):
								_ide = f"{method},{base},{test},{proportion},{sim_metrics}"
								valid_queue.add(_ide)

		taskFolder = os.path.join( self.tmpDir, 'tasks')
		if( os.path.isdir( taskFolder ) ):
			shutil.rmtree( taskFolder )
		os.mkdir( taskFolder )

		index = 1
		for v in valid_queue:
			taskFile = os.path.join(self.tmpDir, 'tasks', f"task_{index}.csv")
			with open( taskFile, 'w' ) as f:
				f.write( v + '\n' )
			index += 1

		stepout = os.path.join( self.tmpDir, 'ada.ready' )
		inn = os.listdir( os.path.join( self.tmpDir, 'tasks') )
		out = os.listdir( os.path.join( self.tmpDir, 'ready') )
		if(inn != out):
			if( os.path.isfile(stepout) ):
				os.remove(stepout)
		print('ok')

	def run(self, args):
		if( args.execution_mode == 1 ):
			self.split_screening_tasks(args.task_mode)

		execo = None
		outfolder = ''
		if( args.execution_mode >= 2 ):
			task_id = args.setup_instance.split('/')[-1].split('.')[0]
			method, base, test, proportion, sim_metrics = open( args.setup_instance ).read().split('\n')[0].split(',')
			proportion = int(proportion)
			sim_metrics = sim_metrics.split(';')

			execo = ADAnalysis( self.dataDir )
				
		if( args.execution_mode == 2 ):
			execo.run_ada( method, base, test, proportion, sim_metrics )
			self._mark_as_ready('ada', task_id)
		
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description='Applicability Domain Analysis module', formatter_class=RawTextHelpFormatter)

parser.add_argument("-em", "--execution_mode", action="store", help="1 - Run Task Splitting\n\
2 - Run Applicability domain analysis", type=int)

parser.add_argument("-dataDir", "--data_directory", type=str, default="", help="Folder to store the result files (ex.: /home/user/experiment/ )\n")
parser.add_argument("-paramFile", "--parameter_file", type=str, default="", help="Configuration file \n")
parser.add_argument("-mode", "--task_mode", type=str, default="", help="Goal of the analysis (train, test or ada) \n")
parser.add_argument("-setup", "--setup_instance", type=str, default="", help="Csv file with task execution instance ( mode, dataset, method and, if in test mode, the test set path) \n")

args = parser.parse_args()

o = ADAManager(args)
o.run(args)