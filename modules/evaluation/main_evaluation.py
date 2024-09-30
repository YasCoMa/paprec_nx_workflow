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

sys.path.insert(0, os.path.join( workflow_path, 'modules', "evaluation") )
from pipeline_model_evaluation import PipelineRankModel

import argparse

class Evaluation:

	def __init__(self, args):
		self.dataDir = args.data_directory
		
		self.ds = Dataset( self.dataDir )
		self.training_datasets = list(self.ds.datasets)
		self.load_available_methods()

		with open( args.parameter_file, 'r' ) as g:
			self.config = json.load(g)

		tempFolder = os.path.join( self.dataDir, 'tmp' )
		if( not os.path.isdir(tempFolder) ):
			os.makedirs( tempFolder )

		tempFolderStep = os.path.join( tempFolder, 'evaluation' )
		self.tmpDir = tempFolderStep

		if( not os.path.isdir(tempFolderStep) ):
			os.mkdir( tempFolderStep )
			os.mkdir( os.path.join(tempFolderStep, 'tasks') )
			os.mkdir( os.path.join(tempFolderStep, 'ready') )
		
	def load_available_methods(self):
		# Identifier | ImplementationClass Name
		self.methods = {
			"aln_free_e_descriptors": "Implementation_vaxijen",
			"aln_free_aaindex_descriptors": "Implementation_vaxijenModified",
			"aln_free_esm2_embedding": "Implementation_es2Embedding"
		}

	def _mark_as_ready(self, prefix, task_id):
		outtask = open( os.path.join(self.tmpDir, 'ready', f"{prefix}_{task_id}.ready" ), 'w' )
		outtask.close()
		outtask = open( os.path.join(self.tmpDir, 'evaluation.ready' ), 'w' )
		outtask.close()

	def split_screening_tasks(self, mode):
		exps = self.config["experiment_combinations"]
		
		valid_queue = set()
		for e in exps:
			if( (e['task'] == 'train') and (mode=='train') ):
				all_all = False
				if( 'is_all_against_all' in e ):
					all_all = e['is_all_against_all']

				imbalance_method = 'downsampling'
				if( 'imbalance_method' in e ):
					if( e['imbalance_method'] in ['smote', 'downsampling'] ):
						imbalance_method = e['imbalance_method']

				rank_metric = 'mcc'
				if( 'rank_metric' in e ):
					if( e['rank_metric'] in ['accuracy', 'auc', 'recall', 'precision', 'f1', 'kappa', 'mcc'] ):
						rank_metric = e['rank_metric']

				exec_method = set()
				if( all_all ):
					exec_dss = self.config["datasets"]
					exec_dss = list( filter( lambda x: x in self.ds.datasets, exec_dss ))

					exec_methods = self.config["methods"]
					exec_methods = list( filter( lambda x: x in self.methods, exec_methods ))
					exec_methods = set(exec_methods)
					if( mode == 'train' ):
						for d in exec_dss:
							for m in exec_methods:
								_ide = f"{d},{m},,{imbalance_method},{rank_metric}"
								valid_queue.add(_ide)

				else:
					for p in e['pairs']:
						d = p['dataset']
						m = p['method']
						if( (d in self.ds.datasets) and  (m in self.methods) ):
							exec_methods.add(m)
							_ide = f"{d},{m},,{imbalance_method},{rank_metric}"
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

		stepout = os.path.join( self.tmpDir, 'evaluation.ready' )
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
		imbalance_method = 'downsampling'
		rank_metric = 'mcc'
		if( args.execution_mode >= 2 ):
			task_id = args.setup_instance.split('/')[-1].split('.')[0]
			dataset, method, dataset_path, imbalance_method, rank_metric = open( args.setup_instance ).read().split('\n')[0].split(',')
			self.ds.initialize_dataset( self.dataDir, dataset )
			outfolder = os.path.join( self.ds.dataset_folder, method )
			if( not os.path.isdir(outfolder) ):
				sys.exit( f"Dataset/method - {dataset}/{method} combination was not processed yet")
			
			best_models_folder = os.path.join(self.dataDir, 'best_trained_models')
			if( not os.path.isdir(best_models_folder) ):
				os.mkdir(best_models_folder)

			execo = PipelineRankModel()
				
		if( args.execution_mode == 2 ):
			execo.make_feature_selection_general_dataset( outfolder, imbalance_method)
			self._mark_as_ready('feature-selection', task_id)
		
		if( args.execution_mode == 3 ):
			global_model, model_id = execo.build_models_and_ranking( outfolder, imbalance_method, rank_metric)
			
			prefix = model_id.replace('-global_best', '')
			info = { "dataset": dataset, "method": method, "subfolder": prefix, "model_path": global_model }
			model_export = os.path.join( best_models_folder, f'{dataset},{method},{prefix}.json' )
			with open( model_export, 'w' ) as g:
				json.dump(info, g)

			self._mark_as_ready('ml-pipeline-ranking', task_id)
		
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description='Evaluation module', formatter_class=RawTextHelpFormatter)

parser.add_argument("-em", "--execution_mode", action="store", help="1 - Run Task Splitting\n\
2 - Run feature selection on general dataset\n\
3 - Run ML pipeline to build models and rank them", type=int)

parser.add_argument("-dataDir", "--data_directory", type=str, default="", help="Folder to store the result files (ex.: /home/user/experiment/ )\n")
parser.add_argument("-paramFile", "--parameter_file", type=str, default="", help="Configuration file \n")
parser.add_argument("-mode", "--task_mode", type=str, default="", help="Goal of the analysis (train or test) \n")
parser.add_argument("-setup", "--setup_instance", type=str, default="", help="Csv file with task execution instance ( mode, dataset, method and, if in test mode, the test set path) \n")

args = parser.parse_args()

o = Evaluation(args)
o.run(args)