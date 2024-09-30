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

sys.path.insert(0, os.path.join( workflow_path, 'modules', "feature_extraction_methods", "methods") )
from aln_free_e_descriptors import Implementation_vaxijen
from aln_free_aaindex_descriptors import Implementation_vaxijenModified
from aln_free_esm2_embedding import Implementation_es2Embedding

import argparse

class FeatureExtraction:

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

		tempFolderStep = os.path.join( tempFolder, 'tasks_feature_extraction' )
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

	def split_screening_tasks(self, mode):
		exps = self.config["experiment_combinations"]
		
		combinations_test_file = {}
		valid_queue = set()
		for e in exps:
			all_all = False
			if( 'is_all_against_all' in e ):
				all_all = e['is_all_against_all']

			imbalance_method = 'downsampling'
			if( e['task'] == 'train' and 'imbalance_method' in e ):
				if( e['imbalance_method'] in ['smote', 'downsampling'] ):
					imbalance_method = e['imbalance_method']

			exec_method = set()
			exec_dss = set()
			if( all_all ):
				exec_dss = self.config["datasets"]
				exec_dss = list( filter( lambda x: x in self.ds.datasets, exec_dss ))

				exec_methods = self.config["methods"]
				exec_methods = list( filter( lambda x: x in self.methods, exec_methods ))
				exec_methods = set(exec_methods)
				if( mode == 'train' ):
					for d in exec_dss:
						for m in exec_methods:
							_ide = f"train,{d},{m},,{imbalance_method}"
							valid_queue.add(_ide)

			else:
				for p in e['pairs']:
					d = p['dataset']
					m = p['method']
					if( (d in self.ds.datasets) and  (m in self.methods) ):
						exec_dss.add(d)
						exec_methods.add(m)
						if( mode == 'train' ):
							_ide = f"train,{d},{m},,{imbalance_method}"
							valid_queue.add(_ide)

			if( (e['task'] == 'test') and (mode=='test') ):
				if( 'test_data' in e ):
					for tds in e['test_data']:
						d = tds['identifier']
						combinations_test_file[d] = { 'datasets': list(exec_dss), 'methods': list(exec_methods) }
						
						dss = [d]
						paths = []
						source = tds['source']
						if(source == "fasta"):
							paths.append( tds['sequence_file'] )
						
						if(source == 'predictor'):
							dss = []
							paths = []
							for t in ['epitope', 'protein']:
								selected_file = os.path.join( self.dataDir, d, 'data_parsing_selection', f'selected_{t}s.fasta')
								if( os.path.isfile(selected_file) ):
									dss.append( f"{d}_{t}" )
									paths.append(selected_file)
						i=0
						for d in dss:
							for m in exec_methods:
								_ide = f"test,{d},{m},{paths[i]},"
								valid_queue.add(_ide)
							i+=1

		if( len(combinations_test_file) > 0 ):
			with open( os.path.join(self.dataDir, 'testset_combinations.json'), 'w' ) as g:
				json.dump(combinations_test_file, g)

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

		stepout = os.path.join( self.tmpDir, 'transformation.ready' )
		inn = os.listdir( os.path.join( self.tmpDir, 'tasks') )
		out = os.listdir( os.path.join( self.tmpDir, 'ready') )
		if(inn != out):
			if( os.path.isfile(stepout) ):
				os.remove(stepout)
		print('ok')

	def transform_datasets(self, task_id, mode, dataset, method, testPath, imbalance_method):
		self.ds.initialize_dataset( self.dataDir, dataset )

		outfolder = os.path.join( self.ds.dataset_folder, method )
		if( not os.path.isdir(outfolder) ):
			os.mkdir( outfolder )

			ds_input = None
			if( mode == 'train'):
				if(imbalance_method == 'smote'):
					ds_input = self.ds.get_processed_data()
				else:
					ds_input = self.ds.get_random_balanced_processed_data()

			if( mode == 'test'):
				ds_input = self.ds.get_processed_data(testPath)

			ds_input = ds_input.dropna()	
			executor_method = self.methods[method]
			extractor = eval( f"{ executor_method }()" )
			extractor.build_numerical_datasets( ds_input, outfolder )

		outtask = open( os.path.join(self.tmpDir, 'ready', f"{task_id}.ready" ), 'w' )
		outtask.close()
		outtask = open( os.path.join(self.tmpDir, 'transformation.ready' ), 'w' )
		outtask.close()

	def run(self, args):
		if( args.execution_mode == 1 ):
			self.split_screening_tasks(args.task_mode)

		if( args.execution_mode == 2 ):
			task_id = args.setup_instance.split('/')[-1].split('.')[0]
			mode, dataset, method, dataset_path, imbalance_method = open( args.setup_instance ).read().split('\n')[0].split(',')
			self.transform_datasets( task_id, mode, dataset, method, dataset_path, imbalance_method )
		
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description='Feature Extraction module', formatter_class=RawTextHelpFormatter)

parser.add_argument("-em", "--execution_mode", action="store", help="1 - Run Task Splitting\n\
2 - Run Dataset Transformation", type=int)

parser.add_argument("-dataDir", "--data_directory", type=str, default="", help="Folder to store the result files (ex.: /home/user/experiment/ )\n")
parser.add_argument("-paramFile", "--parameter_file", type=str, default="", help="Configuration file \n")
parser.add_argument("-mode", "--task_mode", type=str, default="", help="Goal of the analysis (train or test) \n")
parser.add_argument("-setup", "--setup_instance", type=str, default="", help="Csv file with task execution instance ( mode, dataset, method and, if in test mode, the test set path) \n")

if __name__ == '__main__':
	args = parser.parse_args()

	o = FeatureExtraction(args)
	o.run(args)