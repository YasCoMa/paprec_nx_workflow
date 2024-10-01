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

sys.path.insert(0, os.path.join( workflow_path, 'modules', "prediction") )
from prediction_posAnalysis import PredictionAnalysis

import argparse

class Prediction:

	def __init__(self, args):
		self.dataDir = args.data_directory

		self.ds = Dataset( self.dataDir )

		with open( args.parameter_file, 'r' ) as g:
			self.config = json.load(g)

		tempFolder = os.path.join( self.dataDir, 'tmp' )
		if( not os.path.isdir(tempFolder) ):
			os.makedirs( tempFolder )

		tempFolderStep = os.path.join( tempFolder, 'prediction' )
		self.tmpDir = tempFolderStep

		if( not os.path.isdir(tempFolderStep) ):
			os.mkdir( tempFolderStep )
			os.mkdir( os.path.join(tempFolderStep, 'tasks') )
			os.mkdir( os.path.join(tempFolderStep, 'ready') )

	def _mark_as_ready(self, prefix, task_id):
		outtask = open( os.path.join(self.tmpDir, 'ready', f"{prefix}_{task_id}.ready" ), 'w' )
		outtask.close()
		outtask = open( os.path.join(self.tmpDir, 'prediction.ready' ), 'w' )
		outtask.close()

	def split_screening_tasks(self, mode):
		exps = self.config["experiment_combinations"]
		
		valid_queue = set()
		for e in exps:
			if( (e['task'] == 'test') and (mode=='test') ):
				if( 'test_data' in e ):
					for tds in e['test_data']:
						d = tds['identifier']
						
						configs = {}
						if('goldensets' in e and 'compare_to_goldenset' in e):
							if( e['compare_to_goldenset'] ):
								goldensets = e['goldensets']
								configs = { 'epitope': [], 'protein': [] }
								for gs in goldensets:
									configs[ gs['target'] ].append(gs)

						dss = [d]
						source = tds['source']
						if(source == 'predictor'):
							dss = []
							for t in ['epitope', 'protein']:
								selected_file = os.path.join( self.dataDir, d, 'data_parsing_selection', f'selected_{t}s.fasta')
								if( os.path.isfile(selected_file) ):
									dss.append( f"{d}_{t}" )

								if(len(configs) > 0):
									comparison_file = os.path.join( self.dataDir, f"{d}_{t}", 'golden_config.json')
									with open( comparison_file, 'w' ) as g:
										json.dump( configs[t] )
						else:
							if(len(configs) > 0):
								allc = configs['protein'] + configs['epitope']
								comparison_file = os.path.join( self.dataDir, d, 'golden_config.json')
								with open( comparison_file, 'w' ) as g:
									json.dump( allc )

						i=0
						for d in dss:
							_ide = f"{d},{source}"
							valid_queue.add(_ide)
							i+=1

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
		if( args.execution_mode >= 2 ):
			task_id = args.setup_instance.split('/')[-1].split('.')[0]
			dataset, source = open( args.setup_instance ).read().split('\n')[0].split(',')
			
			best_models_folder = os.path.join(self.dataDir, 'best_trained_models')

			execo = PredictionAnalysis( self.dataDir, best_models_folder, dataset, source)
				
		if( args.execution_mode == 2 ):
			execo.perform_prediction(  )
			self._mark_as_ready('prediction', task_id)
		
		if( args.execution_mode == 3 ):
			comparison_file = os.path.join( self.dataDir, d, 'golden_config.json')
			if( os.path.isfile( comparison_file ) ):
				global_model, model_id = execo.perform_comparison( goldenset_path, source )
			self._mark_as_ready('comparison', task_id)
		
from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description='Prediction module', formatter_class=RawTextHelpFormatter)

parser.add_argument("-em", "--execution_mode", action="store", help="1 - Run Task Splitting\n\
2 - Run prediction on new test dataset\n\
3 - Run performance comparison with scores from external tool or database", type=int)

parser.add_argument("-dataDir", "--data_directory", type=str, default="", help="Folder to store the result files (ex.: /home/user/experiment/ )\n")
parser.add_argument("-paramFile", "--parameter_file", type=str, default="", help="Configuration file \n")
parser.add_argument("-mode", "--task_mode", type=str, default="", help="Goal of the analysis (train or test) \n")
parser.add_argument("-setup", "--setup_instance", type=str, default="", help="Csv file with task execution instance ( mode, dataset, method and, if in test mode, the test set path) \n")

args = parser.parse_args()

o = Prediction(args)
o.run(args)