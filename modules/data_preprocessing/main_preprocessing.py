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

sys.path.insert(0, os.path.join( workflow_path, 'modules', "data_preprocessing") )
from parser_selection import ParserCuration

import argparse

class DataPreprocessing:

	def __init__(self, args):
		self.dataDir = args.data_directory
		
		self.ds = Dataset( self.dataDir )

		with open( args.parameter_file, 'r' ) as g:
			self.config = json.load(g)

		tempFolder = os.path.join( self.dataDir, 'tmp' )
		if( not os.path.isdir(tempFolder) ):
			os.makedirs( tempFolder )

		tempFolderStep = os.path.join( tempFolder, 'tasks_data_preprocessing' )
		self.tmpDir = tempFolderStep

		if( not os.path.isdir(tempFolderStep) ):
			os.mkdir( tempFolderStep )
			os.mkdir( os.path.join(tempFolderStep, 'tasks') )
			os.mkdir( os.path.join(tempFolderStep, 'ready') )

	def _mark_as_ready(self, prefix, task_id):
		outtask = open( os.path.join(self.tmpDir, 'ready', f"{prefix}_{task_id}.ready" ), 'w' )
		outtask.close()
		outtask = open( os.path.join(self.tmpDir, 'preprocessing.ready' ), 'w' )
		outtask.close()

	def split_processing_tasks(self, mode):
		exps = self.config["experiment_combinations"]
		
		valid_queue = set()
		for e in exps:

			if( (e['task'] == 'test') and (mode=='test') ):
				if( 'test_data' in e ):
					for tds in e['test_data']:
						if( tds['source'] == 'predictor' ):
							d = tds['identifier']
							protein_file = tds['proteins_file']
							raw_file = tds['raw_prediction_file']
							
							default = ['t', 1, 30, 2]
							cparams = ['cell_type', 'threshold_sim_iedb', 'threshold_alleles', 'threshold_rank']
							if( 'parameters_curation' in tds ):
								i=0
								for p in cparams:
									if p in tds['parameters_curation']:
										default[i] = tds['parameters_curation'][p]
									i+=1
							cell_type, threshold_sim_iedb, threshold_alleles, threshold_rank = tuple( default )
							cell_type = cell_type.lower()

							_ide = f"{d},{protein_file},{raw_file},{cell_type},{threshold_sim_iedb},{threshold_alleles},{threshold_rank}"
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

		stepout = os.path.join( self.tmpDir, 'preprocessing.ready' )
		inn = os.listdir( os.path.join( self.tmpDir, 'tasks') )
		out = os.listdir( os.path.join( self.tmpDir, 'ready') )
		if(inn != out):
			if( os.path.isfile(stepout) ):
				os.remove(stepout)
		print('ok')

	def run(self, args):
		if( args.execution_mode == 1 ):
			self.split_processing_tasks(args.task_mode)

		parser = None
		df = None
		cell_type = None
		if( args.execution_mode >= 2 ):
			task_id = args.setup_instance.split('/')[-1].split('.')[0]
			dataset_id, protein_file, raw_file, cell_type, threshold_sim_iedb, threshold_alleles, threshold_rank = open( args.setup_instance ).read().split('\n')[0].split(',')
			parser = ParserCuration( self.dataDir, dataset_id, protein_file)

		if( args.execution_mode == 2 ):
			parser.parse_prediction_results( raw_file )
			self._mark_as_ready('parsing', task_id)

		if( cell_type == 't'):
			if( args.execution_mode == 3 ):
				parser.filter_rank( float(threshold_rank) )
				self._mark_as_ready('rank', task_id)

			if( args.execution_mode == 4 ):
				parser.check_allele_promiscuity( float(threshold_alleles) )
				self._mark_as_ready('promiscuity', task_id)

		if( args.execution_mode == 5 ):
			parser.check_overlapping_epis_violinet( )
			self._mark_as_ready('filter-violinet', task_id)
		
		if( args.execution_mode == 6 ):
			parser.filter_iedb_epitopes( float(threshold_sim_iedb) )
			self._mark_as_ready('filter-iedb', task_id)
		
		if( args.execution_mode == 7 ):
			parser.check_human_homology( )
			self._mark_as_ready('human-proteome', task_id)

		if( args.execution_mode == 8 ):
			parser.prepare_final_sequences( )
			self._mark_as_ready('final-files', task_id)

from argparse import RawTextHelpFormatter

parser = argparse.ArgumentParser(description='Feature Extraction module', formatter_class=RawTextHelpFormatter)

parser.add_argument("-em", "--execution_mode", action="store", help="1 - Run Task Splitting\n\
2 - Parse results from netMHCPan tool\n\
3 - Filter by rank percentile BA\n\
4 - Filter by allele promiscuity\n\
5 - Filter by overlapping in violinet db\n\
6 - Filter by overlapping in iedb db\n\
7 - Filter by overlapping in human proteins\n\
8 - Write sequence selection and reports", type=int)

parser.add_argument("-dataDir", "--data_directory", type=str, default="", help="Folder to store the result files (ex.: /home/user/experiment/ )\n")
parser.add_argument("-paramFile", "--parameter_file", type=str, default="", help="Configuration file \n")
parser.add_argument("-mode", "--task_mode", type=str, default="", help="Goal of the analysis (train or test) \n")
parser.add_argument("-setup", "--setup_instance", type=str, default="", help="Csv file with task execution instance ( mode, dataset, method and, if in test mode, the test set path) \n")

args = parser.parse_args()

o = DataPreprocessing(args)
o.run(args)