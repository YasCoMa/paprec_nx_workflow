import os
import random
import pandas as pd

class Dataset:
	
	def __init__(self, dataDir, dataset_id = None):
		workflow_path = os.environ.get('paprecPath')
		if(workflow_path[-1]=='/'):
		    workflow_path = workflow_path[:-1]

		self.get_available_datasets()

		self.raw_folder = os.path.join( workflow_path, 'raw_training_datasets' )
		
		if( dataset_id != None ):
			self.initialize_dataset(dataDir, dataset_id)

	def initialize_dataset(self, dataDir, dataset_id):
		self.is_training = False

		self.dataset_folder = os.path.join( dataDir, dataset_id )
		if( not os.path.isdir(self.dataset_folder) ):
			os.makedirs( self.dataset_folder )

		if( dataset_id in self.datasets ):
			self.is_training = True

			self.raw_dataset_folder = os.path.join( self.raw_folder, dataset_id )

			self.path_positive = f"{self.raw_dataset_folder}/{ self.datasets[dataset_id]['positive_file'] }"
			self.path_negative = f"{self.raw_dataset_folder}/{ self.datasets[dataset_id]['negative_file'] }"

		self.processed_dataset_folder = os.path.join( self.dataset_folder, 'processed_datasets' )
		if( not os.path.isdir(self.processed_dataset_folder) ):
			os.makedirs( self.processed_dataset_folder )

	def get_available_datasets(self):
		self.datasets = {
			"hla": { "target": "epitope", 'positive_file': 'dataset_pos.fasta', 'negative_file': 'dataset_neg.fasta' },
			"bcipep": { "target": "epitope", 'positive_file': 'dataset_pos.fasta', 'negative_file': 'dataset_neg.fasta' },
			"gram+_epitope": { "target": "epitope", 'positive_file': 'dataset_pos.fasta', 'negative_file': 'dataset_neg.fasta' },
			"gram-_epitope": { "target": "epitope", 'positive_file': 'dataset_pos.fasta', 'negative_file': 'dataset_neg.fasta' },
			"gram+_protein": { "target": "protein", 'positive_file': 'dataset_pos.fasta', 'negative_file': 'dataset_neg.fasta' },
			"gram-_protein": { 'positive_file': 'dataset_pos.fasta', 'negative_file': 'dataset_neg.fasta' },
			"allgram_epitope": { "target": "protein", "target": "epitope", 'positive_file': 'dataset_pos.fasta', 'negative_file': 'dataset_neg.fasta' },
			"allgram_protein": { "target": "protein", 'positive_file': 'dataset_pos.fasta', 'negative_file': 'dataset_neg.fasta' },
		}

	
	def load_sequences(self, fileSeq):
		dc = {}
		f = open( fileSeq, 'r' )
		for line in f:
			l = line.replace('\n','')
			if( l.startswith('>') ):
				_id = l[1:]
			else:
				dc[_id] = l
		f.close()

		return dc

	def get_processed_data(self, testSequenceFile = ""):
		outfile = os.path.join( self.processed_dataset_folder, 'processed_data.tsv' )

		if( not os.path.isfile(outfile) ):
			ids = []
			sequences = []
			classes = []

			if( (not self.is_training) and os.path.isfile(testSequenceFile) ):
				test = self.load_sequences( testSequenceFile )
				ids = list( test.keys() ) 
				sequences = list( test.values() ) 
				classes = ([-1] * len(ids))
			else:
				pos = self.load_sequences( self.path_positive )
				neg = self.load_sequences( self.path_negative )

				ids = list( pos.keys() ) + list( neg.keys() )
				sequences = list( pos.values() ) + list( neg.values() )
				classes = ([1] * len(pos)) + ([0] * len(neg))

			if( len(ids) > 0 ):
				df = pd.DataFrame()
				df['id'] = ids
				df['sequence'] = sequences
				df['label'] = classes
				df.to_csv( outfile, sep = '\t', index = None )
		else:
			df = pd.read_csv( outfile, sep='\t' )

		return df

	def get_random_balanced_processed_data(self):
		outfile = os.path.join( self.processed_dataset_folder, 'balanced_processed_data.tsv' )
		
		if( not os.path.isfile(outfile) ):
			df = self.get_processed_data()
			pos = df[ df['label'] == 1 ]
			neg = df[ df['label'] == 0 ]

			posidx = list(pos.index)
			negidx = list(neg.index)

			if( len(pos) < len(neg) ):
				negidx = random.sample( negidx, len(pos) )
			else:
				posidx = random.sample( posidx, len(neg) )

			indexes = posidx + negidx
			df = df[ df.index.isin(indexes) ]

			df.to_csv( outfile, sep = '\t', index = None )
		else:
			df = pd.read_csv( outfile, sep='\t' )

		return df

