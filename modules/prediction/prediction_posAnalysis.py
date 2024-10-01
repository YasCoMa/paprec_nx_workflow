import os
import json
import joblib
import pandas as pd

import sys
workflow_path = os.environ.get('paprecPath')
if(workflow_path[-1]=='/'):
    workflow_path = workflow_path[:-1]

sys.path.insert(0, os.path.join( workflow_path, 'modules', "datasets") )
from treat_dataset import Dataset

from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, roc_auc_score, matthews_corrcoef, auc, precision_recall_curve
from sklearn.metrics import make_scorer

def pr_auc_score(ytrue, ypred, **kwargs):
    prec, rec, _ = precision_recall_curve(ytrue, ypred)
    prauc = auc(prec, rec)
    return prauc

class PredictionAnalysis:

	def __init__(self, rootFolder, modelsFolder, dataset_id, source ):
		self.ds = Dataset( rootFolder )

		self.root = rootFolder
		self.modelsPath = modelsFolder
		self.datasetFolder = os.path.join(rootFolder, dataset_id)
		self.resultsFolder = os.path.join( self.datasetFolder, 'results')
		if( not os.path.isdir(self.resultsFolder) ):
			os.mkdir(self.resultsFolder)

		with open( os.path.join(rootFolder, 'testset_combinations.json'), 'r' ) as g:
			auxid = dataset_id
			if( source == 'predictor'):
				auxid = '_'.join( dataset_id.split('_')[:-1] )
			self.pairs = json.load(g)[ auxid ]

	def _load_local_models_list(self):
		methods = []
		for d in os.listdir( self.datasetFolder ):
			if( os.path.isdir( os.path.join( self.datasetFolder, d ) ) and (d in self.pairs['methods'] ) ):
				methods.append(d)
		return methods

	def _load_test_dataset(self, method, subfolder, train_db, model_features):
		infile = os.path.join( self.datasetFolder, method, subfolder, 'computed_dataset.tsv' )
		df = pd.read_csv( infile, sep='\t' )
		
		test_cols = df.columns
		inter_features = []
		for f in model_features:
			if( f in test_cols ):
				inter_features.append(f)
		
		final_input = df[ inter_features ].fillna(0)
		self.x = final_input
		
		self.data_info = df[ ['item_id', 'item_sequence'] ]
		self.mapp_seq = dict( zip( df['item_sequence'].tolist(), df['item_id'].tolist() ) )

		"""
		metacols = ['item_id', 'item_sequence', 'label']
		feat_cols = list( filter( lambda x: (x not in metacols), df.columns ))
		selected_features = feat_cols
		if( subfolder == "all-features"):
			infile = os.path.join( self.root, train_db, method, subfolder, 'features_after.txt' )
			selected_features = open( infile ).read().split('\n')
		"""


	def _load_model(self, dataset, method):
		model = None
		subfolder = None
		target = f'{dataset},{method}'
		for f in os.listdir(self.modelsPath):
			if(f.startswith(target)):
				path = os.path.join( self.modelsPath, f )
				with open( path, 'r') as g:
					infom = json.load(g)

				infom['model_path'] = infom['model_path'].replace('paprec_dummy', 'paprec_data')
				subfolder = infom['model_path'].split('/')[-1].replace("-global_best","")
				if( os.path.isfile(infom['model_path']) ):
					model = joblib.load( infom['model_path'] )
				else:
					dirout = '/'.join(self.root.split('/')[:-1])
					if( self.root .endswith('/') ):
						dirout = '/'.join(self.root[:-1].split('/')[:-1])
						path = os.path.join( dirout, infom['model_path'] )
						model = joblib.load( path )
				#subfolder = infom['subfolder']
		return subfolder, model

	def _load_predictions(self):
		preds_by_sequence = {}
		preds_by_model = {}
		outfileS = os.path.join( self.datasetFolder, 'all_predictions_by_sequence.json')
		outfileM = os.path.join( self.datasetFolder, 'all_predictions_by_model.json')
		if( not os.path.isfile(outfileS) or not os.path.isfile(outfileM) ):
			for d in self.pairs['datasets']:
				for m in self.pairs['methods']:
					subfolder, model = self._load_model(d, m)
					model_features = model.feature_names_in_

					self._load_test_dataset( m, subfolder, d, model_features )
					
					labels = model.predict( self.x )
					labels = labels.tolist()
					try:
						probs = model.predict_proba( self.x )
						probs = probs.tolist()
					except:
						pass

					idx = 0
					_id = f"{d},{m}"
					preds_by_model[_id] = {}
					for i in self.data_info.index:
						seqId = str(self.data_info.loc[i, 'item_id'])
						seq = self.data_info.loc[i, 'item_sequence']
						
						if( not seqId in preds_by_sequence ):
							preds_by_sequence[seqId] = { 'sequence': seq, 'all': {}, 'highest_score': '' }
						
						preds_by_sequence[seqId]['all'][_id] = [ labels[idx], -1 ]
						if(len( probs) > 0):
							preds_by_sequence[seqId]['all'][_id] = [ labels[idx], probs[idx] ]

						aux = f"{seqId},{seq}"
						preds_by_model[_id][aux] = [ labels[idx], probs[idx] ]

						idx += 1
			
			for k in preds_by_sequence:
				scores = preds_by_sequence[k]['all']
				scores = dict( sorted( scores.items(), key=lambda item: item[1][0], reverse=True) )
				key = list(scores)[0]
				preds_by_sequence[k]['all'] = scores
				preds_by_sequence[k]['highest'] = { key: preds_by_sequence[k]['all'][key] }

			with open( outfileS, 'w') as g:
				json.dump(preds_by_sequence, g)
			with open( outfileM, 'w') as g:
				json.dump(preds_by_model, g)

		else:
			with open( outfileS, 'r') as g:
				preds_by_sequence = json.load( g)
			with open( outfileM, 'w') as g:
				preds_by_model = json.load(g)

		return preds_by_sequence, preds_by_model

	def perform_prediction(self ):
		predictionsS, predictionsM = self._load_predictions()

		selected_file = os.path.join( self.resultsFolder, 'selected_antigenic_items.tsv')
		g = open( selected_file, 'w' )
		g.write('item_id\titem_sequence\tratio_agreement\tmean_probability\n')
		
		outfile = os.path.join( self.resultsFolder, 'summary_table_predictions.tsv')
		f = open( outfile, 'w' )
		f.write('item_id\tdataset\tmethod\tlabel\tprobability\n')
		for item_id in predictionsS:
			item_sequence = predictionsS[item_id]['sequence']

			scores = predictionsS[item_id]['all']
			all_probs = []
			labels = []
			total = len(scores)
			for s in scores:
				dataset, method = s.split(',')
				label = scores[s][0]
				probability = scores[s][1][label]

				labels.append(label)
				if(label==1):
					all_probs.append(probability)

				f.write( f"{item_id}\t{dataset}\t{method}\t{label}\t{probability}\n")

			ratio = sum(labels)/total
			mean_prob = sum(all_probs)/total
			if( ratio > 0 ):
				g.write( f"{item_id}\t{item_sequence}\t{ratio}\t{mean_prob}\n")
		f.close()
		g.close()

	def _load_goldenset(self, goldenset):
		df = pd.read_csv(goldenset, sep='\t')
		info = dict( zip( df['sequence'].tolist(), df['label'].tolist() ))

		return info

	def perform_comparison(self, goldenset):
		comparison_config_file = os.path.join( self.datasetFolder, 'golden_config.json' )
		if( os.path.isfile(comparison_config_file) ):
			metrics = ['accuracy_score', 'f1_score', 'precision_score', 'recall_score', 'roc_auc_score', 'pr_auc_score','matthews_corrcoef']
			predictionsS, predictionsM = self._load_predictions()

			outfile = os.path.join( self.resultsFolder, 'models_performance.tsv')
			f = open( outfile, 'w' )
			header = [ 'source_comparison', 'dataset', 'method' ] + metrics
			header = '\t'.join(header)
			f.write( header+'\n' )

			with open(comparison_config_file) as g:
				config = json.load(g)

			for cnf in config:
				goldenset = cnf['path']
				source_id = cnf['identifier']

				gs = self._load_goldenset(goldenset)
				
				for mo in predictionsM:
					d, m = mo.split(',')

					target = self.ds.datasets[d]['target']
					if( target == cnf['target'] ):
						preds = []
						y = []
						for item in predictionsM[mo]:
							i, s = item.split(',')
							if( s in gs):
								y.append( gs[s] )
								preds.append( predictionsM[mo][item][0] )

						ms = {}
						for mt in metrics:
							ms[mt] = eval( f"{mt}(y, preds)" )

						values = '\t'.join( [ source_id, d, m] + [ str(v) for v in ms.values() ] )
						f.write( values+'\n' )
			
			f.close()