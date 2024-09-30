import os
import joblib
import shutil
import pandas as pd
from pycaret.classification import *
# solved: pip3 install pycaret==3.3.1

class PipelineRankModel:

	def __init__(self):
		self.classifier_universe = { 'Extra Trees Classifier': 'et', 'Random Forest Classifier': 'rf', 'Ada Boost Classifier': 'ada', 'Gradient Boosting Classifier': 'gbc', 'Decision Tree Classifier': 'dt', 'SVM - Linear Kernel': 'svm', 'Naive Bayes': 'nb', 'CatBoost Classifier': 'catboost', 'Light Gradient Boosting Machine': 'lightgbm', 'Extreme Gradient Boosting': 'xgboost' }
		self.variation_classifiers = ['et', 'rf', 'ada', 'gbc', 'dt', 'svm', 'nb']
		self.xai_classifiers = ['dt', 'catboost', 'lightgbm', 'et', 'rf', 'xgboost']

		self.metrics = ['acc', 'auc', 'recall', 'precision', 'f1', 'kappa', 'mcc']
		self.metric_titles = {'acc': 'Accuracy', 'auc': 'AUC', 'recall': 'Recall', 'precision': 'Prec.', 'f1': 'F1', 'kappa': 'Kappa', 'mcc': 'MCC'}
	
	def _initialize_summary_table(self, outfolder):
		summary_file = os.path.join( outfolder, 'summary_best_models.tsv' )
		
		header = ['type_analysis', 'model' ] + self.metrics
		f = open( summary_file, 'w' )
		f.write( ('\t'.join(header))+'\n' )
		f.close()

		return summary_file

	def _initialize_experiment(self, wd, imbalanceMethod = 'downsampling', has_feat_selection = False ):
		fix_classes = False 
		if( imbalanceMethod == 'smote'):
			fix_classes = True

		exp = ClassificationExperiment()
		
		infile = os.path.join( wd, "computed_dataset.tsv")
		data = pd.read_csv( infile, sep='\t' )
		
		metacols = ['item_id', 'item_sequence']
		feat_cols = list( filter( lambda x: (x not in metacols), data.columns ))
		data = data[ feat_cols ]

		exp.setup(data = data, target = 'label', fix_imbalance = fix_classes, feature_selection = has_feat_selection, session_id = 42)
		
		if( has_feat_selection ):
			features_before = data.columns
			features_after = exp.dataset_transformed.columns
			for state in ['before', 'after']:
				outfile = os.path.join( wd,  f"features_{state}.txt"  )
				f = open(outfile, 'w')
				f.write( '\n'.join( list( eval( f"features_{state}") ) )  )
				f.close()

		return exp

	def _export_subfolder_model(self, wd, exp, include_explainability):
		model_file = os.path.join( wd, 'best_scored')
		best = None
		if( include_explainability):
			best = exp.compare_models( include = self.xai_classifiers )
		else:
			best = exp.compare_models( include = self.variation_classifiers )
		joblib.dump(best, model_file)

		return best, model_file

	def _export_subfolder_scores(self, wd, subfolder, exp, rankMetric):
		metrics_file = os.path.join( wd, 'metric_scores.tsv')

		df_metrics = exp.get_leaderboard()
		
		cols = list(df_metrics.columns)	
		sel_cols = [ cols[0] ] + cols[2:]
		df_metrics = df_metrics.sort_values( by = self.metric_titles[rankMetric], ascending = False)[ sel_cols ].reset_index(drop=True)
		
		leader = [ subfolder ] + df_metrics.loc[0, sel_cols].values.tolist()
		
		model_names = df_metrics['Model Name']
		model_abbv = [ self.classifier_universe[v] for v in model_names ]
		df_metrics['model_codes'] = model_abbv

		df_metrics.to_csv( metrics_file, sep='\t', index=None )
		
		return df_metrics, leader

	def _choose_save_best_global_model(self, outfolder, rank_analysis, rankMetric):
		sdf = dict( sorted( rank_analysis.items(), key=lambda item: item[1]['metrics'][rankMetric], reverse=True) )
		best_submodel = list(sdf.keys())[0]
		dirout = '/'.join(self.dataDir.split('/')[:-1])
		if( self.dataDir .endswith('/') ):
				dirout = '/'.join(self.dataDir[:-1].split('/')[:-1])
		global_best_model = os.path.join( dirout, outfolder, f'{best_submodel}-global_best' )
		shutil.copy( sdf[best_submodel]['source'], global_best_model )

		return global_best_model, f'{best_submodel}-global_best'

	def _save_plots(self, wd, exp, best, include_explainability = False):
		plotDir = os.path.join( wd, 'plots')
		if( not os.path.isdir(plotDir) ):
			os.mkdir(plotDir)

		dir_bkp = os.getcwd()
		os.chdir( plotDir )
		for p in ['auc', 'pr', 'confusion_matrix', 'feature']:
			try:
				exp.plot_model( best, plot=p, save=True )
			except:
				pass 

		if( include_explainability):
			exp.interpret_model( best, plot='summary', save=True )
		os.chdir(dir_bkp)

	def _run_default_ml_pipeline(self, wd, subfolder, imbalanceMethod = 'downsampling', rankMetric='mcc', has_feat_selection=False, include_explainability=False ):
		exp = self._initialize_experiment( wd, imbalanceMethod, has_feat_selection )
		best, model_file = self._export_subfolder_model( wd, exp, include_explainability)
		df_metrics, leader = self._export_subfolder_scores( wd, subfolder, exp, rankMetric)
		exp_file = os.path.join( wd, 'experiment.pkl')
		exp.save_experiment(exp_file)
		self._save_plots( wd, exp, best, include_explainability)

		return leader, model_file

	def make_feature_selection_general_dataset(self, outfolder, imbalanceMethod = 'downsampling'):
		subfolder = 'all-features'
		wd = os.path.join( outfolder, subfolder )
		if( os.path.isdir(wd) ):
			targetFile = os.path.join( wd, "metric_scores.tsv" )
			if( not os.path.isfile(targetFile) ):
				self._run_default_ml_pipeline( wd, subfolder, imbalanceMethod = imbalanceMethod, has_feat_selection = True, include_explainability = True )

	def build_models_and_ranking(self, outfolder, imbalanceMethod = 'downsampling', rankMetric = 'mcc'):
		rank_analysis = {}
		
		summary_file = self._initialize_summary_table(outfolder)
		for subfolder in os.listdir( outfolder ):
			wd = os.path.join( outfolder, subfolder)
			if( os.path.isdir( wd ) ):
				leader, model_file = self._run_default_ml_pipeline( wd, subfolder, imbalanceMethod = imbalanceMethod, has_feat_selection = False, include_explainability = False )
				
				col_values = [ str(v) for v in leader ]
				with open( summary_file, 'a' ) as g:
					g.write( ('\t'.join(col_values))+'\n' )

				rank_analysis[subfolder] = {}
				rank_analysis[subfolder]['metrics'] = dict( zip( self.metrics, leader[2:] ) )
				rank_analysis[subfolder]['source'] = model_file

		global_model, model_id = self._choose_save_best_global_model( outfolder, rank_analysis, rankMetric)

		return global_model, model_id