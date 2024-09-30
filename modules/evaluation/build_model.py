import os
import joblib
import shutil
import pandas as pd
from pycaret.classification import *

class BuildBaseModel:
	
	def _initialize_summary_table(self, outfolder):
		summary_file = os.path.join( outfolder, 'summary_best_models.tsv' )
		metrics = ['accuracy', 'auc', 'recall', 'precision', 'f1', 'kappa', 'mcc']
		header = ['type_analysis', 'model' ] + metrics
		f = open( summary_file, 'w' )
		f.write( ('\t'.join(header))+'\n' )
		f.close()

		return summary_file, metrics

	def _initialize_experiment(self, outfolder, subfolder, imbalanceMethod = 'downsampling', has_feat_selection = False ):
		fix_classes = False 
		if( imbalanceMethod == 'smote'):
			fix_classes = True

		exp = ClassificationExperiment()
		infile = os.path.join( outfolder, subfolder, "computed_dataset.tsv")
		data = pd.read_csv( infile, sep='\t' )
		exp.setup(data = data, target = 'label', fix_imbalance = fix_classes, feature_selection = has_feat_selection session_id = 42)
		
		if( has_feat_selection ):
			features_before = data.columns
			features_after = exp.dataset_transformed.columns
			for state in ['before', 'after']:
				outfile = os.path.join( outfolder, subfolder, f"features_{state}.txt"  )
				f = open(outfile, 'w')
				f.write( '\n'.join( list( eval( f"features_{state}") ) )  )
				f.close()
		return exp

	def _export_subfolder_model(self, outfolder, subfolder, exp):
		model_file = os.path.join(outfolder, subfolder, 'best_scored.model')
		best = exp.compare_models( exclude = ['lr', 'ridge', 'catboost', 'dummy'] )
		joblib.dump(best, model_file)

		return best, model_file

	def _export_subfolder_scores(self, outfolder, subfolder, exp):
		metrics_file = os.path.join(outfolder, subfolder, 'metric_scores.tsv')
		df_metrics = exp.get_leaderboard()
		cols = list(df_metrics.columns)
		sel_cols = [ cols[0] ] + cols[2:]
		df_metrics = df_metrics[ sel_cols ]
		df_metrics.to_csv( metrics_file, sep='\t', index=None )
		leader = [ subfolder ] + df_metrics.loc[0, sel_cols].values.tolist()

		return df_metrics, leader

	def _choose_save_best_global_model(self, outfolder, rank_analysis, rankMetric):
		sdf = dict( sorted( rank_analysis.items(), key=lambda item: item[1]['metrics'][rankMetric], reverse=True) )
		best_submodel = list(sdf.keys())[0]
		global_best_model = os.path.join( outfolder, f'{best_submodel}-global_best_model.model' )
		shutil.copy( sdf[best_submodel]['source'], global_best_model )

	def _save_plots(self, exp, best, include_explainability = False):
		plotDir = os.path.join( outfolder, subfolder, 'plots')
		if( not os.path.isdir(plotDir) ):
			os.mkdir(plotDir)

		dir_bkp = os.getcwd()
		os.chdir( plotDir )
		for p in ['auc', 'pr', 'confusion_matrix', 'feature']:
			exp.plot_model( best, plot=p, save=True )

		if( include_explainability):
			exp.interpret_model( plot='summary', save=True )
		os.chdir(dir_bkp)

	def _run_default_ml_pipeline(self, outfolder, subfolder, imbalanceMethod = 'downsampling', has_feat_selection=False, include_explainability=False ):
		exp = self._initialize_experiment(outfolder, subfolder, imbalanceMethod, has_feat_selection )
		best, model_file = self._export_subfolder_model( outfolder, subfolder, exp)
		df_metrics, leader = self._export_subfolder_scores( outfolder, subfolder, exp)
		exp_file = os.path.join(outfolder, subfolder, 'experiment.pkl')
		exp.save_experiment(exp_file)
		self._save_plots(exp, best, include_explainability)

		return leader, model_file

	def make_feature_selection_general_dataset(self, outfolder, imbalanceMethod = 'downsampling'):
		subfolder = 'all-features'
		wd = os.path.join( outfolder, subfolder )
		if( os.path.isdir(wd) ):
			self._run_default_ml_pipeline( wd, imbalanceMethod = imbalanceMethod, has_feat_selection = True, include_explainability = True )


	def build_models_and_ranking(self, outfolder, imbalanceMethod = 'downsampling', rankMetric = 'mcc'):
		rank_analysis = {}
		
		summary_file, metrics = self._initialize_summary_table(outfolder)

		for subfolder in os.listdir( outfolder ):
			if( os.path.isdir( os.path.join( outfolder, subfolder) ) ):
				leader, model_file = self._run_default_ml_pipeline( outfolder, subfolder, imbalanceMethod = imbalanceMethod, has_feat_selection = False, include_explainability = False )
				
				col_values = [ str(v) for v in leader ]
				with open( summary_file, 'a' ) as g:
					g.write( ('\t'.join(col_values))+'\n' )

				rank_analysis[subfolder] = {}
				rank_analysis[subfolder]['metrics'] = dict( zip( metrics, leader[2:] ) )
				rank_analysis[subfolder]['source'] = model_file

		self._choose_save_best_global_model( outfolder, rank_analysis, rankMetric)