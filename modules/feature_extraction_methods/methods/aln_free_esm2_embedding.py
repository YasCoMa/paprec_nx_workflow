import os
import torch
from transformers import EsmTokenizer, EsmModel

import sys
workflow_path = os.environ.get('paprecPath')
if(workflow_path[-1]=='/'):
    workflow_path = workflow_path[:-1]

from method import MetaExtractionMethod

class Implementation_es2Embedding(MetaExtractionMethod):
    
    def __init__(self):
        self.tokenizer = EsmTokenizer.from_pretrained("facebook/esm2_t6_8M_UR50D")
        self.model = EsmModel.from_pretrained("facebook/esm2_t6_8M_UR50D")
    
    def _calculate_features(self, seqs):
        inputs = self.tokenizer(seqs, return_tensors = "pt", padding = True, truncation = True, max_length = 2048)
        outputs = self.model(**inputs)
        last_hidden_states = outputs.last_hidden_state

        x = last_hidden_states.detach()
        result = [ v.mean(axis=0) for v in x ]

        return result
    
    def build_numerical_datasets(self, df_input, outfolder):
        modes = [ 'es2_model' ]
        for mode in modes:
            ids = [ str(x).replace('\t','').replace(' ','').replace("'",'').replace('[','').replace(']','').replace(':','').replace('-','_') for x in df_input['id'].tolist() ]
            sequences = df_input['sequence'].tolist()
            labels = df_input['label'].tolist()

            i = 0
            for s in sequences:
                _id = ids[i]
                sequence = [s]
                label = labels[i]

                identifier = mode
                features = self._calculate_features(sequence)[0].tolist()
                n_features = len(features)
                outfile = self._handle_subanalysis_file( outfolder, identifier, n_features = n_features)

                features = [ str(v) for v in features ]
                features = '\t'.join(features)

                with open( outfile, 'a' ) as f:
                    f.write( "%s\t%s\t%i\t%s\n" %( _id, sequence, label, features) )
                i+=1

        #self._build_combined_ds_for_feature_selection( outfolder )
