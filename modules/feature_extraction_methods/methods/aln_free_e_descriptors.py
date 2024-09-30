import os

import sys
workflow_path = os.environ.get('paprecPath')
if(workflow_path[-1]=='/'):
    workflow_path = workflow_path[:-1]

from method import MetaExtractionMethod

class Implementation_vaxijen(MetaExtractionMethod):
    
    def __init__(self):
        self.descriptors = self.build_descriptors()
    
    def build_descriptors(self):
        descriptors={}
        f=open( f"{workflow_path}/modules/pipeline_data/data_features_extraction/e_table_descriptors.tsv","r")
        for line in f:
            l=line.replace("\n","").split(" ")
            aa=l[0].upper()
            zds=[]
            for z in l[1:]:
                zds.append(float(z))
            descriptors[aa]=zds
        f.close()
        
        return descriptors
    
    def _calculate_features(self, seq, mode):
        seq = seq.upper()
        aac={}
        
        n=len(seq)
        
        l=8
        lags=list(range(1, l+1))
        a=[]
        for l in lags:
            aac["l"+str(l)]=[]
            
            for i in range(5):
                for j in range(5):
                    if(mode=='auto'):
                        cond=(i==j)
                    if(mode=='cross'):
                        cond=(i!=j)
                        
                    if(cond):
                        s=0
                        
                        k=0
                        while k<(n-l):
                            aa=seq[k]
                            flag=False
                            if(aa in self.descriptors.keys()):
                                e=self.descriptors[aa]
                                flag=True
                            
                            aal=seq[k+l]
                            flag2=False
                            if(aal in self.descriptors.keys()):
                                el=self.descriptors[aal]
                                flag2=True
                            
                            if(flag and flag2):
                                s+= (e[j]*el[j]) / (n-l)
                            
                            k+=1
                        
                        aac["l"+str(l)].append(s)
        return aac

    def build_numerical_datasets(self, df_input, outfolder):
        modes = [ 'auto', 'cross' ]
        for mode in modes:
            for i in df_input.index:
                _id = str( df_input.loc[i, 'id'] ).replace('\t','').replace(' ','').replace(':','').replace('-','_')
                sequence = df_input.loc[i, 'sequence']
                label = df_input.loc[i, 'label']

                data = self._calculate_features(sequence, mode)
                for lag in data:
                    identifier = f"{mode}-{lag}"
                    features = data[lag]
                    n_features = len(features)
                    outfile = self._handle_subanalysis_file( outfolder, identifier, n_features = n_features)
                    
                    features = [ str(v) for v in features ]
                    features = '\t'.join(features)

                    with open( outfile, 'a' ) as f:
                        f.write( "%s\t%s\t%i\t%s\n" %( _id, sequence, label, features) )

        self._build_combined_ds_for_feature_selection( outfolder )
