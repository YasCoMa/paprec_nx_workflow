import os

import sys
workflow_path = os.environ.get('paprecPath')
if(workflow_path[-1]=='/'):
    workflow_path = workflow_path[:-1]

from method import MetaExtractionMethod

class Implementation_vaxijenModified(MetaExtractionMethod):
    
    def __init__(self):
        self.descriptors = self.build_descriptors()
    
    def build_descriptors(self):
        descriptors={}
        f=open( f"{workflow_path}/modules/pipeline_data/data_features_extraction/descriptors_pmc5549711.tsv","r")
        for line in f:
            l=line.replace("\n","").split("\t")
            aa=l[0].upper()
            zds=[]
            for z in l[1:]:
                zds.append(float(z))
            descriptors[aa]=zds
        f.close()
        
        return descriptors
        
    def _get_optimal_lag(self, seqs):
        n = len(seqs)
        i=30
        init=1
        while(init!=n):
            init=0
            for s in seqs:
                if(len(s)>=i):
                    init+=1
            i-=1    
        if(i<4):
            i=4
        return i
    
    def _calculate_features(self, seq, max_lag):
        seq = seq.upper()
        aac={}
        
        n=len(seq)
        
        l=max_lag
        l = 8
        
        lags=list(range(1, l+1))
        a=[]
        for l in lags:
            aac["l"+str(l)]=[]
            
            for j in range(6):
                mean_ = 0
                for aa in seq:
                    if( aa in self.descriptors ):
                        el = self.descriptors[aa]
                        mean_ += el[j]
                mean_=mean_/n
                
                s=0
                i=0
                while i<(n-l):
                    aa = seq[i]
                    flag = False
                    if(aa in self.descriptors.keys()):
                        e = self.descriptors[aa]
                        va = e[j] - mean_
                        flag = True
                    
                    aal = seq[i+l]
                    flag2 = False
                    if(aal in self.descriptors.keys()):
                        el = self.descriptors[aal]
                        vb = el[j] - mean_
                        flag2 = True
                    
                    if(flag and flag2):
                        s+= (va * vb) / (n-l)
                    
                    i+=1
                        
                aac["l"+str(l)].append(s)
        return aac
    
    def build_numerical_datasets(self, df_input, outfolder):
        seqs = set( df_input['sequence'].tolist() )
        max_lag = self._get_optimal_lag(seqs)

        modes = [ 'unique_with_lags' ]
        for mode in modes:
            for i in df_input.index:
                _id = str( df_input.loc[i, 'id'] ).replace('\t','').replace(' ','').replace(':','').replace('#','').replace('|','_').replace('-','_').replace('>','')
                sequence = df_input.loc[i, 'sequence']
                label = df_input.loc[i, 'label']

                data = self._calculate_features(sequence, max_lag)
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
