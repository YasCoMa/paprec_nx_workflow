import os
import json
import pandas as pd
import numpy as np
import Levenshtein
import statistics as st
from tqdm import tqdm

import sys
workflow_path = os.environ.get('paprecPath')
if(workflow_path[-1]=='/'):
    workflow_path = workflow_path[:-1]

sys.path.insert(0, os.path.join( workflow_path, 'modules', "datasets") )
from treat_dataset import Dataset

"""
inputs :
    - object: { 'saureus': { 'folder_in': 'saureus_seqs.faa', 'path_result_il': 't' } }
    - threshold rank
    - threshold number of alleles for promiscuity
    - threshold score similarity with iedb epitopes
    
scrap il inducer:
    https://webs.iiitd.edu.in/raghava/ifnepitope/predict.php
    https://webs.iiitd.edu.in/raghava/il4pred/predict.php
    https://webs.iiitd.edu.in/raghava/il10pred/predict3.php
    http://metagenomics.iiserb.ac.in/IL17eScan/pred.php
"""

class ParserCuration:
    def __init__(self, outfolder, dataset_id, protein_file):
        self.ds = Dataset( outfolder, dataset_id )
        self.wd = os.path.join( self.ds.dataset_folder, 'data_parsing_selection')
        if(not os.path.isdir(self.wd) ):
            os.mkdir(self.wd)

        self._load_proteins(protein_file)
        self._initialize_report_files()

    def _load_proteins(self, protein_file):
        dt = {}
        g = open( protein_file,'r')
        for line in g:
            l=line.replace('\n','')
            if(l.startswith('>')):
                k = l.replace('>','')
                k = k.lower().replace('/','').replace(':','').replace('+','').replace('(','').replace(')','').replace('[','').replace(']','').replace('|','_').replace(' ','_')
                aux = list( filter( lambda x: x!='', k.split('_') ) )
                k = aux[0]
                if( len(aux) > 1):
                    k+='_'+aux[1]
                dt[k]=""
            else:
                dt[k]+=l
        g.close()

        self.proteins = dt

    def _initialize_report_files(self):
        self.log_file = os.path.join( self.wd, 'log_curation.tsv' )
        self.removed_file = os.path.join( self.wd, 'removed_sequences.tsv' )

        if( not os.path.isfile(self.log_file) ):
            f = open( self.log_file, 'w')
            f.write('step\titems_removed\n')
            f.close()

        if( not os.path.isfile(self.removed_file) ):
            f = open( self.removed_file, 'w')
            f.close()

    def _update_log_file(self, step, n_removed):
        with open( self.log_file, 'a') as g:
            g.write( f'{step}\t{n_removed}\n' )

    def _update_removed_file(self, peptides):
        gone = self._load_removed_peptides()

        epitopes =  set(peptides) - gone
        if( len(epitopes) > 0 ):
            epitopes = '\n'.join(epitopes)
            with open( self.removed_file, 'a') as g:
                g.write( f'{epitopes}\n' )

    def _load_removed_peptides(self):
        gone = set()
        if( os.path.isfile(self.removed_file) ):
            gone = set( open(self.removed_file).read().split('\n') )
        return gone

    def parse_prediction_results(self, raw_file):
        dir_out = self.wd
        
        outfile = os.path.join( dir_out, 'table_parser_results.tsv' )
        gone=set()
        """
        if( os.path.isfile( outfile ) ):
            i=0
            gf = open( outfile,'r')
            for line in f:
                if(i>0):
                    l=line.split('\t')
                    ide = f'{ l[0] }-{ l[1] }'
                    gone.add( ide )
                i+=1
            gf.close()
        else:
        """
        
        if( not os.path.isfile( outfile ) ):
            gf = open( outfile,'w')
            gf.write('protein\tmhc\tpeptide\tcore\tscore_ba\taffinity\trank%\n')
            gf.close()
            
            g = open( raw_file, 'r')
            for line in g:
                line=line.replace('\n','')
                if(line.replace(' ','').lower().find('<=wb')!=-1):
                    l = list( filter( lambda x: x!='', line.split(' ') ))

                    k = l[7]
                    k = k.lower().replace('/','').replace(':','').replace('+','').replace('(','').replace(')','').replace('[','').replace(']','').replace('|','_').replace(' ','_')
                    aux = list( filter( lambda x: x!='', k.split('_') ) )
                    k = aux[0]
                    if( len(aux) > 1):
                        k += '_'+aux[1]
                    prot = k
                    mhc = l[1]
                    pep = l[2]
                    core = l[4]
                    scoreba = l[11]
                    aff = l[13]
                    rank = l[12]
                    ide = f'{prot}-{mhc}-{pep}'

                    if(not ide in gone):
                        with open( outfile,'a') as gf:
                            gf.write(f'{prot}\t{mhc}\t{pep}\t{core}\t{scoreba}\t{aff}\t{rank}\n')
                    
            g.close()

        df = pd.read_csv( outfile, sep='\t' )

        return df
            
    def filter_rank(self, cutoff):
        dir_out = self.wd
        outfile = os.path.join( dir_out, 'table_parser_results.tsv' )
        df = pd.read_csv( outfile, sep='\t' )

        step = "Filter by rank percentile BA"
        epi_out=set()
        
        outfile = os.path.join( dir_out, 'filtered_by_rank.tsv' )
        if( not os.path.isfile( outfile ) ):
            f = open( outfile, 'w' )
            f.write( f"sequence\trank\n" )
            for i in df.index:
                epi = df.loc[i, 'peptide']
                rank = df.loc[i, 'rank%']
                if( rank > cutoff ):
                    epi_out.add(epi)
                    f.write( f"{epi}\t{rank}\n" )
            f.close()

            self._update_log_file(step, len(epi_out) )
            self._update_removed_file(epi_out)

    def check_allele_promiscuity(self, cutoff):
        dir_out = self.wd
        outfile = os.path.join( dir_out, 'table_parser_results.tsv' )
        df = pd.read_csv( outfile, sep='\t' )

        step = "Filter by allele promiscuity"
        
        outfile = os.path.join( dir_out, 'promiscuous_epitopes.tsv' )
        if( not os.path.isfile( outfile ) ):
            epi_out=set()
            
            dt={}
            for i in df.index:
                epi=df.loc[i, 'peptide']
                mhc=df.loc[i, 'mhc']
                if( not epi in dt):
                    dt[epi]=set()
                dt[epi].add(mhc)
            
            f=open( outfile, "w")
            f.write('sequence\tnum_alleles\talleles\n')
            for e in dt:
                n = len(dt[e])
                if( n > cutoff ):
                    epi_out.add(e)
                    mhcs = ','.join( list(dt[e]) )
                    f.write(f"{e}\t{len(dt[e])}\t{mhcs}\n")
            f.close() 
            
            self._update_log_file(step, len(epi_out) )
            self._update_removed_file(epi_out)
    
    def check_overlapping_epis_violinet(self):
        dir_out = self.wd
        outfile = os.path.join( dir_out, 'table_parser_results.tsv' )
        df = pd.read_csv( outfile, sep='\t' )

        step = "Filter by overlapping in violinet db"
        
        outfile = os.path.join( dir_out, 'overlapping_violinet.tsv' )
        if( not os.path.isfile( outfile ) ):
            epi_out=set()
            
            epis = set( df['peptide'].tolist() )
            gone = self._load_removed_peptides()
            epis =  set(epis) - gone
            
            gone=set()
            seq=''
            g=open( outfile,'w')
            g.write('protein_id\tprotein_sequence\tsequence\n')
            
            f=open( f"{workflow_path}/modules/pipeline_data/data_filters_curation/protegen-bacterium.faa","r")
            for line in f:
                if(line!='\n'):
                    l=line.replace('\n','').replace('>','')
                    if(line.startswith('>')):
                        if(seq!=''):
                            for e in epis:
                                ide = e+'--'+id_
                                if(seq.find(e)!=-1 and not ide in gone ):
                                    gone.add(ide)
                                    epi_out.add(e)
                                    g.write('%s\t%s\t%s\n' %(id_, seq, e) )
                        id_=line.split('|')[0]+'_'+line.split('|')[1]
                        seq=''
                    else:
                        seq+=l
            f.close()
            
            if(seq!=''):
                for e in epis:
                    ide = e+'--'+id_
                    if(seq.find(e)!=-1 and not ide in gone ):
                        gone.add(ide)
                        epi_out.add(e)
                        g.write('%s\t%s\t%s\n' %(id_, seq, e) )
            g.close()
            
            self._update_log_file(step, len(epi_out) )
            self._update_removed_file(epi_out)
        
    def filter_iedb_epitopes(self, cutoff):
        dir_out = self.wd
        outfile = os.path.join( dir_out, 'table_parser_results.tsv' )
        df = pd.read_csv( outfile, sep='\t' )

        step = "Filter by overlapping in iedb db"
        
        outfile = os.path.join( dir_out, 'overlapping_violinet.tsv' )
        if( not os.path.isfile( outfile ) ):
            epi_out=set()
            
            epis = set( df['peptide'].tolist() )
            gone = self._load_removed_peptides()
            epis =  set(epis) - gone

            pubepis={}
            i=0
            f=open( f"{workflow_path}/modules/pipeline_data/data_filters_curation/epitope_table_iedb.tsv","r")
            for line in f:
                l=line.replace("\n","").split("\t")
                if(i>0):
                    ide = l[0].split('/')[-1]
                    epi = l[2]
                    org = l[-4].replace('"','')
                    pubepis[epi] = [ide, org]
                i+=1
            f.close()

            report={}
            details={}
            for epi in epis:
                info=[]

                similar=0
                for target in pubepis.keys():
                    identity = Levenshtein.ratio(epi, target)
                    if(identity >= cutoff):
                        epi_out.add(epi)
                        
                        similar+=1
                        inf = pubepis[target]
                        inf.append( str(identity*100) )
                        
                        info.append( '-'.join( inf ) )

                report[epi]=similar
                details[epi]=info

            sorted_list = sorted( report.items(), key=lambda kv: kv[1], reverse=True )  
            
            g = open( outfile, "w")
            g.write('sequence\tnumber_matches\tepitope_iedb_info\n') 
            g.close()     
            for epi in sorted_list:
                with open( outfile, "a") as gf:
                    gf.write("%s\t%s\t%s\n" %(epi[0], epi[1], ",".join(details[epi[0]]) ) )
            
            self._update_log_file(step, len(epi_out) )
            self._update_removed_file(epi_out)
                
    def check_human_homology(self):
        dir_out = self.wd
        outfile = os.path.join( dir_out, 'table_parser_results.tsv' )
        df = pd.read_csv( outfile, sep='\t' )
        epitopes = set( df['peptide'].tolist() )

        step = "Filter by overlapping in human proteins"
        
        outfile = os.path.join( dir_out, 'result_overlapping_human_genome.tsv' )
        """
        else:
            df = pd.read_csv( outfile, sep='\t' )
            epitopes = epitopes - set( df['sequence'].tolist() )
        """
        if( not os.path.isfile( outfile ) ):
            f = open( outfile, "w")
            f.write('sequence\tnumber_proteins\tid_human_protein\n') 
            f.close()
        
            epi_out=set()
            
            gone = self._load_removed_peptides()
            epitopes =  set(epitopes) - gone
            
            humanprot = set()
            aux = {}
            f = open( f"{workflow_path}/modules/pipeline_data/data_filters_curation/uniprotkb_proteome_UP000005640_2024_09_26.fasta","r")
            for line in f:
                l = line.replace("\n","")
                if( l.startswith(">") ):
                    id_ = l.replace('>','')
                    aux[id_] = ''
                else:
                    aux[id_] += l
            f.close()

            mp = {}
            for k in aux:
                seq = aux[k]
                if( not seq in humanprot ):
                    humanprot.add( seq )
                    mp[seq] = set()
                mp[seq].add(k)
                
            epis={}
            for l in epitopes:
                if(not l in epis):
                    epis[l]=set()

            for e in tqdm(epitopes):
                for hs in humanprot:
                    if(e in hs):
                    #if(hs.find(e)!=-1):
                        epi_out.add(e)
                        #epis[e].add( mp[hs] )
                        for k in mp[hs]:
                            with open( outfile, "a") as f:
                                f.write("%s\t%i\t%s\n" %( e, -1, k ) )
                        
            """
            for ep in epis.keys():
                if( len(epis[ep]) > 0):
                    epi_out.add(ep)
                    ids = ','.join( list(epis[ep]) )
                    f.write("%s\t%i\t%s\n" %( ep, len(epis[ep]), ids ) )
            f.close()  
            """

            self._update_log_file(step, len(epi_out) )
            self._update_removed_file(epi_out)
                
    def check_human_homology_blast(self):
        dir_out = self.wd
        outfile = os.path.join( dir_out, 'table_parser_results.tsv' )
        df = pd.read_csv( outfile, sep='\t' )

        step = "Filter by overlapping in human proteins"
        
        outfile = os.path.join( dir_out, 'result_overlapping_human_genome.tsv' )

        epi_out=set()
        
        epitopes = set( df['peptide'].tolist() )
        gone = self._load_removed_peptides()
        epitopes =  set(epitopes) - gone
        
        f = open( outfile, "w"  )
        f.write('sequence\tnumber_human_proteins\tid_human_proteins\n')
        for l in epitopes:
            g = open( f'{self.wd}/query.fasta', 'w')
            g.write( f'>query\n{ l }\n' )
            g.close()
            
            os.system( f"blastp -task blastp-short -db {workflow_path}/modules/pipeline_data/data_filters_curation/human -query {self.wd}/query.fasta -out {self.wd}/temp.out -outfmt 6")
            matches_epis = set()
            x = open( f"{self.wd}/temp.out","r")
            for line in x:
                l1 = line.replace("\n","").split("\t")
                if( float(l1[2]) > 99 ):
                    matches_epis.add( l1[1] )
            x.close()
            
            if( len(matches_epis) > 0 ):
                epi_out.add(l)
                ids = ','.join( list( matches_epis ) )
                f.write("%s\t%i\t%s\n" %( l, len( matches_epis ), ids ) )
        f.close()  
        
        self._update_log_file(step, len(epi_out) )
        self._update_removed_file(epi_out)
        
    def _write_selected_sequences(self, df, final_epis):
        dir_out = self.wd

        edc = {}

        outfile = os.path.join( dir_out, 'selected_epitopes.fasta' )
        j=1
        f=open( outfile, "w")
        for e in final_epis:
            f.write( f">peptide_{j}\n{e}\n")
            edc[e] = f"peptide_{j}"

            j+=1
        f.close()
        
        outfile = os.path.join( dir_out, 'selected_proteins.fasta' )
        f=open( outfile, "w")  
        prots=set()
        for i in df.index:
            epi = df.loc[i, 'peptide']
            protein = df.loc[i, 'protein']
            
            if( (epi in final_epis) and (not protein in prots) ):
                prots.add(protein)
                f.write( f">{protein}\n{ self.proteins[protein] }\n")
        f.close()  

        return edc 

    def prepare_final_sequences(self):
        dir_out = self.wd
        outfile = os.path.join( dir_out, 'table_parser_results.tsv' )
        df = pd.read_csv( outfile, sep='\t' )
        
        epitopes = set( df['peptide'].tolist() )
        gone = self._load_removed_peptides()
        final_epis =  set(epitopes) - gone
            
        epidc = self._write_selected_sequences( df, final_epis)
        
        outfile = os.path.join( dir_out, 'table_final_proteins_epitopes.tsv' )
        dff = {}
        for i in df.index:
            peptide = df.loc[i, 'peptide']

            if(peptide in final_epis):
                if(not peptide in dff):
                    dff[peptide] = { 'mhc': set(), 'protein': set() }
                    
                protein = df.loc[i, 'protein']
                dff[peptide]['protein'].add(protein)
                
                allele = df.loc[i, 'mhc']
                dff[peptide]['mhc'].add(allele)
        
        f=open( outfile, "w")
        f.write("epitope_id\tepitope\tprotein\tmhc_alleles\n")
        for e in final_epis:
            mhc = ','.join( list(dff[e]['mhc']) )
            protein = ','.join( list(dff[e]['protein']) )
            f.write( f"{epidc[e]}\t{e}\t{protein}\t{mhc}\n")
        f.close()
                