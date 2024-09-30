import os
import shutil
import pandas as pd

class MetaExtractionMethod:
    
    def _handle_subanalysis_file(self, outfolder, identifier, n_features = 0, columns = [] ):
        subfolder = os.path.join(outfolder, identifier)
        outfile = os.path.join( subfolder, "computed_dataset.tsv")

        if( not os.path.isdir( subfolder ) ):
            os.mkdir( subfolder )

            header = [ "item_id", "item_sequence", "label" ]
            if( n_features > 0 ):
                for i in range( 1, n_features+1 ):
                    header.append( f"feat_{identifier}_{i}" )
                
            elif( len(columns) != 0 ):
                header += columns
            header = '\t'.join(header)

            f = open( outfile, 'w' )
            f.write( f"{ header }\n" )
            f.close()

        return outfile

    def _get_all_feature_columns(self, outfolder):
        feature_columns = [  ]
        for subfolder in os.listdir( outfolder ):
            if( os.path.isdir( os.path.join( outfolder, subfolder) ) ):
                infile = os.path.join( outfolder, subfolder, "computed_dataset.tsv")
                header = open( infile, 'r').readline().replace('\n','').split('\t')
                feature_columns += list(header)[3:]
        return feature_columns

    def _build_combined_ds_for_feature_selection(self, outfolder):
        nfile = 0
        feature_columns = self._get_all_feature_columns(outfolder)
        outfile = self._handle_subanalysis_file( outfolder, "all-features", columns = feature_columns)
        tmpTransfer = os.path.join( outfolder, "all_feats.tsv")

        for subfolder in os.listdir( outfolder ):
            if( os.path.isdir( os.path.join( outfolder, subfolder) ) and (subfolder != 'all-features') ):
                idx = 0
                info = {}
                new_col_names = []
                infile = os.path.join( outfolder, subfolder, "computed_dataset.tsv")
                f = open( infile, 'r')
                for line in f:
                    l = line.replace('\n','').split('\t')
                    if(idx == 0):
                        if(nfile == 0):
                            new_col_names = l
                        else:
                            new_col_names = l[3:]
                    else:
                        _id = l[0]
                        label = l[2]
                        sk = f"{_id}#{label}"
                        info[sk] = '\t'.join(l[3:])

                    idx+=1
                f.close()

                tempf = open( tmpTransfer, 'w' )
                if( nfile == 0):
                    colnames = ['item_id', 'label'] + new_col_names
                    newline = ('\t'.join(colnames))+'\n'
                    tempf.write(newline)

                    for sk in info:
                        _id, label = sk.split('#')
                        newline = ('\t'.join( [_id, label, info[sk] ] ))+'\n'
                        tempf.write(newline)
                else:
                    idx = 0
                    f = open(outfile, 'r')
                    for line in f:
                        l = line.replace('\n','').split('\t')
                        if(idx==0):
                            colnames = l + new_col_names
                            newline = ('\t'.join(colnames))+'\n'
                            tempf.write(newline)
                        else:
                            _id = l[0]
                            label = l[1]
                            sk = f"{_id}#{label}"
                            newline = ('\t'.join( l + [ info[sk] ] ))+'\n'
                            tempf.write(newline)
                        idx+=1
                    f.close()
                tempf.close()

                shutil.move( tmpTransfer, outfile )

                nfile += 1

        """
        info = {}
        feature_columns = [  ]
        for subfolder in os.listdir( outfolder ):
            infile = os.path.join( outfolder, subfolder, "computed_dataset.tsv")
            df = pd.read_csv( infile, sep='\t' )

            feature_columns += list(df.columns)[3:]

            for i in df.index:
                _id = str( df.loc[i, 'item_id'] )
                sequence = df.loc[i, 'item_sequence']
                label = df.loc[i, 'label']

                sk = f"{_id}#{sequence}#{label}"
                if( not sk in ides ):
                    ides.add(sk)
                    info[ sk ] = []
                for v in list(df.columns)[3:]:
                	info[ sk ].append( df.loc[i, v] )

        outfile = self._handle_subanalysis_file( outfolder, "all-features", columns = feature_columns)
        for key in info:
            print(key)
            _id, sequence, label = key.split('#')
            features = info[key]
            features = [ str(v) for v in features ]
            features = '\t'.join(features)

            with open( outfile, 'a' ) as f:
                f.write( "%s\t%s\t%s\t%s\n" %( _id, sequence, label, features) )
        """
       