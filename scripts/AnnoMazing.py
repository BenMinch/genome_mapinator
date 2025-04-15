###Annotation Tool

import sys,argparse,os,subprocess,re
import pandas as pd
from collections import defaultdict
from Bio import SearchIO
#setup flags
argparser = argparse.ArgumentParser(description='''Script to use Annotate Proteins with multiple databases''')
argparser.add_argument('-i', '--input', help='''Input protein file''', required=True)
argparser.add_argument('-o', '--output', help='''Output name''', required=True)
argparser.add_argument('-db', '--databases', help='''database names separated by a comma''', required=True)
argparser.add_argument('-annot', '--annotation', help='''Annotation files (csv), separated by comma (make sure identifyer is "id")''', required=True)


args = argparser.parse_args()
input_file = args.input
output_file = args.output
databases = args.databases
annotation = args.annotation

#separate databases
databases=databases.split(',')
number_of_databases=len(databases)
annotation=annotation.split(',')
number_of_annotation=len(annotation)
if number_of_databases != number_of_annotation:
    print('Number of databases and annotation files do not match')
    sys.exit()
#Do Hmmscan for databases
for i in range(0,number_of_databases):
    database1=databases[i]
    hi= str(i)
    print('Running HMMscan for '+database1)
    hmmscan='hmmscan -E 0.00001 --cpu 12 --tblout '+output_file+'_'+hi+'.hmmtable '+ database1+ ' '+input_file
    subprocess.call(hmmscan, shell=True)

#Parse best hits
print('Parsing best hits')
for file in os.listdir('.'):
    if file.endswith('.hmmtable'):
        new_name=re.sub('.hmmtable', '_parsed.tsv', file)
        parse="awk '!x[$3]++' "+file+' > '+new_name
        subprocess.call(parse, shell=True)

for file in os.listdir('.'):
    if file.endswith('_parsed.tsv'):
        new_name=re.sub('_parsed.tsv', '_parsinator.tsv', file)
        attribs = ['accession', 'bias', 'bitscore', 'description', 'cluster_num', 'domain_exp_num',  'domain_included_num', 'domain_obs_num', 'domain_reported_num', 'env_num', 'evalue', 'id', 'overlap_num', 'region_num']
        hits = defaultdict(list)
        with open(file) as handle:
            for queryresult in SearchIO.parse(handle, 'hmmer3-tab'):
      #put query id in the dictionary
                hits['query_id'].append(queryresult.id)
      #print(queryresult.accession)
      #print(queryresult.description)
                for hit in queryresult.hits:
                    for attrib in attribs:
                        hits[attrib].append(getattr(hit, attrib))

        pd.DataFrame.from_dict(hits).to_csv(new_name, index=False)

for i in range(0,number_of_databases):
    database1=databases[i]
    annotation1=annotation[i]
    for file in os.listdir('.'):
        if file.endswith('_parsinator.tsv'):
            if str(i) in file:
                new_name=re.sub('_parsinator.tsv','_annotated.tsv', file)
                input_df=pd.read_csv(file)
                #determine if annotation_df is tsv or csv
                if annotation1.endswith('.tsv'):
                    annotation_df=pd.read_csv(annotation1, sep='\t')
                else:
                    annotation_df=pd.read_csv(annotation1)
                merged_df=pd.merge(input_df, annotation_df, on='id', how='left')
                #remove columns 2-11 from merged_df
                merged_df=merged_df.drop(merged_df.columns[[1,2,3,4,5,6,7,8,9,10,13,14]], axis=1)
                merged_df.to_csv(new_name,index=False)

#combining all the files
import glob
file_list=glob.glob('*_annotated.tsv')
df_list=[]

merged_df=pd.read_csv(file_list[0])
for csv_file in file_list[1:]:
    df=pd.read_csv(csv_file)
    merged_df=merged_df.merge(df, on='query_id', how='outer')
merged_df.to_csv(output_file+'_final.csv',index=False)


#Remove all the intermediate files
remove='rm *parsed.tsv *parsinator.tsv *hmmtable *annotated.tsv'
subprocess.call(remove, shell=True)
