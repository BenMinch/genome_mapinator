import numpy as np
import pandas as pd
import os, sys, argparse, subprocess
import re

def parse_args():
    parser = argparse.ArgumentParser(description="This script will edit multiple gff files with fasta information and add sequence-region headers.")
    parser.add_argument("-i", "--input_dir", type=str, help="Directory containing Genomes")
    parser.add_argument("-o", "--output", type=str, help="Directory to save output")
    return parser.parse_args()
args = parse_args()
input_dir = args.input_dir
output_dir = args.output
os.makedirs(output_dir)
#Make a function to process the input directory to make all contig headers be the file name
def process_fasta_headers(input_dir):
    for filename in os.listdir(input_dir):
        if filename.endswith(".fna") or filename.endswith(".fasta"):
            file_path = os.path.join(input_dir, filename)
            with open(file_path, 'r') as file:
                lines = file.readlines()
            # Modify the header lines
            new_lines = []
            for line in lines:
                if line.startswith('>'):
                    # Replace the header with the filename without extension
                    new_header = f">{os.path.splitext(filename)[0]}\n"
                    new_lines.append(new_header)
                else:
                    new_lines.append(line)
            # Write the modified lines back to the file
            with open(file_path, 'w') as file:
                file.writelines(new_lines)
process_fasta_headers(input_dir)
def run_prodigal(input_dir, output_dir):
    os.mkdir(output_dir+"_Prodigal")
    #run prodigal script from scripts folder
    prodigal= 'python scripts/prodigal_launcher.py '+input_dir+' '+output_dir+'_Prodigal'
    subprocess.call(prodigal, shell=True)
print('Running Prodigal!')
run_prodigal(input_dir, output_dir)

#Combine all .faa files
combine= 'cat '+output_dir+'_Prodigal/*.faa > '+output_dir+'_Prodigal/combined.faa'
subprocess.call(combine, shell=True)
#move these to the output directory
move= 'mv '+output_dir+'_Prodigal/combined.faa '+output_dir
subprocess.call(move, shell=True)
#Create a folder called gff_files
os.makedirs(output_dir+"_gff_files")
#move all gff files to the gff_files folder
move= 'mv '+output_dir+'_Prodigal/*.gff '+output_dir+'_gff_files'
subprocess.call(move, shell=True)

##Run annomazing on the combined.faa file
print('Running AnnoMazing!')
annomazing='python scripts/AnnoMazing.py -i '+output_dir+'/combined.faa -o '+output_dir+'_annomazing -db hmms/Pfam.hmm -annot resources/Pfam_annotations.csv'
subprocess.call(annomazing, shell=True,stderr=open("anomazing.err", "w"),stdout=open("anomazing.out", "w"))

##Running the gff modificaiton script##
#Create a csv file with gff filenames, fasta filenames and output filename which is the same as the gff filename with _final after

def create_gff_csv(input_dir, output_dir):
    # Create a list to store the data
    data = []
    # Iterate through the files in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith(".gff"):
            # Get the corresponding fasta filename
            fasta_filename = filename.replace(".gff", "")
            # Create the output filename
            output_filename = filename.replace(".gff", "_final.gff")
            # Append the data to the list
            data.append([filename, fasta_filename, output_filename])
    # Create a DataFrame from the data
    df = pd.DataFrame(data, columns=["gff", "fasta", "output"])
    # Save the DataFrame to a CSV file
    csv_filename = os.path.join(output_dir, "gff_fasta_mapping.csv")
    df.to_csv(csv_filename, index=False)
    return csv_filename
gff_folder =output_dir+"_gff_files"
create_gff_csv(gff_folder, output_dir)

#Run the gffcombiner script
gff_out= output_dir+'_gff_files_final'
os.makedirs(gff_out)
gffcomb='python scripts/gff_combiner.py -csv '+output_dir+'/gff_fasta_mapping.csv -gff_dir '+ gff_folder+' -fasta_dir '+input_dir+ ' -output_dir '+gff_out
subprocess.call(gffcomb, shell=True)

###Run Lovis4U with the gff files###
lovis_1= 'lovis4u -gff '+gff_out+ ' --set-category-colour -c A4p2 -o lovis_out_1 -alip -hl -ufid'
subprocess.call(lovis_1, shell=True)

##Matching feature annotations
feature_annot=pd.read_csv('lovis_out_1/feature_annotation_table.tsv', sep='\t')
pfam_annot=pd.read_csv(output_dir+'_annomazing_final.csv')
#remove .fna_final-1 from feature_id column
feature_annot['feature_id'] = feature_annot['feature_id'].str.replace('.fna_final-1', '')
feature_annot['feature_id'] = feature_annot['feature_id'].str.replace('.fasta_final-1', '')
#map pfam DE column to feature annotation[name]
feature_annot['name']=feature_annot['feature_id'].map(pfam_annot.set_index('query_id')['DE'])
#If the name column is not NA, set show_label to 1
feature_annot['show_label'] = np.where(feature_annot['name'].notna(), 1, 0)
#Generate a random color for each unique name
unique_names = feature_annot['name'].unique()
#remove NA values from unique_names
unique_names = unique_names[~pd.isna(unique_names)]
colors = {}
for name in unique_names:
    colors[name] = "#{:06x}".format(np.random.randint(0, 0xFFFFFF))
# Map the colors to the names
feature_annot['fill_colour'] = feature_annot['name'].map(colors)
#replace Na in fill color with default
feature_annot['fill_colour'] = feature_annot['fill_colour'].fillna('default')
#save tsv
#add back in .fna_final-1 to feature id before last underscore
feature_annot2=pd.read_csv('lovis_out_1/feature_annotation_table.tsv', sep='\t')
feature_annot['feature_id']=feature_annot2['feature_id']

feature_annot.to_csv('feature_annotation_table.tsv', sep='\t', index=False)

#Run lovis4u with the feature annotation table
lovis_2= 'lovis4u -gff '+gff_out+ ' --set-category-colour -c A4p2 -o lovis_out_2 -alip -hl -ufid -faf feature_annotation_table.tsv'
subprocess.call(lovis_2, shell=True)

#Delete lovis_out_1
delete= 'rm -r lovis_out_1'
subprocess.call(delete, shell=True)
#move GFF folders to output directory
move= 'mv lovis_out_2/lovis4u.pdf '+output_dir+'/Genome_map.pdf'
subprocess.call(move, shell=True)
#move annomazing files to output directory
move= 'mv '+output_dir+'_annomazing_final.csv '+output_dir
subprocess.call(move, shell=True)
#move gff files to output directory
move= 'mv '+gff_out+' '+output_dir
subprocess.call(move, shell=True)
#move prodigal files to output directory
move= 'mv '+output_dir+'_Prodigal '+output_dir
subprocess.call(move, shell=True)
#remove feature_annotation_table.tsv
delete= 'rm feature_annotation_table.tsv'
#subprocess.call(delete, shell=True)
#remove .err and .out files
delete= 'rm *.err'
subprocess.call(delete, shell=True)
delete= 'rm *.out'
subprocess.call(delete, shell=True)
#remove lovis_out_2 folder
delete= 'rm -r lovis_out_2'
subprocess.call(delete, shell=True)
#remove gff_files folder
delete= 'rm -r '+gff_folder
subprocess.call(delete, shell=True)
