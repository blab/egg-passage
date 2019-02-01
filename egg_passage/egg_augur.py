"""
Force all available h3n2 sequences that were passaged in eggs to be included by augur. Also include
any available unpassaged- or cell-passaged sequences from these strains (i.e. the virus was
sequenced before and after passaging in eggs)
"""

import argparse
import sys, os, subprocess, shutil
from Bio import SeqIO

def find_egg_seqs(fasta_file, summary_file):

    summary = open(summary_file, 'w+')

    egg_seqs = []
    force_include = []
    cell_matches = 0
    unpassaged_matches = 0
    both_matches = 0

    #Find all egg-passaged sequences
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        virus_name = str.split(seq_record.id, '|')[0]

        if 'egg' in virus_name:
            egg_seqs.append(virus_name)

    #Find all pairs for egg-passaged sequences and add to list of viruses to include during augur sampling
    for egg_seq in egg_seqs:
        seq_name = str.split(egg_seq,'-egg')[0]
        cell=0
        unpassaged=0

        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            if seq_name in seq_record.id:
                force_include.append(str.split(seq_record.id, '|')[0])

                if '-egg' not in seq_record.id:
                    if '-cell' in seq_record.id:
                        cell+=1
                    else:
                        unpassaged+=1

        cell_matches+=cell
        unpassaged_matches+=unpassaged
        if cell!=0 and unpassaged!=0:
            both_matches+=1

    summary.write('Number of egg-passaged sequences: %d \n' % len(egg_seqs))
    summary.write('Egg-passaged seqs with cell-passaged match: %d \n' % cell_matches)
    summary.write('Egg-passaged seqs with unpassaged match: %d \n' % unpassaged_matches)
    summary.write('Egg-passaged seqs with unpassaged AND cell-passaged match: %d \n' % both_matches)
    summary.close()

    return force_include

def run_augur(fasta_file, summary_file, resolution, titers):

    force_include = find_egg_seqs(fasta_file, summary_file)

    #Navigate to nextstrain
    os.chdir('../../nextstrain/augur/builds/flu')

    #Add force_include sequences to flu_info.py
    with open('flu_info.py', 'a') as f:
        f.write("reference_viruses['h3n2']=reference_viruses['h3n2']+" + str(force_include))

    #Run flu.prepare.py
    if titers != 'notiter':
        if titers == 'hi':
            titer_path = '../../../../egg-passage/egg_passage/input_data/h3n2_who_hi_concat_titers.tsv'
        if titers == 'fra':
            titer_path = '../../../../egg-passage/egg_passage/input_data/h3n2_who_fra_concat_titers.tsv'
        prep= subprocess.Popen(["python", "flu.prepare.py", "--sequences", "../../../fauna/data/h3n2_ha.fasta", "-r", str(resolution), "--titers", str(titer_path)], stdout=subprocess.PIPE)
        #print prepare.py output in terminal
        print(prep.stdout.read())
    else:
        prep= subprocess.Popen(["python", "flu.prepare.py", "--sequences", "../../../fauna/data/h3n2_ha.fasta", "-r", str(resolution)], stdout=subprocess.PIPE)
        print(prep.stdout.read())

    #Run flu.process.py
    proc= subprocess.Popen(["python", "flu.process.py", "--json", "prepared/flu_seasonal_h3n2_ha_"+str(resolution)+".json", "--titers_export", "--export_translations"], stdout=subprocess.PIPE)
    #print prepare.py output in terminal
    print(proc.stdout.read())

    #Copy augur output to egg_passage
    auspice_prefix = 'flu_seasonal_h3n2_ha_'
    if titers != 'notiter':
        #Change names of augur output files to contain titer tag
        if titers not in ['hi', 'fra']:
            titer_tag = 'titer'
        else:
            titer_tag = titers
    if titers == 'notiter':
        titer_tag = 'notiter'

    #Copy all necessary files
    for json_file in os.listdir("./auspice/"):
        if json_file.startswith(auspice_prefix+str(resolution)):
            os.rename('auspice/'+str(json_file), 'auspice/'+str(auspice_prefix+str(resolution)+'_'+str(titer_tag)+json_file.split(auspice_prefix+str(resolution))[1]))
    subprocess.call('cp auspice/flu_seasonal_h3n2_ha_'+str(resolution)+'_'+str(titer_tag)+'_* ../../../../egg-passage/egg_passage/augur/' + str(titer_tag), shell=True)


    #Remove force_include from flu_info
    with open("flu_info.py","r+") as f:
        new_f = f.readlines()
        f.seek(0)
        for line in new_f:
            if "reference_viruses['h3n2']=reference_viruses['h3n2']+" not in line:
                f.write(line)
        f.truncate()
    #Clean augur directories
    clean = subprocess.Popen(["sh", "clean_directory.sh"], stdout=subprocess.PIPE)
    print(clean.stdout.read())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Selectively include egg-passaged viruses during subsampling in augur flu build")
    parser.add_argument('-f', '--fasta_file', default = '../../nextstrain/fauna/data/h3n2_ha.fasta', help="filepath of h3n2 sequence fasta file")
    parser.add_argument('-s', '--summary_file', default = 'egg_summary_file.txt', help = "name of file to write summary data to")
    parser.add_argument('--resolution', choices=['2y', '3y', '6y', '12y'], default='6y', type = str,  help = "resolution for flu.prepare.py (default: 6y)")
    parser.add_argument('--titers', default= 'hi', type = str,  help = ".tsv file with titers data. Specify 'hi' or 'fra' to used combined egg and cell titers from WHO. Specify 'notiter' to run without titers. Otherwise specify filepath to titers data")
    args = parser.parse_args()
    run_augur(fasta_file = args.fasta_file, summary_file = args.summary_file, resolution = args.resolution, titers = args.titers)
