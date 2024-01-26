import os
import pandas as pd
import sys
import glob
from collections import defaultdict
gencode = { 'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W' }
segments = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]

sample = sys.argv[1]

#Select allAlleles.txt files
file_list = glob.glob('A_*allAlleles.txt')

# exclude files that are empty by checking the variants.txt files
variant_files = glob.glob('A_*variants.txt')
notEmptyFiles = []
for var_file in variant_files:
   with open(var_file, 'r') as abc:
      lines = [i for i in abc.readlines()]
      if len(lines) > 1:
        file= var_file.split("-")[0] + "-allAlleles.txt"
        notEmptyFiles.append(file)
#create an empty list to store the name of the segments/files that encode extra proteins.
files_ExtraProt = []

for file in notEmptyFiles:
   if (segments[0]) in file or (segments[3]) in file or (segments[4]) in file or (segments[5]) in file:
      continue
   else:
      files_ExtraProt.append(file)

#read the file and specify the column names

def define_df(file_nv, pos0, pos1):
    cons_min_var=pd.read_table(file_nv, header=0)
    df = cons_min_var[["Reference_Name","Position","Allele", "Count", "Frequency", "Average_Quality", "Allele_Type"]]
    #filter the data
    df_filt_min = df.loc[(df["Average_Quality"] >= 30) & (df["Frequency"] > 0.02) & (df["Count"] > 15) &(df["Allele_Type"] == "Minority")] 
    df_filt_cons =df.loc[(df["Average_Quality"] > 20) & (df["Allele_Type"] == "Consensus")] 
    #find minority nt positions
    min_positions = df_filt_min["Position"].tolist()
    #create a list with all the alleles and merge them into a string IF the mutations are present between pos0 and pos1
    for pos in min_positions:
        if (pos0 <= pos <= pos1):
            cons = df_filt_cons["Allele"].tolist()
            #create dictionary to connect position and allele of the min variants.
            min_dic = dict (zip(df_filt_min['Position'] -1, df_filt_min['Allele']))
            cons_dic = dict (zip(df_filt_cons['Position'] -1, df_filt_cons['Allele']))
            #merge the 2 dictionaries
            all_dic = defaultdict(list)
            for d in (cons_dic, min_dic): # you can list as many input dicts as you want here
                for key, value in d.items():
                    all_dic[key].append(value)
            #get the minority consensus list
            minority = defaultdict(list)
            minority2 = []
            for key, values in all_dic.items():
                if len(values) == 1:
                    minority[key].append(values[0])
                    minority2.append(values[0])
                elif len(values) == 2:
                    minority[key].append(values[1])
                    minority2.append(values[1])
            consensus_maj = ''.join(cons)
            #make a list of starting position of each codon. 
            start_codons = []
            for i in range(0, len(consensus_maj), 3): 
                start_codons.append(i)
            #determine the consensus codon for the sequence
            cons_codons = []
            for i in start_codons:
                b = i+3
                x = consensus_maj[i:b]
                cons_codons.append(x)
            #determine the amino acid list and join them in a string.
            cons_aa = []
            for codon in cons_codons:
                for key, value in gencode.items():
                    if codon == key:
                        cons_aa.append(value)
            #get the minority consensus fasta
            minority_cons = ''.join(minority2)
            #Find the position of the codon for each mutations, find the amino acid and if the mutation is synonymous or nonsynonymous and write it in a file
            seg_name= file_nv.split("-")[0]
            output_file = sample + "_Mutations_" + file_nv.split("-")[0] + ".txt"
            with open(output_file, "w") as out:
                title = ["Sample", "Segment", "Position_real", "Position_py", "Major_Codon", "Major_aa", "Minor_Codon", "Minor_aa", "Mutation", '\n']
                out.write('\t'.join(title))
                for pos in start_codons:
                    x = start_codons.index(pos)
                    for minpos in min_positions:
                        pyminpos = minpos -1
                        if pos <= pyminpos & pyminpos < start_codons[x+1]:
                            for key, value in gencode.items():
                                maj_cod = consensus_maj[start_codons[x]:start_codons[x+1]]
                                min_cod = minority_cons[start_codons[x]:start_codons[x+1]]
                                if maj_cod == key:
                                    if gencode[min_cod] == value: 
                                        case1 = [str(sample), str(seg_name), str(minpos), str(pyminpos), maj_cod, value, min_cod, value, "Synonymous", '\n']
                                        out.write('\t'.join(case1))
                                    elif gencode[min_cod] != value:
                                        case2 = [str(sample), str(seg_name), str(minpos), str(pyminpos), maj_cod, value, min_cod, gencode[min_cod], "Non-Synonymous", '\n']
                                        out.write('\t'.join(case2))

def define_df2(file_yv, pos0, pos1, pos2, pos3):
    df = pd.read_table(file_yv, header=0)[["Reference_Name","Position","Allele", "Count", "Frequency", "Average_Quality", "Allele_Type"]]
    min_positions = df.loc[(df["Average_Quality"] >= 30) & (df["Frequency"] > 0.02) & (df["Count"] > 15) & (df["Allele_Type"] == "Minority")]["Position"].tolist()
    for pos in min_positions:
        if (pos0 <= pos <= pos1) | (pos2 <= pos <= pos3):
            cons = df.loc[(df["Average_Quality"] >= 20) & (df["Allele_Type"] == "Consensus")] ["Allele"].tolist()
            min_dic = dict (zip(df.loc[(df["Average_Quality"] >= 30) & (df["Frequency"] > 0.02) & (df["Count"] > 15) & (df["Allele_Type"] == "Minority")]['Position'] -1, df.loc[(df["Average_Quality"] >= 30) & (df["Count"] > 15) & (df["Frequency"] > 0.02) & (df["Allele_Type"] == "Minority")]['Allele']))
            cons_dic = dict (zip(df.loc[(df["Average_Quality"] > 20) & (df["Allele_Type"] == "Consensus")] ['Position'] -1, df.loc[(df["Average_Quality"] > 20) & (df["Allele_Type"] == "Consensus")] ['Allele']))
            all_dic = defaultdict(list)
            for d in (cons_dic, min_dic): 
                for key, value in d.items():
                    all_dic[key].append(value)
            minority = defaultdict(list)
            minority2 = []
            for key, values in all_dic.items():
                if len(values) == 1:
                    minority[key].append(values[0])
                    minority2.append(values[0])
                elif len(values) == 2:
                    minority[key].append(values[1])
                    minority2.append(values[1]) 
            
            if pos2 == pos3 == 0:
                consensus_maj_1 = ''.join(cons[pos0:pos1]) 
                minority_cons_1 = ''.join(minority2[pos0:pos1])  

            else:
                consensus_maj_1 = ''.join(cons[pos0:pos1]) + ''.join(cons[pos2:pos3])
                minority_cons_1 = ''.join(minority2[pos0:pos1]) + ''.join(minority2[pos2:pos3])
            
            start_codons_1 = []
            for i in range(0, len(consensus_maj_1), 3): 
                start_codons_1.append(i)
            cons_codons_1 = []
            for i in start_codons_1:
                cons_codons_1.append(consensus_maj_1[i:i+3])
            cons_aa_1 = []
            for codon in cons_codons_1:
                for key, value in gencode.items():
                    if codon == key:
                        cons_aa_1.append(value)

                
            
            seg_name= file_yv.split("-")[0] + "_2"
            output_file_1 = sample + "_Mutations_" + file_yv.split("-")[0] + "_prot2.txt"
            if pos2 == pos3 == 0:
                with open(output_file_1, "w") as out:
                    title = ["Sample", "Segment", "Position_real", "Position_py", "Major_Codon", "Major_aa", "Minor_Codon", "Minor_aa", "Mutation", '\n']
                    out.write('\t'.join(title))
                    for pos in start_codons_1:
                        x = start_codons_1.index(pos)
                        for minpos in min_positions:
                            pyminpos = minpos -1
                            if pos0 <= pyminpos <= pos1:
                                new_pyminpos = pyminpos - pos0
                                if pos <= new_pyminpos < start_codons_1[x+1]:
                                    for key, value in gencode.items():
                                        maj_cod = consensus_maj_1[start_codons_1[x]:start_codons_1[x+1]]
                                        min_cod = minority_cons_1[start_codons_1[x]:start_codons_1[x+1]]
                                        if maj_cod == key:
                                            if gencode[min_cod] == value: 
                                                case = [str(sample), str(seg_name), str(minpos), str(pyminpos), maj_cod, value, min_cod, value, "Synonymous", '\n']
                                                out.write('\t'.join(case))
                                            elif gencode[min_cod] != value:
                                                case = [str(sample), str(seg_name), str(minpos), str(pyminpos), maj_cod, value, min_cod, gencode[min_cod], "Non-Synonymous", '\n']
                                                out.write('\t'.join(case))
            else:
                with open(output_file_1, "w") as out:
                    title = ["Sample", "Segment", "Position_real", "Position_py", "Major_Codon", "Major_aa", "Minor_Codon", "Minor_aa", "Mutation", '\n']
                    out.write('\t'.join(title))
                    for pos in start_codons_1:
                        x = start_codons_1.index(pos)
                        for minpos in min_positions:
                            pyminpos = minpos -1
                            if pos0 <= pyminpos <= pos1:
                                new_pyminpos = pyminpos - pos0
                                if pos <= new_pyminpos < start_codons_1[x+1]:
                                    for key, value in gencode.items():
                                        maj_cod = consensus_maj_1[start_codons_1[x]:start_codons_1[x+1]]
                                        min_cod = minority_cons_1[start_codons_1[x]:start_codons_1[x+1]]
                                        if maj_cod == key:
                                            if gencode[min_cod] == value: 
                                                case1 = [str(sample),str(seg_name), str(minpos), str(new_pyminpos), maj_cod, value, min_cod, value, "Synonymous", '\n']
                                                out.write('\t'.join(case1))
                                            elif gencode[min_cod] != value:
                                                case2 = [str(sample),str(seg_name), str(minpos), str(new_pyminpos), maj_cod, value, min_cod, gencode[min_cod], "Non-Synonymous", '\n']
                                                out.write('\t'.join(case2))
                            elif pos2 <= pyminpos <= pos3:
                                new_pyminpos = pyminpos - pos2 + pos1
                                if pos <= new_pyminpos < start_codons_1[x+1]:
                                    for key, value in gencode.items():
                                        maj_cod = consensus_maj_1[start_codons_1[x]:start_codons_1[x+1]]
                                        min_cod = minority_cons_1[start_codons_1[x]:start_codons_1[x+1]]
                                        if maj_cod == key:
                                            if gencode[min_cod] == value: 
                                                case1 = [str(sample),str(seg_name), str(minpos), str(new_pyminpos), maj_cod, value, min_cod, value, "Synonymous", '\n']
                                                out.write('\t'.join(case1))
                                            elif gencode[min_cod] != value:
                                                case2 = [str(sample),str(seg_name), str(minpos), str(new_pyminpos), maj_cod, value, min_cod, gencode[min_cod], "Non-Synonymous", '\n']
                                                out.write('\t'.join(case2))

#create file for all the segments (1st prot)
for file in notEmptyFiles: 
    if (segments[0]) in file:
        define_df(file, 0, 2277)
    elif (segments[1]) in file:
        define_df(file, 0, 2274)
    elif (segments[2]) in file:
        define_df(file, 0, 2151)
    elif (segments[3]) in file:
        define_df(file, 0, 1701)
    elif (segments[4]) in file:
        define_df(file, 0, 1497)     
    elif (segments[5]) in file:
        define_df(file, 0, 1410)
    elif (segments[6]) in file:
        define_df(file, 0, 759)
    elif (segments[7]) in file:
        define_df(file, 0, 660)


# create second file for segments PB1, PA, MP, NS that have extra proteins 
for file in files_ExtraProt:
    if (segments[1]) in file:
        define_df2(file, 94, 367, 0, 0)
    elif (segments[2]) in file:
        define_df2(file, 0, 760, 0, 0)
    elif (segments[6]) in file:
        define_df2(file, 0, 26, 714, 982)     
    elif (segments[7]) in file:
        define_df2(file, 0, 30, 502, 838)


