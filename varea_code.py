### THIS IS THE FULL VERSION OF THE TOOL:
### from patient's VARIANTS TSV to EVE ANNOTATED & RANKED tsv

# import packages:
import pandas as pd
import numpy as np
import sys
import re
import argparse, textwrap
import glob
import chardet


### ARGUMENTS
parser = argparse.ArgumentParser(description='VAREA - Variant Automatic Ranking and EVE Annotation\nVAREA is a Python command-line tool developed to facilitate routine variant analysis by performing the variant prioritisation step automatically.', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--input', '-i',
    dest='input',
    type=str,
    help="Patient(s) variants in a tsv format file. Or a txt file with filenames with paths to patients' tsv files if you wish to merge them before.",
    required=True)
parser.add_argument('--merge_files', '-merge',
    action='store_true',
    help='Flag if you provide a file with filenames of samples to be merged before annotation and ranking.\nBeware! Takes the first number (1-3 digits) from the file name (the whole path) as a sample number',
    required=False)
parser.add_argument('--gene_map', '-gm',
    dest='gmap',
    type=str,
    help='Gene names mapping file in the format: \ntab separated txt file consisting of a column with HGNC/GeneCards gene names and a column with matching UniProtKB IDs. \nSkips first row as a header. \n "/ngc/projects/gm_ext/valgof/mthesis/metadata/hgnc_to_uniprot_expanded.txt" ',
    required=True)
parser.add_argument('--filter_genes', '-fg',
    action='store_true',
    help='Flag if you wish to filter out variants outside the genes listed in a file provided with --gene_map/-gm.',
    required=False)
parser.add_argument('--iformat', '-if',
    dest='iformat',
    type=str,
    choices=['vcf', 'varseq'],
    default='varseq',
    help='Input format: vcf or varseq, default: varseq',
    required=False)
parser.add_argument('--varseq_af_column', '-af_col',
    dest='af_col',
    type=int,
    help='Index (1-based) of a column that contains gnom_AD allele frequency of the variants. Required if input file is from VarSeq.',
    required=False)
parser.add_argument('--af_bins', '-af_bins',
    dest='af_bins',
    type=str,
    nargs='+',
    default=[0,0.001,0.005, 0.01, 1],
    help='Specify AF binning type, default: [0,0.001,0.005,1] \nOther possibilites: \n1) a single number N to divide scores into N equal bins; \n2) a list, specifying binning ranges, rightmost included, e.g. -af_bins 0 0.1 0.5 1 produces 3 bins: (0,0.1]<(0.1,0.5]<(0.5,1]; \n3)"-af_bins quartile": divides AF into equal-sized bins based on sample quartiles.',
    required=False)
parser.add_argument('--add_score', '-score',
    dest='score',
    help='If you wish to add additional gene prioritisation score, provide a file in the format: \nfisrt column - gene names, second column - score; with a header. \nVariants will be prioritised by this score prior to Consequences and AF. Use flag -nas to avoid sorting on this score.',
    required=False)
parser.add_argument('--score_bins', '-s_bins',
    dest='score_bins',
    type=str,
    nargs='+',
    default=None,
    help='Specify binning type. \nDefault: None, uses existing score values as they are. \nOther possibilites: \n1) a single number N to divide scores into N equal bins; \n2) a list, specifying binning ranges, rightmost included, e.g. -s_bins 0 0.1 0.5 1 divides into 3 bins: (0,0.1]<(0.1,0.5]<(0.5,1];\n3) "-s_bins quartile": divides score into equal-sized bins based on sample quartiles.',
    required=False)
parser.add_argument('--ascending_score', '-a',
    action='store_true',
    help='Flag if optional gene prioritisation score should be sorted in ascending order. Default: descending.',
    required=False)
parser.add_argument('--no_add_sorting', '-nas',
    action='store_true',
    help='Flag if optional gene prioritisation score should not be used for additional sorting.',
    required=False)
parser.add_argument('--output', '-o',
    dest='out',
    type=str,
    help='Desired name for output file(s). Uses input name if argument is not given.',
    required=False)
parser.add_argument('--missense', '-m',
    action='store_true',
    help='Flag if you wish to produce a separate file for ranked missense variants.',
    required=False)

args = parser.parse_args()

### ______________________________________________________________________
### SETUP, ERRORS AND MESSAGES
# to disable false warnings
pd.options.mode.chained_assignment = None
# print out help when no arguments are passed
if len(sys.argv) == 1: parser.print_help(); sys.exit(1)
# error message if column index for AF is not given for VarSeq input file
if (args.iformat == 'varseq') & (args.af_col is None): print('Please specify the index (1-based) of a column that contains variant AF with -af_col option.'); exit()
# print warning that Optional gene prioritisation Score is sorted in descending order, if the user does not specify the order; just in case
optional_score = args.score != None
if optional_score & (args.ascending_score == False) & (args.no_add_sorting == False): print('Please, note: variants will be sorted using the provided additional gene prioritisation score in DESCENDING ORDER.\nIf you wish to change this, add a flag: -a (--ascending_score)')
if (optional_score == False) & ((args.ascending_score == True) | (args.no_add_sorting == True)): print('Additional gene prioritisation score was not provided (use: -score(--add_score)). Proceeding to ranking without it.')
if (args.no_add_sorting == True) & ((args.ascending_score == True) | (args.score_bins != None)) : print('Warning: additional gene prioritisation score is not used for sorting (remove -nas(--no_add_sorting) flag if you wish to change this).')


### ______________________________________________________________________
print('Processing input files...')
### READ IN VCF/VarSeq
if args.merge_files == True:
    samplefiles = []
    with open(args.input, "r") as f:
        samplefiles = f.read().splitlines()
    allsamples = []
    for samplename in samplefiles:
        # Check encoding
        with open(samplename, 'rb') as f:
            res = chardet.detect(f.read())
        sample_i = pd.read_csv(samplename, sep='\t', header=0, dtype=str, na_values='.', encoding=res['encoding'])
        # get sample number from the file name (1,2 or 3 digits)
        match = re.search('[\d]{1,3}', samplename)
        sample_i['Sample'] = str(match.group())
        allsamples.append(sample_i)
    df_sample = pd.concat(allsamples)
    print('files successfully merged')
else:
    # Check encoding
    with open(args.input, 'rb') as f:
        res = chardet.detect(f.read())
        # 'ISO-8859-1'
    df_sample = pd.read_csv(args.input, sep='\t', header=0, dtype=str, na_values='.', encoding=res['encoding'])
### ______________________________________________________________________


# assigning conditions/some values to variables for further convenience
if args.iformat == 'vcf':
    conseq_col = 'Consequence'
    af_col = 'gnomAD_AF'
    gene_col = 'SYMBOL'
    hgvsp = 'HGVSp'
elif args.iformat == 'varseq':
    af_col = df_sample.columns[args.af_col-1]
    for c in df_sample.columns:
        match_conseq = re.search('Sequence[ _]Ontology[ _]\(Clinically[ _]Relevant\)', c)
        if match_conseq: conseq_col = match_conseq.group()
        match_gene = re.search('Gene[ _]Names', c)
        if match_gene: gene_col = match_gene.group()
        match_hgvsp = re.search('HGVS[ _]p.[ _]\(Clinically[ _]Relevant\)', c)
        if match_hgvsp: hgvsp = match_hgvsp.group()
print('column names retrieved')
### OPTIONAL
### Filter out genes not present in the list
if args.filter_genes == True:
    gene_map = pd.read_csv(args.gmap, sep='\t',header=0, dtype=str)
    df_sample = df_sample.loc[df_sample[gene_col].isin(gene_map.iloc[:,0])]
    print('genes filtered')

### ______________________________________________________________________
### Pipeline starts with annotating missense variants with EVE scores.
### 0. Subselect missense variants only
df_sample_missense=df_sample[df_sample[conseq_col].str.contains('missense_variant')]


### ______________________________________________________________________
### GETTING EVE SCORES FOR MISSENSE
print('Annotating missense variants with EVE scores...')
### ______________________________________________________________________
### 1. Converting HGNC into UniProtKB ID in the sample df using provided mapping file
def convert_gene_names(df, how='h2u'):
    """ Convert gene names in the sample file
    if how == 'h2u': from HGNC to UniProt_ID
    if how == 'u2h': from UniProt_ID to HGNC
    format == 'vcf'  for tsv originating from VCF (format = args.iformat when using a function)
    format == 'varseq' for tsv originating from VarSeq
    ________________________________________________
    ?: ask to input a mapping file OR where to keep the mapping file for all genes
    ?: check why there're 1957 genes mapped to 1965 UniProt id's
    """
    # reading mapping file as a dictionary
    hgnc, uniprot = np.genfromtxt(args.gmap, delimiter='\t', unpack=True, dtype=str, skip_header=1)
    map_dict = {}
    if how == 'h2u':
        for hgnc, uniprot in zip(hgnc, uniprot):
            map_dict[hgnc] = uniprot
    elif how == 'u2h':
        for hgnc, uniprot in zip(hgnc, uniprot):
            map_dict[uniprot] = hgnc
    # adding a UniProtKB ID column into df_sample
    df['UniProt_ID'] = df[gene_col].map(map_dict).fillna('.')
    print('UniProt IDs retrieved. Getting EVE scores...')

convert_gene_names(df_sample_missense, how='h2u')

### ______________________________________________________________________
### 2. Extracting list of genes from sample df that have EVE score
# 1) get unique gene IDs
eve_gene_names = df_sample_missense['UniProt_ID'].unique()
# 2) format list so it becomes a filenames list
eve_gene_files = ['/ngc/shared/EVE/variant_files/' + x + '.csv' for x in eve_gene_names]
# 3) check if genes are in EVE folder
# 3.1) load all file names from eve folder
evefolder_files = glob.glob('/ngc/shared/EVE/variant_files/*.csv')
# 3.2) check presence and drop gene from list if not present in eve folder
for g in eve_gene_files[:]: # using a copy of list here so indexes don't shift when removing elements
    if g in evefolder_files:
        continue
    else:
        eve_gene_files.remove(g)
# check that this list is not empty
if bool(eve_gene_files) == False:
    #df_sample['Substitution_EVE','Protein_pos_EVE','EVE_score','EVE_uncertainty','EVE_class_75']="."
    #df_sample.to_csv(args.out, sep='\t', index=False, na_rep='.')
    print('All variants are in genes with no EVE score. Missense variants will be prioritized by Allele Frequency only')

### ______________________________________________________________________
### 3. Concat all required EVE files
print('Adding EVE scores...')
temp_filename_list = []
for f in eve_gene_files:
    # read file
    df_eve_temp = pd.read_csv(f, header=0, usecols=['position', 'wt_aa', 'mt_aa', 'EVE_scores_ASM',
    'uncertainty_ASM', 'EVE_classes_75_pct_retained_ASM'],dtype=str)
    # get UniProt ID from filename
    df_eve_temp['UniProt_ID'] = f[30:-4]
    temp_filename_list.append(df_eve_temp)
df_eve = pd.concat(temp_filename_list, axis=0, ignore_index=True)

### Add Substitutuion field in format Aaa<position>Aaa
aa_dict_1to3 = {'A': 'Ala', 'C': 'Cys', 'D': 'Asp', 'E': 'Glu', 'F': 'Phe', 'G': 'Gly',
'H': 'His', 'I': 'Ile', 'K': 'Lys', 'L': 'Leu', 'M': 'Met', 'N': 'Asn', 'P': 'Pro',
'Q': 'Gln', 'R': 'Arg', 'S': 'Ser', 'T': 'Thr', 'V': 'Val', 'W': 'Trp', 'Y': 'Tyr'}
# convert A to Aaa format in EVE df:
for col in ['wt_aa','mt_aa']:
    df_eve[col] = df_eve[col].map(aa_dict_1to3).fillna('.')
# adding custom 'substitution' field - to later match with df_sample substitution field
df_eve.insert(0, 'Substitution_EVE', df_eve.wt_aa+df_eve.position+df_eve.mt_aa)
# dropping unnecessary mt and wt
df_eve = df_eve.drop(['wt_aa', 'mt_aa', 'position'], axis=1)


### ______________________________________________________________________
### ADDING EVE DATA TO MISSENSE SAMPLE DF
### ______________________________________________________________________
# 1. Getting 'Substitution' column for df_sample - from HGVSp in the format: Aaa<position>Aaa
df_sample_missense['Substitution_sample'] = df_sample_missense[hgvsp]
df_sample_missense.Substitution_sample = df_sample_missense.Substitution_sample.str.replace(r'.*:p.', r'')

### Merging based on Substitution: 'Aaa<position>Aaa' and Protein names: 'UniProt_ID'
df_missense_merged = df_sample_missense.merge(df_eve, how='left', left_on=['Substitution_sample', 'UniProt_ID'], right_on=['Substitution_EVE', 'UniProt_ID'])
df_missense_merged.rename(columns={'EVE_scores_ASM': 'EVE_score', 'uncertainty_ASM': 'EVE_uncertainty', 'EVE_classes_75_pct_retained_ASM': 'EVE_class_75'}, inplace=True)
df_missense_merged = df_missense_merged.drop('Substitution_sample', axis=1)


### ______________________________________________________________________
### RANKING BASED ON CONSEQUENCES, EVE AND AF
### ______________________________________________________________________
### 1. subselect rows with all the rest variants (not missense)
### 2. merge back missense variants to the rest so that df_all also has EVE columns
### 3. sort all on CONSEQS, EVE and AF
### 4. OPTIONAL (for ibd pipeline): annotate with gene annotations counts

variant_severity = {'transcript_ablation[&,]?.*': '57','splice_acceptor_variant[&,]?.*': '56','splice_donor_variant[&,]?.*': '55','stop_gained[&,]?.*': '54','frameshift_variant[&,]?.*': '53','stop_lost[&,]?.*': '52','start_lost[&,]?.*': '51','transcript_amplification[&,]?.*': '50',
'inframe_insertion[&,]?.*': '44','inframe_deletion[&,]?.*': '43','missense_variant[&,]?.*': '42','protein_altering_variant[&,]?.*': '41', 'regulatory_region_ablation[&,]?.*':'40',
'splice_region_variant[&,]?.*': '37','splice_donor_5th_base_variant[&,]?.*': '36','splice_donor_region_variant[&,]?.*': '35','splice_polypyrimidine_tract_variant[&,]?.*': '34','incomplete_terminal_codon_variant[&,]?.*': '33','start_retained_variant[&,]?.*': '32','stop_retained_variant[&,]?.*': '31','synonymous_variant[&,]?.*': '30',
'coding_sequence_variant[&,]?.*': '27', 'mature_miRNA_variant[&,]?.*': '26','5_prime_UTR_variant[&,]?.*': '25','3_prime_UTR_variant[&,]?.*': '24','non_coding_transcript_exon_variant[&,]?.*': '23','non_coding_exon_variant[&,]?.*': '23','intron_variant[&,]?.*': '22','NMD_transcript_variant[&,]?.*': '21','non_coding_transcript_variant[&,]?.*': '20','upstream_gene_variant[&,]?.*': '19','downstream_gene_variant[&,]?.*': '18','TFBS_ablation[&,]?.*': '17','TFBS_amplification[&,]?.*': '16','TF_binding_site_variant[&,]?.*': '15','regulatory_region_amplification[&,]?.*': '14','feature_elongation[&,]?.*': '13','regulatory_region_variant[&,]?.*': '12','feature_truncation[&,]?.*': '11','intergenic_variant[&,]?.*': '10'}


# joining all back
df_sample_rest=df_sample[df_sample[conseq_col].str.contains('missense_variant')==False]
#df_sample_rest.reset_index(drop=True, inplace=True)
# if occasionally initial df had already EVE annotations, drop duplicates:
# df_missense_merged.drop(columns=df_missense_merged.columns[df_missense_merged.columns.duplicated(keep='last')], inplace=True)
df_all_ranked = df_missense_merged.append(df_sample_rest, ignore_index=True)
df_all_ranked.replace('.', np.NaN)
print('Missense variants successfully annotated with EVE scores.\n Ranking in process...')

#_____________________________________________________________________
### Check that AF col has only 1 value per entry. Otherwise take the highest and print out warning
counts=0
afs = []
for n in df_all_ranked[af_col].astype('string').fillna('.'):
    if ',' in n:
        counts+=1
        afs.append(n)
mask = df_all_ranked[af_col].isin(afs)
df_rows= df_all_ranked[mask]
df_all_ranked.loc[mask, af_col] = df_rows[af_col].str.split(',').apply(lambda x: max(x))
print(counts, "variants have multiple allele frequencies given. The ranking of these variants will be based on the highest one, the other omitted.\n")
#_____________________________________________________________________

# AF binning
if args.af_bins[0] == "quartile":
    # print('af qartile')
    df_all_ranked['AF_bins']=pd.qcut(df_all_ranked[af_col].astype('float', errors='ignore'), 4)
elif len(args.af_bins) != 1: # a custom list, rightmost edge included, incl. default
    # print('af custom bins')
    df_all_ranked['AF_bins']=pd.cut(df_all_ranked[af_col].astype('float', errors='ignore'), bins=list(np.array(args.af_bins, dtype=float)))
else: # n of equal bins
    # print('af n equal bins')
    df_all_ranked['AF_bins']=pd.cut(df_all_ranked[af_col].astype('float', errors='ignore'), bins=int(args.af_bins[0]))


# optional score
if optional_score:
    df_ann = pd.read_csv(args.score, sep='\t',header=0, dtype=str)
    df_all_ranked=df_all_ranked.merge(df_ann,how='left',left_on=[gene_col], right_on=[df_ann.columns[0]])
    df_all_ranked.drop(columns=df_ann.columns[0],inplace=True)
    if args.no_add_sorting == True: # no sorting based on additional gene prioritisation score
        print('no_sorting=true')
        df_all_ranked.sort_values(by=[conseq_col, 'AF_bins', 'EVE_score'], key=lambda x: pd.to_numeric(x.replace(variant_severity, regex=True), errors='ignore'), ascending=(False,True,False), inplace=True)
    elif args.score_bins == None: # if binning type not specified - use score values to sort as they are, without binning
        print('score_bins=none')
        df_all_ranked.sort_values(by=[df_ann.columns[1], conseq_col, 'AF_bins', 'EVE_score'], key=lambda x: pd.to_numeric(x.replace(variant_severity, regex=True), errors='ignore'), ascending=(args.ascending_score,False,True,False), inplace=True)
    elif args.score_bins[0] == "quartile":
        print('quartile')
        df_all_ranked['Additional_score_bins']=pd.qcut(df_all_ranked[df_ann.columns[1]].astype('float', errors='ignore'), 4)
        df_all_ranked.sort_values(by=['Additional_score_bins', conseq_col, 'AF_bins', 'EVE_score'], key=lambda x: pd.to_numeric(x.replace(variant_severity, regex=True), errors='ignore'), ascending=(args.ascending_score,False,True,False), inplace=True)
    elif len(args.score_bins) != 1: # a custom list, rightmost edge included
        print('custom bins')
        df_all_ranked['Additional_score_bins']=pd.cut(df_all_ranked[df_ann.columns[1]].astype('float', errors='ignore'), bins=list(np.array(args.score_bins, dtype=float)))
        df_all_ranked.sort_values(by=['Additional_score_bins', conseq_col, 'AF_bins', 'EVE_score'], key=lambda x: pd.to_numeric(x.replace(variant_severity, regex=True), errors='ignore'), ascending=(args.ascending_score,False,True,False), inplace=True)
    else: # n of equal bins
        print('n equal bins')
        df_all_ranked['Additional_score_bins']=pd.cut(df_all_ranked[df_ann.columns[1]].astype('float', errors='ignore'), bins=int(args.score_bins[0]))
        df_all_ranked.sort_values(by=['Additional_score_bins', conseq_col, 'AF_bins', 'EVE_score'], key=lambda x: pd.to_numeric(x.replace(variant_severity, regex=True), errors='ignore'), ascending=(args.ascending_score,False,True,False), inplace=True)
else:
    # sort without arranging AF into bins:
    # df_all_ranked.sort_values(by=[conseq_col, af_col, 'EVE_score'], key=lambda x: pd.to_numeric(x.replace(variant_severity, regex=True), errors='ignore'), ascending=(False,True,False))
    df_all_ranked.sort_values(by=[conseq_col, 'AF_bins', 'EVE_score'], key=lambda x: pd.to_numeric(x.replace(variant_severity, regex=True), errors='ignore'), ascending=(False,True,False), inplace=True)

#filter out canonical (already done usually)
# df_all_ranked=df_all_ranked.loc[df_all_ranked.CANONICAL=='YES']
print('Variants ranked. Saving output...')

### saving to output
if args.out: df_all_ranked.to_csv(args.out+'.tsv', sep='\t', index=False, na_rep='.')
else: df_all_ranked.to_csv(args.input[:-4]+'_ranked.tsv', sep='\t', index=False, na_rep='.')

### saving missense separately
if args.missense == True:
    df_final_missense=df_all_ranked[df_all_ranked[conseq_col].str.contains('missense_variant')]
    if args.out:
        df_final_missense.to_csv(args.out+'_missense.tsv', sep='\t', index=False, na_rep='.')
    else:
        df_final_missense.to_csv(args.input[:-4]+'_ranked_missense.tsv', sep='\t', index=False, na_rep='.')

print('Done.')
