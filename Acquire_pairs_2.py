
# coding: utf-8

# In[1]:



# <nbformat>3.0</nbformat>

# <markdowncell>
# <codecell>



import pandas as pd
import numpy as np
import os 
from collections import Counter
from scipy import spatial
import sys
import networkx as nx
sys.setrecursionlimit(10000)


def str_to_bool(s):
    if s == 'True':
         return True
    elif s == 'False':
         return False
    else:
         raise ValueError # evil ValueError that doesn't tell you what the wrong value was



###############################################################################
###
###   Compute random genes matching for distance and chromosome
###
###############################################################################


def negative_eQTL_oneSNP(eQTL, distance_threshold, egene_list):
    '''
    match the eQTL using chr and distance
    '''
    eSNP_pos = eQTL['SNP_POS']
    echr = eQTL['SNP_CHR']
    egene = eQTL['GENE']

    t = GRCh37_genes.apply(lambda x: abs(x['gene_start'] - eQTL['GENE_START_POS']) < distance_threshold, axis=1)
    neg_gene_pool = GRCh37_genes[t][GRCh37_genes[t].apply(lambda x: x['GENE'] not in egene_list,axis=1)]
    
    if sum(neg_gene_pool['expression_mean'] > 0 ) > 0:
        neg_gene_pool = neg_gene_pool[neg_gene_pool['expression_mean'] > 0]
    if sum(neg_gene_pool['gene_description'] == 'protein_coding') > 0:
        neg_gene_pool = neg_gene_pool[neg_gene_pool['gene_description'] == 'protein_coding']
    
    if len(neg_gene_pool) > 0 :
        neg_gene = neg_gene_pool.sample(1)
    else:
        print "No matched random gene for %s" % eQTL['SNP']

    return list(neg_gene['GENE'])[0]
    




def negative_eQTL_onecell(eQTL_list, distance_threshold = 1000000):

    print 'Compute negative SNP-gene pairs match for distance and chromosome'
    
    eQTL_list = eQTL_list.loc[eQTL_list.groupby('GENE')['P-VALUE'].idxmin()]
    eQTL_list = eQTL_list.reset_index(drop=True)
    negative_genes = eQTL_list.apply(lambda x: negative_eQTL_oneSNP(x,distance_threshold, list(eQTL_list['GENE'])),axis=1)
    negative_genes = GRCh37_genes.loc[np.array(negative_genes)].reset_index(drop=True)

    eQTL_nega_pergene = eQTL_list.merge(negative_genes, left_index=True, right_index=True)
    eQTL_nega_pergene = eQTL_nega_pergene[['SNP', 'SNP_CHR', 'SNP_POS', 'GENE_y', 'gene_start','expression_mean', 'gene_description']]
    eQTL_nega_pergene.columns = [['SNP', 'SNP_CHR', 'SNP_POS', 'GENE', 'GENE_START_POS','expression_mean', 'gene_description']]
    eQTL_nega_pergene['DISTANCE'] = abs(eQTL_nega_pergene['SNP_POS'] - eQTL_nega_pergene['GENE_START_POS'])

    return eQTL_nega_pergene



###############################################################################
###
###   Retrive information from PC-HiC data
###
###############################################################################



######  Compute degree of contacts (0/1/2) from genes to the SNPs


def map_to_1based(anarray):
    return dict(zip(set(anarray), xrange(len(anarray))))


def bait_names(mapped_1based_idx, names):
    return dict(zip(mapped_1based_idx, names))

def edges_inbetween(D, n_nodes):
    ### input: a graph object
    ### return the degree of contacts between nodes: 1/2/0
    goal_matrix = np.ones([n_nodes,n_nodes]) * 100
    for n,nbrs in D.adjacency_iter():
        # assign 1 to the direct linkage
        goal_matrix[n][D[n].keys()] = 1
        # assign 2 to the indirect linkage
        for t in D[n].keys():
            goal_matrix[n][D[t].keys()] = map(lambda x: min(x), 
                                              zip(list(goal_matrix[n][D[t].keys()]),
                                                  [goal_matrix[n][t]+1] * len(goal_matrix[n][D[t].keys()])))
    goal_matrix[goal_matrix==100] = 0
    np.fill_diagonal(goal_matrix,0)
    return goal_matrix



def contacting_to_degree(data,edge_threshold=0):
    ## obtain the nodes and the weights edges
    if edge_threshold != 0:
        data = data.iloc[np.where(data['obs_exp']>edge_threshold)[0]]
    fragment_len = max(max(data['i']),max(data['j'])) + 1
    nodes = range(fragment_len)
    weighted_edges = zip(list(data['i']),list(data['j']),list(data['obs_exp']))
    print fragment_len
    
    ## construct the graph
    G = nx.Graph()
    G.add_nodes_from(nodes)
    G.add_weighted_edges_from(weighted_edges)
    
    ## obtain number of edges between the bins
    number_of_edges_inbetween = edges_inbetween(G, fragment_len)
    return G, number_of_edges_inbetween



def compute_frag_by_frag(d, celltype):
    anarray = list(d['baitID']) + list(d['oeID'])
    t = map_to_1based(anarray)
    d['baitID'] = map(lambda x: t[x], list(d['baitID']))
    d['oeID'] = map(lambda x: t[x], list(d['oeID']))
    baitsnames = bait_names(list(d['baitID']) + list(d['oeID']), list(d['baitName']) + list(d['oeName']))

    df = d[['baitID','oeID',celltype]]
    df.columns = ['i','j','obs_exp']

    df = df.reset_index(drop=True)

    G,frag_by_frag = contacting_to_degree(df,edge_threshold=0)
    frag_by_frag = pd.DataFrame(frag_by_frag)
    frag_by_frag.index = baitsnames.values()
    frag_by_frag.loc['NOGENE'] = [0]*len(frag_by_frag)
    
    frag_by_frag = frag_by_frag.drop(["."])
    frag_by_frag = frag_by_frag.to_dict(orient = 'index')
    
    return frag_by_frag


def find_the_key(frag_by_frag, gg):
    '''
    To deal with the baits' names, which has several gene names in one string
    '''
    gg_in_frag = np.where(map(lambda x: gg in str(x), frag_by_frag.keys()))[0]
    if len(gg_in_frag) > 0:
        temp = np.array(frag_by_frag.keys())[np.array(gg_in_frag)]
        return np.random.choice(temp)
    else:
        return 'NOGENE' 


def degree_of_contact_SNP_gene(eQTL_df, frag_by_frag, PIR_window=50000):

    if compute_chromosome_pos_LD:
        SNPs_PIR = chromosome_pos_LD[list(eQTL_df['SNP_POS'].astype('int') - 1)]
    else:
        ### find the other fragments where eSNP locate in 
        left = np.array(pairs['oeStart']) - PIR_window
        right = np.array(pairs['oeEnd']) + PIR_window
        SNPs_PIR = eQTL_df.apply(lambda x: [x['SNP_POS'] in range(t[0],t[1]) for t in zip(left,right)], axis=1)
        SNPs_PIR = SNPs_PIR.apply(lambda x: list(np.where(x)[0]))   ## a list of lists

    eQTL_df['SNP_contacting_df_idx'] = list(SNPs_PIR)

    ### find the degree of contact by using frag_by_frag dictionary
    degree_contact = eQTL_df.apply(lambda eq: [frag_by_frag[find_the_key(frag_by_frag, eq['GENE'])][pairs.iloc[t]['oeID']] 
                                               for t in eq['SNP_contacting_df_idx']], axis=1)
    ## print statistic 
    temp = map(lambda x:len(x),list(SNPs_PIR))
    print 'There are %i SNPs among %i that fall in at least one fragment (50kb window)' % (sum(np.array(temp) >0),len(temp))

    return degree_contact






######   LD blocks that overlap fragments

def overlapping_part(x1,x2,y1,y2):
    common_part = min(x2,y2) - max(x1,y1)
    if common_part > 0:
        return common_part
    else:
        return 0

def LD_in_fragments_one_SNP(eQTL, gene_th=500000):
    '''
    Allow a window for the genes
    '''
    temp = pairs.apply(lambda x: overlapping_part(x['oeStart'],x['oeEnd'],eQTL['LD_START'],eQTL['LD_END']),axis=1)
    if sum(temp) == 0:
        return 0
    else:
        ld_in_frag = pairs.iloc[list(temp).index(max(list(temp)))]
        if (np.abs(ld_in_frag['baitStart'] - eQTL['GENE_START_POS']) < gene_th ) or (np.abs(ld_in_frag['baitEnd'] - eQTL['GENE_START_POS']) < gene_th ):
            return 1
        else:
            return 0






######   Annotate the SNP-gene pairs


def eQTL_in_fragments_one_cell(eQTL_list):
    print 'Compute whether the eQTL pair fall in any promoter-PIR pair'
    
    ## approach 1: compute the degree of contacts
    frag2 = compute_frag_by_frag(pairs,celltype) 
    degree_contact = degree_of_contact_SNP_gene(eQTL_list, frag2)

    contact_vector = np.zeros(len(eQTL_list))
    contact_vector[np.where(map(lambda x: 2 in x, degree_contact))[0]] = 2
    contact_vector[np.where(map(lambda x: 1 in x, degree_contact))[0]] = 1
    eQTL_list['contacting2'] = contact_vector
    print Counter(contact_vector)
        
    ## approach 2: use the LD block to represent SNP
    contact_vector = eQTL_list.apply(lambda x: LD_in_fragments_one_SNP(x),axis=1)
    eQTL_list['contacting'] = contact_vector

    print Counter(contact_vector)
                
    return eQTL_list





######   Save the bed files for intersection with TFB motif

def save_bed_files(result_df, SNP_window, promoter_window, negative_set = False):

    SNP_bed = []
    SNP_bed.append(['chr%i'% x for x in (result_df['SNP_CHR'])])
    SNP_bed.append(list(np.array(result_df['SNP_POS']) - SNP_window))
    SNP_bed.append(list(np.array(result_df['SNP_POS']) + SNP_window))
    SNP_bed.append(list(result_df['SNP']))
    SNP_bed = pd.DataFrame(SNP_bed).transpose()
    if negative_set:
        SNP_bed.to_csv(os.path.join(DATA_DIR, 'intermediate/%s/PCHiC_peak_matrix_cutoff5_%i_eSNP_negativeset.bed'% (celltype, chr)),
                   sep='\t',index=False, header=False)
    else:
        SNP_bed.to_csv(os.path.join(DATA_DIR, 'intermediate/%s/PCHiC_peak_matrix_cutoff5_%i_eSNP.bed' % (celltype, chr)),
                   sep='\t',index=False, header=False)
    
    ## 
    promoter_bed = []
    promoter_bed.append(['chr%i'% x for x in (result_df['SNP_CHR'])])
    promoter_bed.append(list(np.array(result_df['GENE_START_POS']) - promoter_window))
    promoter_bed.append(list(np.array(result_df['GENE_START_POS']) + promoter_window))
    promoter_bed.append(list(result_df['GENE']))
    promoter_bed = pd.DataFrame(promoter_bed).transpose()
    if negative_set:
        promoter_bed.to_csv(os.path.join(DATA_DIR, 'intermediate/%s/PCHiC_peak_matrix_cutoff5_%i_egene_negativeset.bed'% (celltype, chr)),
                   sep='\t',index=False, header=False)
    else:
        promoter_bed.to_csv(os.path.join(DATA_DIR, 'intermediate/%s/PCHiC_peak_matrix_cutoff5_%i_egene.bed'% (celltype, chr)),
                   sep='\t',index=False, header=False)
        
        


###############################################################################
###
###   Scan the bed files for known motifs using fimo (in bash)
###
###############################################################################

#### run the following codes in DATA_DIR
#### bash TFBS.sh ${celltype} eSNP
#### bash TFBS.sh ${celltype} gene



###############################################################################
###
###   Annotate the pairs with TF motif information
###
###############################################################################


def merge_TF_motif(result_df, negative_set=False):
    if negative_set:
        f1 = 'fimo.output.%s.eSNP.%i.negativeset.bed' % (celltype,int(chr))
        f2 = 'fimo.output.%s.egene.%i.negativeset.bed' % (celltype,int(chr))
    else:
        f1 = 'fimo.output.%s.eSNP.%i.bed' % (celltype,int(chr))
        f2 = 'fimo.output.%s.egene.%i.bed' % (celltype,int(chr))
        
    # read in eQTLs with TF motif annotated
    eSNP_tf = pd.read_csv(os.path.join(DATA_DIR,'intermediate/%s/motif' % celltype,f1),sep='\t',header=None)
    egene_tf = pd.read_csv(os.path.join(DATA_DIR,'intermediate/%s/motif' % celltype,f2),sep='\t',header=None)
    
    eSNP_tf.columns = ['SNPID','chr','start','end','motif_SNP']
    egene_tf.columns = ['geneID','chr','start','end','motif_gene']
    
    eSNP_tf = pd.DataFrame(eSNP_tf.groupby('SNPID')['motif_SNP'].apply(lambda x:list(set(x))))
    egene_tf = pd.DataFrame(egene_tf.groupby('geneID')['motif_gene'].apply(lambda x:list(set(x))))
    
    result_df = result_df.merge(eSNP_tf, left_on='SNP',right_index=True, how='left')
    result_df = result_df.merge(egene_tf, left_on='GENE',right_index=True, how='left')

    return result_df






def merge_all_TF_motif(result_df, negative_set=False):
    if negative_set:
        f1 = '%s.eSNP.%i.negativeset.bed' % (celltype,int(chr))
        f2 = '%s.egene.%i.negativeset.bed' % (celltype,int(chr))
    else:
        f1 = '%s.eSNP.%i.bed' % (celltype,int(chr))
        f2 = '%s.egene.%i.bed' % (celltype,int(chr))
        
    # read in eQTLs with TF motif annotated
    eSNP_tf = pd.read_csv(os.path.join(DATA_DIR,'intermediate/%s/motif' % celltype,f1),sep='\t',header=None)
    egene_tf = pd.read_csv(os.path.join(DATA_DIR,'intermediate/%s/motif' % celltype,f2),sep='\t',header=None)
    
    eSNP_tf.columns = ['SNPID','chr','start','end','motif_SNP_all']
    egene_tf.columns = ['geneID','chr','start','end','motif_gene_all']
    
    eSNP_tf = pd.DataFrame(eSNP_tf.groupby('SNPID')['motif_SNP_all'].apply(lambda x:list(set(x))))
    egene_tf = pd.DataFrame(egene_tf.groupby('geneID')['motif_gene_all'].apply(lambda x:list(set(x))))
    
    result_df = result_df.merge(eSNP_tf, left_on='SNP',right_index=True, how='left')
    result_df = result_df.merge(egene_tf, left_on='GENE',right_index=True, how='left')

    return result_df




###############################################################################
###
###   Retrive information from ATAC-seq data
###
###############################################################################



def readin_ATACseq(filename):
    ATACseq = pd.read_csv(filename,sep='\t', header=None)
    ATACseq.columns = ['chr','Start','End','Peak','count']
    return ATACseq



def eQTL_in_ATACseq_one_chr(eQTL_pos, ATAC_regions, ATAC_window = 1000):

    N = len(ATAC_regions)
    eQTL_tree = spatial.KDTree(eQTL_pos)

    frag_start = np.reshape(ATAC_regions['Start'],[N,1])
    eQTL_near_fragStart = eQTL_tree.query_ball_point(frag_start, ATAC_window)
    frag_end = np.reshape(ATAC_regions['End'],[N,1])
    eQTL_near_fragEnd = eQTL_tree.query_ball_point(frag_end, ATAC_window)

    eQTL_near_frag = [list(set(eQTL_near_fragStart[i] + eQTL_near_fragEnd[i])) for i in xrange(N)]
    eQTLID = [x for x in eQTL_near_frag if len(x) > 0]
    eQTLID = list(set([a for b in eQTLID for a in b]))
    
    x = np.zeros(len(eQTL_pos))
    x[eQTLID] = 1
    return x



def eQTL_in_ATACseq_one_cell(eQTL, ATACseq, SNP_ATAC_window = 1000, gene_ATAC_window = 2000):
    
    ATAC = ATACseq.copy()
    ATAC = ATAC[ATAC['chr'] == 'chr%i'%chr].reset_index(drop=True)
    
    if "ATAC_SNP" in eQTL.columns:
        del eQTL['ATAC_SNP']
        del eQTL['ATAC_gene']
    if 1:
        result_df = pd.DataFrame()
        for name, g in eQTL.groupby('SNP_CHR'):
            ## eSNP
            eSNP_pos = np.reshape(g['SNP_POS'],[len(g),1])
            ATAC_SNP_vector = eQTL_in_ATACseq_one_chr(eSNP_pos,ATAC, SNP_ATAC_window)
            ## egene
            egene_pos = np.reshape(g['GENE_START_POS'],[len(g),1])
            ATAC_gene_vector = eQTL_in_ATACseq_one_chr(egene_pos,ATAC, gene_ATAC_window)
            print 'Openness for gene promoters:', Counter(ATAC_gene_vector)
        
            ## ATAC profile
            ATAC_eQTL = pd.DataFrame({'ATAC_SNP':ATAC_SNP_vector,'ATAC_gene':ATAC_gene_vector})
            g = pd.concat((g.reset_index(drop=True),ATAC_eQTL),axis=1)
            result_df = result_df.append(g)        
    
    return result_df






###############################################################################
###
###   Main scripts
###
###############################################################################



######   Global parameters

# DATA_DIR = '/Users/Yuan/Documents/BLab/Predict_target_genes/data'
DATA_DIR = '/scratch1/battle-fs1/heyuan/Predict_target_gene'

chr = int(os.environ['chr'])    ### do it chromosome by chromosome, because of the huge contact matrix
celltype = 'Mon'
# celltype = os.environ(['celltype'])

compute_negative_set = str_to_bool(os.environ['compute_negative_set'])

compute_contacting_degree = str_to_bool(os.environ['compute_contacting_degree'])
compute_chromosome_pos_LD = str_to_bool(os.environ['compute_chromosome_pos_LD'])

SNP_window_for_motif = int(os.environ['SNP_window_for_motif'])           ### 1kb
promoter_window_for_motif = int(os.environ['promoter_window_for_motif']) ### 2kb

merge_ATAC_data=str_to_bool(os.environ['merge_ATAC_data'])
SNP_ATAC_window = int(os.environ['SNP_ATAC_window'])           ### 1kb
gene_ATAC_window = int(os.environ['gene_ATAC_window'])           ### 2kb

merge_TF_motifs = str_to_bool(os.environ['merge_TF_motifs'])

######   compute negative sets

if compute_negative_set:
    match_distance = int(os.environ['match_distance'])
    ### read in all genes' infomation
    GRCh37_genes = pd.read_csv(os.path.join(DATA_DIR, 'eQTL/GRCh37_genes_annotated.txt'),sep='\t')
    # remove the HGxxx coordinates
    GRCh37_genes = GRCh37_genes[GRCh37_genes['CHR'] == str(chr)]
    GRCh37_genes = GRCh37_genes.drop_duplicates(subset='GENE', keep='first')
    GRCh37_genes = GRCh37_genes.fillna(0)
    GRCh37_genes.index = list(GRCh37_genes['GENE'])

    ### compute negative sets match for chr and distance
    eQTL_list = pd.read_csv(os.path.join(DATA_DIR, 'eQTL/FairFax/%s_cis_fdr05_50kb.csv' % celltype),sep='\t')
    eQTL_list = eQTL_list[eQTL_list['SNP_CHR'] == chr]
    eQTL_list_neg = negative_eQTL_onecell(eQTL_list, match_distance)
    eQTL_list_neg.to_csv(os.path.join(DATA_DIR, 'eQTL/FairFax/%s_cis_%i_fdr05_50kb_negative_annotated.csv' % (celltype,chr)),sep='\t',index=False)




######   annotate the pairs


### SNP LD blocks

def compute_LD_blocks_function(filename, ld_block_window = 50000):
    eQTL_list = pd.read_csv(filename, sep='\t')
    eQTL_list['LD_START'] = eQTL_list['SNP_POS'] - ld_block_window
    eQTL_list['LD_END'] = eQTL_list['SNP_POS'] + ld_block_window
    eQTL_list.to_csv(filename, sep='\t', index=False)




### CP-HiC

def compute_contacting_degree_function(filename, negative_set_flag):
    eQTL_list = pd.read_csv(filename,sep='\t')
    ## compute the overlap of baits/oe with the eQTL list
    if 'negative' not in filename:
        eQTL_list = eQTL_list.loc[eQTL_list.groupby('GENE')['P-VALUE'].idxmin()] 
        print "For eGenes"
    else:
        print "For random genes"
    result_pchic = eQTL_in_fragments_one_cell(eQTL_list)
    ## save results
    result_pchic.to_csv(filename,sep='\t',index=False)
    save_bed_files(result_pchic, SNP_window_for_motif, promoter_window_for_motif, negative_set_flag)


### DNase / ATAC data

def merge_ATAC_data_function(filename):
    data = pd.read_csv(filename, sep='\t')
    result_pchic_atac = eQTL_in_ATACseq_one_cell(data, DNase_data, SNP_ATAC_window, gene_ATAC_window)
    result_pchic_atac.to_csv(filename, sep='\t', index=False)



###### 1. compute the contacting map for pairs
    
if compute_contacting_degree:

    fn_true_eQTL = os.path.join(DATA_DIR, 'eQTL/FairFax/%s_cis_%i_fdr05_50kb_annotated.csv' % (celltype,chr) )
    fn_random_pairs = os.path.join(DATA_DIR, 'eQTL/FairFax/%s_cis_%i_fdr05_50kb_negative_annotated.csv' % (celltype,chr) )

    ### 1) compute the LD blocks for SNPs
    compute_LD_blocks_function(fn_true_eQTL)
    compute_LD_blocks_function(fn_random_pairs)

    ### 2) compute the contacting map

    ###    Now use all contacting pairs (could transfer to contacting E_P pairs later)
    pairs = pd.read_csv(os.path.join(DATA_DIR,'CPHiC/Interactions/PCHiC_peak_matrix_cutoff5.txt'),sep='\t', low_memory=False)
    pairs = pairs[pairs['baitChr'] == pairs['oeChr']]
    pairs = pairs[pairs['baitChr'] == str(chr)]
    pairs = pairs[list(pairs.columns[:11])+[celltype]]  
    pairs = pairs[pairs[celltype] > 0]   

    if compute_chromosome_pos_LD:

        chr_length = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431,
                       135534747, 135006516, 133851895, 115169878, 107349540, 102531392 , 90354753, 81195210, 78077248, 
                       59128983, 63025520, 48129895, 51304566]

        ### create the whole chromomse, each position corresponding to all the LD block(s)
        PIR_window = 50000
        left = np.array(pairs['oeStart']) - PIR_window
        left[left<0] = 0
        right = np.array(pairs['oeEnd']) + PIR_window
        right[right>int(chr_length[chr-1])] = int(chr_length[chr-1])
        LD_blocks = zip(left,right)
        chromosome_pos_LD = [[] for _ in xrange(int(chr_length[chr-1]))]
        for t in xrange(len(LD_blocks)):
            for tx in range(LD_blocks[t][0],LD_blocks[t][1]):
                chromosome_pos_LD[tx].append(t)
        ### convert list to array, to save time
        chromosome_pos_LD = np.array(chromosome_pos_LD)

    compute_contacting_degree_function(fn_true_eQTL, negative_set_flag=False)
    compute_contacting_degree_function(fn_random_pairs, negative_set_flag=True)








###### 2. compute the opennese map for genes

if merge_ATAC_data:
    DNase_fn = os.path.join(DATA_DIR, 'ATACseq/ENCODE/%s/%s_DNAase_GRCh37.bed' % (celltype, celltype))
    DNase_data = readin_ATACseq(DNase_fn)

    fn_true_eQTL = os.path.join(DATA_DIR, 'eQTL/FairFax/%s_cis_%i_fdr05_50kb_annotated.csv' % (celltype,chr) )
    fn_random_pairs = os.path.join(DATA_DIR, 'eQTL/FairFax/%s_cis_%i_fdr05_50kb_negative_annotated.csv' % (celltype,chr) )
    merge_ATAC_data_function(fn_true_eQTL)
    merge_ATAC_data_function(fn_random_pairs)










### SNPs in LD

LD_DIR = '/Users/Yuan/Documents/BLab/Predict_target_genes/data/LD_blocks/VCF_dataset/%s' % celltype
R2_threshold = 0.2

def SNPs_in_LD(G, snp):
    try:
        return list(G[snp])
    except:
        return []

def compute_SNPs_LD_function(chr, filename):

    ld_profile = pd.read_csv(os.path.join(LD_DIR, '%s_cis_%i_fdr05_50kb_eSNP.ld' % (celltype, chr)), sep=' ')
    ld_profile = ld_profile[ld_profile['R2'] > R2_threshold]
    ld_profile = ld_profile[ld_profile['SNP_A'] != ld_profile['SNP_B']]
    all_nodes = list(set(list(ld_profile['SNP_A']) + list(ld_profile['SNP_B'])))
    all_edges = list(ld_profile[['SNP_A','SNP_B']].apply(lambda x: tuple(x), axis=1))
    G = nx.Graph()
    G.add_nodes_from(all_nodes)
    G.add_edges_from(all_edges)

    data = pd.read_csv(os.path.join(DATA_DIR, filename), sep='\t')
    data['SNPs_LD'] = eQTL_pergene['SNP'].apply(lambda x: SNPs_in_LD(G,x))
    data.to_csv(filename, sep='\t', index=False)



compute_SNPs_LD = False

if compute_SNPs_LD:
    filename = 'eQTL/FairFax/%s_cis_%i_fdr05_50kb_annotated.csv' % (celltype,chr)    
    compute_SNPs_LD_function(chr,filename)
    filename = 'eQTL/FairFax/%s_cis_%i_fdr05_50kb_negative_annotated.csv' % (celltype,chr)    
    compute_SNPs_LD_function(chr,filename)



    

### TF motifs sites

def merge_TF_motifs_function(filename, negative_set_flag):
    result_pchic = pd.read_csv(filename, sep='\t')
    result_pchic_tfms = merge_TF_motif(result_pchic, negative_set_flag)
    result_pchic_tfms.to_csv(filename, sep='\t',index=False)

    
if merge_TF_motifs:
    ### read directly from the intermediate results
    ### and add on the TF motif information
    fn_true_eQTL = os.path.join(DATA_DIR, 'eQTL/FairFax/%s_cis_%i_fdr05_50kb_annotated.csv' % (celltype,chr))
    fn_random_pairs = os.path.join(DATA_DIR,'eQTL/FairFax/%s_cis_%i_fdr05_50kb_negative_annotated.csv' % (celltype,chr))
    merge_TF_motifs_function(fn_true_eQTL, False)
    merge_TF_motifs_function(fn_random_pairs, True)



