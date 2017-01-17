
# coding: utf-8

# In[1]:



# <nbformat>3.0</nbformat>

# <markdowncell>
# <codecell>


# In[1]:

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


# In[2]:

def negative_eQTL_oneSNP(eQTL, distance_threshold, egene_list):
    '''
    match the eQTL using chr and distance
    '''
    eSNP_pos = eQTL['SNP_POS']
    echr = eQTL['SNP_CHR']
    eDistance = eQTL['DISTANCE']
    egene = eQTL['GENE']
    
    genes_chr = GRCh37_genes[GRCh37_genes['Chromosome/scaffold name'] == str(echr)]
    t = genes_chr.apply(lambda x: abs(x['Gene Start (bp)'] - eSNP_pos) < distance_threshold, axis=1)
    
    if sum(t)>1:
        neg_gene = genes_chr[t].sample(1)
        while neg_gene.index[0] in egene_list:
            neg_gene = genes_chr[t].sample(1)
    else: 
        print "No negative gene found"
        
    return neg_gene.index[0]
    


# In[3]:

def negative_eQTL_onecell(eQTL, distance_threshold = 1000000):

    print 'Compute negative SNP-gene pairs match for distance and chromosome'
    
    eQTL_nega = pd.DataFrame()
    
    for name, df in eQTL.groupby('SNP_CHR'):
        print name
        eQTL_pergene = df.loc[df.groupby('GENE')['P-VALUE'].idxmin()]    ## One Lead SNP per gene
        eQTL_pergene = eQTL_pergene.reset_index(drop=True)
        negative_genes = eQTL_pergene.apply(lambda x: negative_eQTL_oneSNP(x,distance_threshold, list(eQTL_pergene['GENE'])),axis=1)
        negative_genes = GRCh37_genes.loc[np.array(negative_genes)].reset_index(drop=True)

        eQTL_nega_pergene = eQTL_pergene.merge(negative_genes, left_index=True, right_index=True)
        eQTL_nega_pergene = eQTL_nega_pergene[['SNP', 'SNP_CHR', 'SNP_POS', 'Associated Gene Name', 'Gene Start (bp)','P-VALUE']]
        eQTL_nega_pergene.columns = [['SNP', 'SNP_CHR', 'SNP_POS', 'GENE', 'GENE_START_POS','P-VALUE']]
        eQTL_nega_pergene['DISTANCE'] = abs(eQTL_nega_pergene['SNP_POS'] - eQTL_nega_pergene['GENE_START_POS'])
        
        eQTL_nega = eQTL_nega.append(eQTL_nega_pergene)
        
    return eQTL_nega


# In[4]:

def eQTL_in_fragments_one_chr(eQTL_df, contacting_df, SNP_window, promoter_window, 
                              fragments = ['baitStart', 'baitEnd', 'oeStart', 'oeEnd']):
    '''
    e-SNP near bait/oe, and promoter of the e-gene near oe/bait
    '''
    
    N = len(contacting_df)
    Start1, End1, Start2, End2 = fragments

    # e-SNP
    SNP_positions = np.reshape(eQTL_df['SNP_POS'],[len(eQTL_df),1])
    SNP_tree = spatial.KDTree(SNP_positions)

    frag1_start = np.reshape(contacting_df[Start1],[N,1])
    SNP_near_frag1Start = SNP_tree.query_ball_point(frag1_start, SNP_window)
    frag1_end = np.reshape(contacting_df[End1],[N,1])
    SNP_near_frag1End = SNP_tree.query_ball_point(frag1_end, SNP_window)

    SNP_near_frag1 = [SNP_near_frag1Start[i] + SNP_near_frag1End[i] for i in xrange(N)]


    # e-gene
    gene_start = np.reshape(eQTL_df['GENE_START_POS'], [len(eQTL_df),1])
    gene_tree = spatial.KDTree(gene_start)

    frag2_start = np.reshape(contacting_df[Start2], [N,1])
    gene_near_frag2Start = gene_tree.query_ball_point(frag2_start, promoter_window)
    frag2_end = np.reshape(contacting_df[End2], [N,1])
    gene_near_frag2End = gene_tree.query_ball_point(frag2_end, promoter_window)

    gene_near_frag2 = [gene_near_frag2Start[i] + gene_near_frag2End[i] for i in xrange(N)]
    
    SNPID = np.array(SNP_near_frag1)[np.array([len(np.intersect1d(SNP_near_frag1[t],gene_near_frag2[t]))>0 
                                               for t in xrange(len(SNP_near_frag1))])]
    geneID = np.array(gene_near_frag2)[np.array([len(np.intersect1d(SNP_near_frag1[t],gene_near_frag2[t]))>0 
                                                 for t in xrange(len(SNP_near_frag1))])]
    QTLID = [list(np.intersect1d(SNPID[t], geneID[t])) for t in xrange(len(SNPID))]
    
    return list(set([x for sublist in QTLID for x in sublist]))




# In[6]:

def save_bed_files(result_df, celltype, SNP_window, promoter_window, negative_set = False):
    '''
    save the bed files for intersection with TFB motif
    '''
    if negative_set:
        result_df.to_csv(os.path.join(DATA_DIR, 'eQTL/FairFax/%s_cis_%i_fdr05_50kb_negative_annotated.csv' % (celltype,chr)), sep='\t',index=False)
    else:
        result_df.to_csv(os.path.join(DATA_DIR, 'eQTL/FairFax/%s_cis_%i_fdr05_50kb_annotated.csv' % (celltype,chr)), sep='\t',index=False)        

    SNP_bed = []
    SNP_bed.append(['chr%i'% x for x in (result_df['SNP_CHR'])])
    SNP_bed.append(list(np.array(result_df['SNP_POS']) - SNP_window))
    SNP_bed.append(list(np.array(result_df['SNP_POS']) + SNP_window))
    SNP_bed.append(list(result_df['SNP']))
    SNP_bed = pd.DataFrame(SNP_bed).transpose()
    if negative_set:
        SNP_bed.to_csv(os.path.join(DATA_DIR, 'intermediate/%s/PCHiC_peak_matrix_cutoff5_%i_eSNP_negativeset.bed'% ("Mon", chr)),
                   sep='\t',index=False, header=False)
    else:
        SNP_bed.to_csv(os.path.join(DATA_DIR, 'intermediate/%s/PCHiC_peak_matrix_cutoff5_%i_eSNP.bed' % ("Mon", chr)),
                   sep='\t',index=False, header=False)
    
    ## 
    promoter_bed = []
    promoter_bed.append(['chr%i'% x for x in (result_df['SNP_CHR'])])
    promoter_bed.append(list(np.array(result_df['GENE_START_POS']) - promoter_window))
    promoter_bed.append(list(np.array(result_df['GENE_START_POS']) + promoter_window))
    promoter_bed.append(list(result_df['GENE']))
    promoter_bed = pd.DataFrame(promoter_bed).transpose()
    if negative_set:
        promoter_bed.to_csv(os.path.join(DATA_DIR, 'intermediate/%s/PCHiC_peak_matrix_cutoff5_%i_egene_negativeset.bed'% ("Mon", chr)),
                   sep='\t',index=False, header=False)
    else:
        promoter_bed.to_csv(os.path.join(DATA_DIR, 'intermediate/%s/PCHiC_peak_matrix_cutoff5_%i_egene.bed'% ("Mon", chr)),
                   sep='\t',index=False, header=False)
        
        
def readin_intermediate_result(celltype, SNP_window, promoter_window, negative_set=False):
    '''
    At first used to recover eQTL pairs from bed files. 
    Not used since the annotated eQTLs are saved.
    '''
    if negative_set:
        fn_SNP = 'PCHiC_peak_matrix_cutoff5_%s_eSNP_negative.bed' % chr
        fn_gene = 'PCHiC_peak_matrix_cutoff5_%s_egene_negative.bed' % chr
    else:
        fn_SNP = 'PCHiC_peak_matrix_cutoff5_%s_eSNP.bed' % chr
        fn_gene = 'PCHiC_peak_matrix_cutoff5_%s_egene.bed' % chr
    a = pd.read_csv(os.path.join(DATA_DIR, 'intermediate/%s/%s' % (celltype,fn_SNP)),sep='\t',header=None)
    b = pd.read_csv(os.path.join(DATA_DIR, 'intermediate/%s/%s' % (celltype,fn_gene)),sep='\t',header=None)
    
    inter_result = pd.concat((a,b),axis=1)
    inter_result.columns = ['SNP_CHR','SNP_POS','non','SNP','GENE_chr','GENE_START_POS','NAN','GENE']
    inter_result = inter_result[['SNP_CHR','SNP_POS','SNP','GENE_chr','GENE_START_POS','GENE']]
    inter_result['SNP_POS'] = inter_result['SNP_POS'] + SNP_window
    inter_result['GENE_START_POS'] = inter_result['GENE_START_POS'] + promoter_window
    return inter_result

# In[ ]:

#### Need to scan the bed files for known motifs using fimo
#### run the following codes in DATA_DIR
#### bash TFBS.sh Mon eSNP
#### bash TFBS.sh Mon gene


# In[ ]:

def merge_TF_motif(result_df, celltype, negative_set=False):
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


# In[78]:

def eQTL_in_ATACseq_one_chr(eQTL_pos, ATAC_regions, window = 1000):

    N = len(ATAC_regions)
    eQTL_tree = spatial.KDTree(eQTL_pos)

    frag_start = np.reshape(ATAC_regions['Start'],[N,1])
    eQTL_near_fragStart = eQTL_tree.query_ball_point(frag_start, window)
    frag_end = np.reshape(ATAC_regions['End'],[N,1])
    eQTL_near_fragEnd = eQTL_tree.query_ball_point(frag_end, window)

    eQTL_near_frag = [list(set(eQTL_near_fragStart[i] + eQTL_near_fragEnd[i])) for i in xrange(N)]
    eQTLID = [x for x in eQTL_near_frag if len(x) > 0]
    eQTLID = list(set([a for b in eQTLID for a in b]))
    
    x = np.zeros(len(eQTL_pos))
    x[eQTLID] = 1
    return x


# In[119]:

def eQTL_in_ATACseq_one_cell(eQTL, ATAC, window_SNP = 1000, window_gene = 2000000):
    
    result_df = pd.DataFrame()
    for name, g in eQTL.groupby('SNP_CHR'):
        ## eSNP
        eSNP_pos = np.reshape(g['SNP_POS'],[len(g),1])
        ATAC_SNP_vector = eQTL_in_ATACseq_one_chr(eSNP_pos,ATAC[ATAC['Chr'] == name], window = window_SNP)
        ## egene
        egene_pos = np.reshape(g['GENE_START_POS'],[len(g),1])
        ATAC_gene_vector = eQTL_in_ATACseq_one_chr(egene_pos,ATAC[ATAC['Chr'] == name], window = window_gene)
        
        ## ATAC profile
        ATAC_eQTL = pd.DataFrame({'ATAC_SNP':ATAC_SNP_vector,'ATAC_gene':ATAC_gene_vector})
        g = pd.concat((g.reset_index(drop=True),ATAC_eQTL),axis=1)
        result_df = result_df.append(g)        
        
    return result_df



# In[120]:

#### ATAC-seq

# ATACseq = pd.read_csv(os.path.join(DATA_DIR, 'ATACseq/GSE74912_ATACseq_All_Counts.txt'),sep='\t')



def SNP_in_PIR(eQTL, contacting_df, PIR_window = 1000):
    gene_promoter = eQTL['GENE']
    SNP = eQTL['SNP']
    SNP_POS = eQTL['SNP_POS']

    gene_promoter_idx = np.where([gene_promoter in str(x).split(';') for x in list(contacting_df['baitName'])])[0]
    gene_promoter_window = contacting_df.iloc[gene_promoter_idx]

    # if len(gene_promoter_window)==0:
    #     print gene_promoter

    left = np.array(gene_promoter_window['oeStart']) - PIR_window
    right = np.array(gene_promoter_window['oeEnd']) + PIR_window

    if sum([(t[0] < SNP_POS and t[1] > SNP_POS) for t in zip(left,right)]) > 0:
        return 1
    else:
        return 0
    



def eQTL_in_fragments_one_cell(eQTL, pairs, cell, SNP_window, promoter_window):

    print 'Compute whether the eQTL pair fall in any promoter-PIR pair'
    
    contacting_cell = pairs[list(pairs.columns[:11])+[cell]]  
    contacting_cell = contacting_cell[contacting_cell[cell] > 0]   

    result_df = pd.DataFrame()
    for name, df in eQTL.groupby('SNP_CHR'):
        contacting_pairs = contacting_cell[contacting_cell['baitChr'] == str(name)]
        eQTL_pergene = df.loc[df.groupby('GENE')['P-VALUE'].idxmin()]    ## One Lead SNP per gene
        
        ## approach 1: use KDTree: not right
        # functional_contacting_idx1 = eQTL_in_fragments_one_chr(eQTL_pergene, contacting_pairs, SNP_window, promoter_window)
        # functional_contacting_idx2 = eQTL_in_fragments_one_chr(eQTL_pergene, contacting_pairs, SNP_window, promoter_window,
        #                                                        fragments = ['oeStart', 'oeEnd','baitStart', 'baitEnd'])
        # contact_vector = np.array([0] * len(eQTL_pergene))
        # contact_vector[np.array(list(set(functional_contacting_idx1+functional_contacting_idx2)))] = 1

        ## approach 2: can only compute the direct contacts
        # functional_contacting_idx = eQTL_pergene.apply(lambda x: SNP_in_PIR(x, contacting_pairs),axis=1)
        # contact_vector = np.array(functional_contacting_idx)
        # eQTL_pergene['contacting'] = contact_vector

        ## approach 3: compute the degree of contacts
        frag2 = compute_frag_by_frag(contacting_pairs,cell) 
        degree_contact = degree_of_contact_SNP_gene(eQTL_pergene, contacting_pairs, frag2)

        contact_vector = np.zeros(len(eQTL_pergene))
        contact_vector[np.where(map(lambda x: 2 in x, degree_contact))[0]] = 2
        contact_vector[np.where(map(lambda x: 1 in x, degree_contact))[0]] = 1

        eQTL_pergene['contacting'] = contact_vector

        save_bed_files(eQTL_pergene, 'Mon', SNP_window_for_motif, promoter_window_for_motif, negative_set_flag)
        
        result_df = result_df.append(eQTL_pergene)

        print name
        print Counter(contact_vector)
                
    return result_df



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
    gg_in_frag = np.where(map(lambda x: gg in str(x), frag_by_frag.keys()))[0]
    if len(gg_in_frag) > 0:
        temp = np.array(frag_by_frag.keys())[np.array(gg_in_frag)]
        return np.random.choice(temp)
    else:
        return 'NOGENE' 

def degree_of_contact_SNP_gene(eQTL_df, contacting_df, frag_by_frag, PIR_window=2000):

    ### find the other fragments where eSNP locate in 
    left = np.array(contacting_df['oeStart']) - PIR_window
    right = np.array(contacting_df['oeEnd']) + PIR_window
    SNPs_PIR = eQTL_df.apply(lambda x: [x['SNP_POS'] in range(t[0],t[1]) for t in zip(left,right)], axis=1)
    SNPs_PIR = SNPs_PIR.apply(lambda x: list(np.where(x)[0]))   ## a list of lists
    eQTL_df['SNP_contacting_df_idx'] = list(SNPs_PIR)
    
    ### find the degree of contact by using frag_by_frag dictionary
    degree_contact = eQTL_df.apply(lambda eq: [frag_by_frag[find_the_key(frag_by_frag, eq['GENE'])][contacting_df.iloc[t]['oeID']] 
                                               for t in eq['SNP_contacting_df_idx']], axis=1)

    
    return degree_contact






# In[ ]:

DATA_DIR = '/Users/Yuan/Documents/BLab/Predict_target_genes/data'
# DATA_DIR = '/scratch1/battle-fs1/heyuan/Predict_target_gene'

compute_negative_set = str_to_bool(os.environ['compute_negative_set'])
match_distance = int(os.environ['match_distance'])

from_beginning = str_to_bool(os.environ['from_beginning'])
negative_set_flag = str_to_bool(os.environ['negative_set_flag'])

SNP_window_for_PHiC = int(os.environ['SNP_window_for_PHiC'])
promoter_window_for_PHiC = int(os.environ['promoter_window_for_PHiC'])


chr = int(os.environ['chr'])

SNP_window_for_motif = 1000
promoter_window_for_motif = 2000

if compute_negative_set:
    ### read in all genes infomation
    GRCh37_genes = pd.read_csv(os.path.join(DATA_DIR, 'eQTL/GRCh37_genes.txt'),sep='\t')
    gene_idx = []
    for name, df in GRCh37_genes.groupby(['Associated Gene Name']):
        gene_idx.append(np.random.choice(list(df.index)))
    GRCh37_genes = GRCh37_genes.iloc[np.array(gene_idx)]
    GRCh37_genes.index = GRCh37_genes['Associated Gene Name']
    
    ### compute negative sets match for chr and distance
    mon_eQTL = pd.read_csv(os.path.join(DATA_DIR, 'eQTL/FairFax/monocytes_cis_fdr05_50kb.csv'),sep='\t')
    mon_eQTL = mon_eQTL[mon_eQTL['SNP_CHR'] == chr]
    mon_eQTL_neg = negative_eQTL_onecell(mon_eQTL, match_distance)
    mon_eQTL_neg.to_csv(os.path.join(DATA_DIR, 'eQTL/FairFax/monocytes_cis_%i_fdr05_50kb_negative.csv' % chr),sep='\t',index=False)

    
if from_beginning:
    ### First use all contacting pairs (could transfer to contacting E_P pairs later)
    pairs = pd.read_csv(os.path.join(DATA_DIR,'CPHiC/Interactions/PCHiC_peak_matrix_cutoff5.txt'),sep='\t', low_memory=False)
    pairs = pairs[pairs['baitChr'] == pairs['oeChr']]
    
    ### read in the (negative) eQTL list
    if negative_set_flag:
        mon_eQTL = pd.read_csv(os.path.join(DATA_DIR, 'eQTL/FairFax/monocytes_cis_%i_fdr05_50kb_negative.csv' % chr),sep='\t')
    else:
        mon_eQTL = pd.read_csv(os.path.join(DATA_DIR, 'eQTL/FairFax/monocytes_cis_fdr05_50kb.csv'),sep='\t')
    mon_eQTL = mon_eQTL[mon_eQTL['SNP_CHR'] == chr]

    ### compute the overlap of baits/oe with the eQTL list
    print "Mon"
    Mon_result = eQTL_in_fragments_one_cell(mon_eQTL, pairs, 'Mon',SNP_window_for_PHiC, promoter_window_for_PHiC)
    temp = map(lambda x:len(x)-2,list(Mon_result['SNP_contacting_df_idx']))
    print 'There are %i SNPs among %i that fall in a fragment (2000 window)' % (sum(np.array(temp) >0),len(temp))
    save_bed_files(Mon_result, 'Mon', SNP_window_for_motif, promoter_window_for_motif, negative_set_flag)
    
else:
    ### read directly from the intermediate results
    ### and add on the TF motif information
    if negative_set_flag:
        filename = 'eQTL/FairFax/%s_cis_%i_fdr05_50kb_negative_annotated.csv' % ('Mon',chr)
    else:
        filename = 'eQTL/FairFax/%s_cis_%i_fdr05_50kb_annotated.csv' % ('Mon',chr)
    Mon_result = pd.read_csv(os.path.join(DATA_DIR, filename), sep='\t')
    Mon_result = merge_TF_motif(Mon_result, 'Mon', negative_set_flag)
    Mon_result.to_csv(os.path.join(DATA_DIR, filename), sep='\t',index=False)


### Merge the TF traits with eQTL pairs






