
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
sys.setrecursionlimit(10000)


# In[2]:

def negative_eQTL_oneSNP(eQTL, distance_threshold):
    '''
    match the eQTL using chr and distance
    '''
    eSNP_pos = eQTL['SNP_POS']
    echr = eQTL['SNP_CHR']
    eDistance = eQTL['DISTANCE']
    egene = eQTL['GENE']
    
    genes_chr = GRCh37_genes[GRCh37_genes['Chromosome/scaffold name'] == str(echr)]
    t = genes_chr.apply(lambda x: abs(x['Gene Start (bp)'] - eSNP_pos) < distance_threshold, axis=1)
    neg_gene = genes_chr[t].sample(1)
    
    if sum(t)>1:
        while neg_gene.index[0] == egene:
            neg_gene = genes_chr[t].sample(1)
    else: 
        print "No negative gene found"
        
    return neg_gene.index[0]
    


# In[3]:

def negative_eQTL_onecell(eQTL, distance_threshold = 1000000):
    
    eQTL_nega = pd.DataFrame()
    
    for name, df in eQTL.groupby('SNP_CHR'):
        eQTL_pergene = df.loc[df.groupby('GENE')['P-VALUE'].idxmin()]    ## One Lead SNP per gene
        eQTL_pergene = eQTL_pergene.reset_index(drop=True)
        negative_genes = eQTL_pergene.apply(lambda x: negative_eQTL_oneSNP(x,distance_threshold=distance_threshold),axis=1)
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
    
    QTLID = np.where([len(np.intersect1d(SNP_near_frag1[t],gene_near_frag2[t]))>0 for t in xrange(len(SNP_near_frag1))])[0]

    return list(QTLID)







# In[5]:

def eQTL_in_fragments_one_cell(eQTL, pairs, cell, SNP_window, promoter_window):
    
    contacting_cell = pairs[list(pairs.columns[:11])+[cell]]  
    contacting_cell = contacting_cell[contacting_cell[cell] > 0]   

    result_df = pd.DataFrame()
    for name, df in eQTL.groupby('SNP_CHR'):
        contacting_pairs = contacting_cell[contacting_cell['baitChr'] == str(name)]
        eQTL_pergene = df.loc[df.groupby('GENE')['P-VALUE'].idxmin()]    ## One Lead SNP per gene
        
        functional_contacting_idx1 = eQTL_in_fragments_one_chr(eQTL_pergene, contacting_pairs, SNP_window, promoter_window)
        functional_contacting_idx2 = eQTL_in_fragments_one_chr(eQTL_pergene, contacting_pairs, SNP_window, promoter_window,
                                                               fragments = ['oeStart', 'oeEnd','baitStart', 'baitEnd'])
        contact_vector = np.array([0] * len(eQTL_pergene))
        contact_vector[np.array(list(set(functional_contacting_idx1+functional_contacting_idx2)))] = 1
        eQTL_pergene['contacting'] = contact_vector
        
        result_df = result_df.append(eQTL_pergene)
        
        print 'There are %i contacting eQTL pairs among %i for chr%i' % (sum(contact_vector),len(eQTL_pergene), name)
                
    return result_df



# In[6]:

def save_bed_files(result_df, celltype, SNP_window, promoter_window, negative_set = False):
    '''
    save the bed files for intersection with TFB motif
    '''
    
    SNP_bed = []
    SNP_bed.append(['chr%i'% x for x in (result_df['SNP_CHR'])])
    SNP_bed.append(list(np.array(result_df['SNP_POS']) - SNP_window))
    SNP_bed.append(list(np.array(result_df['SNP_POS']) + SNP_window))
    SNP_bed.append(list(result_df['SNP']))
    SNP_bed = pd.DataFrame(SNP_bed).transpose()
    if negative_set:
        promoter_bed.to_csv(os.path.join(DATA_DIR, 'intermediate/%s/PCHiC_peak_matrix_cutoff5_%s_eSNP_negativeset.bed'% (celltype, celltype)),
                   sep='\t',index=False, header=False)
    else:
        SNP_bed.to_csv(os.path.join(DATA_DIR, 'intermediate/%s/PCHiC_peak_matrix_cutoff5_%s_eSNP.bed' % (celltype, celltype)),
                   sep='\t',index=False, header=False)
    
    ## 
    promoter_bed = []
    promoter_bed.append(['chr%i'% x for x in (result_df['GENE_CHR'])])
    promoter_bed.append(list(np.array(result_df['GENE_START_POS']) - promoter_window))
    promoter_bed.append(list(np.array(result_df['GENE_START_POS']) + promoter_window))
    promoter_bed.append(list(result_df['GENE']))
    promoter_bed = pd.DataFrame(promoter_bed).transpose()
    if negative_set:
        promoter_bed.to_csv(os.path.join(DATA_DIR, 'intermediate/%s/PCHiC_peak_matrix_cutoff5_%s_egene_negativeset.bed'% (celltype, celltype)),
                   sep='\t',index=False, header=False)
    else:
        promoter_bed.to_csv(os.path.join(DATA_DIR, 'intermediate/%s/PCHiC_peak_matrix_cutoff5_%s_egene.bed'% (celltype, celltype)),
                   sep='\t',index=False, header=False)
        
        
def readin_intermediate_result(celltype,SNP_window, promoter_window, negative_set=False):
    if negative_set:
        fn_SNP = 'PCHiC_peak_matrix_cutoff5_%s_eSNP_negativese.bed' % celltype
        fn_gene = 'PCHiC_peak_matrix_cutoff5_%s_egene_negativese.bed' % celltype
    else:
        fn_SNP = 'PCHiC_peak_matrix_cutoff5_%s_eSNP.bed' % celltype
        fn_gene = 'PCHiC_peak_matrix_cutoff5_%s_egene.bed' % celltype
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
#### bash TFBS.sh tCD4 eSNP
#### bash TFBS.sh tCD4 gene


# In[ ]:

def merge_TF_motif(celltype, result_df):
    # read in eQTLs with TF motif annotated
    eSNP_tf = pd.read_csv(os.path.join(DATA_DIR,'intermediate/%s' % celltype,
                                       'fimo.output.%s.eSNP.bed' % celltype),sep='\t',header=None)
    egene_tf = pd.read_csv(os.path.join(DATA_DIR,'intermediate/%s' % celltype,
                                        'fimo.output.%s.egene.bed' % celltype),sep='\t',header=None)
    eSNP_tf.columns = ['SNPID','chr','start','end','motif_SNP']
    egene_tf.columns = ['geneID','chr','start','end','motif_gene']

    # merge
    TF_df = result_df.merge(eSNP_tf[['SNPID','motif_SNP']], left_on="SNP", right_on='SNPID',how='left')
    TF_df = TF_df.merge(egene_tf[['geneID','motif_gene']], left_on="GENE", right_on='geneID',how='left')
    TF_df = TF_df[['SNP','SNP_CHR','GENE','SNP_POS','GENE_START_POS','motif_SNP','motif_gene']]
    return TF_df


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
# ATAC_CD4 = pd.concat((ATACseq[['Chr','Start','End']],
#                      ATACseq[[x for x in ATACseq.columns if 'CD4' in x]].mean(axis=1)),axis=1)
# tCD4_TF_ATAC = eQTL_in_ATACseq_one_cell(tCD4_TF, ATACseq)



# In[ ]:



# In[ ]:

DATA_DIR = '/Users/Yuan/Documents/BLab/Predict_target_genes/data'

compute_negative_set = True
match_distance = 1000000

from_beginning=True
negative_set_flag = True

SNP_window_for_PHiC = 1000
promoter_window_for_PHiC = 1000000

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
    CD4_eQTL = pd.read_csv(os.path.join(DATA_DIR, 'eQTL/Raj/cd4T_cis_fdr05.tsv'),sep='\t')
    mon_eQTL = pd.read_csv(os.path.join(DATA_DIR, 'eQTL/Raj/monocytes_cis_fdr05.tsv'),sep='\t')
    CD4_eQTL_neg = negative_eQTL_onecell(CD4_eQTL, match_distance)
    CD4_eQTL_neg.to_csv(os.path.join(DATA_DIR, 'eQTL/Raj/cd4T_cis_fdr05_negative.tsv'),sep='\t',index=False)
    mon_eQTL_neg = negative_eQTL_onecell(mon_eQTL, match_distance)
    mon_eQTL_neg.to_csv(os.path.join(DATA_DIR, 'eQTL/Raj/monocytes_cis_fdr05_negative.tsv'),sep='\t',index=False)

    
if from_beginning:
    ### First use all contacting pairs (could transfer to contacting E_P pairs later)
    pairs = pd.read_csv(os.path.join(DATA_DIR,'CPHiC/Interactions/PCHiC_peak_matrix_cutoff5.txt'),sep='\t', low_memory=False)
    pairs = pairs[pairs['baitChr'] == pairs['oeChr']]
    
    ### read in the (negative) eQTL list
    if negative_set_flag:
        CD4_eQTL = pd.read_csv(os.path.join(DATA_DIR, 'eQTL/Raj/cd4T_cis_fdr05_negative.tsv'),sep='\t')
        mon_eQTL = pd.read_csv(os.path.join(DATA_DIR, 'eQTL/Raj/monocytes_cis_fdr05_negative.tsv'),sep='\t')
    else:
        CD4_eQTL = pd.read_csv(os.path.join(DATA_DIR, 'eQTL/Raj/cd4T_cis_fdr05.tsv'),sep='\t')
        mon_eQTL = pd.read_csv(os.path.join(DATA_DIR, 'eQTL/Raj/monocytes_cis_fdr05.tsv'),sep='\t')
        
    ### compute the overlap of baits/oe with the eQTL list
    print "tCD4"
    tCD4_result = eQTL_in_fragments_one_cell(CD4_eQTL, pairs, 'tCD4',SNP_window_for_PHiC, promoter_window_for_PHiC)
    save_bed_files(tCD4_result, 'tCD4', SNP_window_for_motif, promoter_window_for_motif, negative_set_flag)
    print "Mon"
    Mon_result = eQTL_in_fragments_one_cell(mon_eQTL, pairs, 'Mon',SNP_window_for_PHiC, promoter_window_for_PHiC)
    save_bed_files(Mon_result, 'Mon', SNP_window_for_motif, promoter_window_for_motif, negative_set_flag)
    
else:
    ### read directly from the intermediate results
    tCD4_result = readin_intermediate_result('tCD4', SNP_window_for_motif, promoter_window_for_motif, negative_set_flag)
    Mon_result = readin_intermediate_result('Mon', SNP_window_for_motif, promoter_window_for_motif, negative_set_flag)




