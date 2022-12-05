# General framework:
* Overall steps: \
    ```https://github.com/aertslab/pySCENIC/blob/c120979378f2e6a1415b1d9e973b19a718531484/src/pyscenic/cli/pyscenic.py```
    * This script overall descirbe the procedure that run the Pyscenic file

* Network inference:
from arboreto.algo import grnboost2 \
For this step, we are running the lightGBM model which consists several shallow trees that use the boosting method to make the prediction for each target genes. Note we will have many target genes for this step, so we are running parallel for many lightGBM model and it will take a while to finish. After that, for each target gene, we compute the feature importance for each transcription factor and take the average for all the rows in the expression data in terms of the importance score. For the detail about how the feature importance score is calculated, please see similar model's explanation: ```https://scikit-learn.org/stable/auto_examples/ensemble/plot_forest_importances.html```. Basically, the features for internal nodes are selected with some criterion, which for classification tasks can be gini impurity or infomation gain, and for regression is variance reduction. We can measure how each feature decrease the impurity of the split (the feature with highest decrease is selected for internal node). For each feature we can collect how on average it decreases the impurity. The average over all trees in the forest is the measure of the feature importance. Also note this step is not deterministic, which means results changes everytime you run the code. A solution they mentioned might be running the code many times and take the average, time-consuming also?

* Module generation:\
a.	Firstly, we generate some target-TF pairs(which they called module in the Pyscenic paper), which solely based on the top N targets for each factor, with a default of 50 targets
```https://github.com/aertslab/pySCENIC/blob/c120979378f2e6a1415b1d9e973b19a718531484/src/pyscenic/utils.py#L219```
b.	Secondly, we calculate the pearson correlation for each row(gene). For the below codes, this is a snippet about how they implement the module "filtering". So they use the previous results which is feature importance csv file, then they extract all the genes including target and TFs. Then they compute the correlation between each gene(column). The result is a correlation matrix which looks like [[1,0.xx,0.xx,...],[0.xx,1,0.xx,...],[...],[...]]. Then they try to delete those negative correlation between target and TF pairs, which they think is meaningless(negative correlation probably means this TF doesn't contribute to the cell's functioning?, see the rhos step for the deletion). At the end of this step, we get the target-TF pairs that could have positive correlation solely based on the expression matrix' correlation, which their paper mentioned also. 

```Genes = list(set(adjacencies[COLUMN_NAME_TF]).union(set(adjacencies[COLUMN_NAME_TARGET]))) 
ex_mtx = ex_mtx[ex_mtx.columns[ex_mtx.columns.isin(genes)]] 
corr_mtx = pd.DataFrame(index=ex_mtx.columns, columns=ex_mtx.columns, data=np.corrcoef(ex_mtx.values.T)) 
rhos = np.array([corr_mtx[s2][s1] for s1, s2 in zip(adjacencies.TF, adjacencies.target)]) 
*Drop negative correlation 
*set some threshold
```
* Pruning(CisTarget):  
a.	For this step, since the previous step is bascially based on the expression data, which might or might not a representation of the true relationship between TF and target pair. So in this step, they add additional step that I consider further extract the meaningful TF-target pairs. Given a pre-definied database of ranking for the gene(n_feature,n_genes,pre-defined), use the np.cumsum(row) for each gene to get the recovery curve. For the recovery curve, see the below code snippet and github link for further information. Basically, given a ranking dataframe, they calculated the occurance frequency for each module(tf-target pair) and finally calculate the cumulative sum for those frequencies.Then they use the ```np.sum``` to get the auc for each gene based on some rank-cut-off threshold. This remain the most challenging step for understanding the Pyscenic, I still have some confusions yet.
```https://github.com/aertslab/ctxcore/blob/e68a0b168b9511efb3140c5cf1dd2b07bbd811cc/src/ctxcore/recovery.py#L64```
```
n_features = rankings.shape[0] 
rccs = np.empty(shape=(n_features, rank_threshold))  # Pre-allocation. 
for row_idx in range(n_features):
    curranking = rankings[row_idx, :]  
    rccs[row_idx, :] = np.cumsum( 
        np.bincount(curranking, weights=weights)[:rank_threshold]) 
return rccs
```


* Auc cell \
a.	Do the same recovery curve calculation as the previous step for each cell each module, which means each module represent the same regulon that control that cell. For this step, ranking are created solely based on the expression matrix also for each cell, just take each row and compare the number on that row (See code below). After we get those rankings for each pairs, then we could measure the cell's functionality based on those TF's enrichment implicitly
```
def create_rankings(ex_mtx: pd.DataFrame, seed=None) -> pd.DataFrame: 
    """
    Create a whole genome rankings dataframe from a single cell expression profile dataframe.
    :param ex_mtx: The expression profile matrix. The rows should correspond to different cells, the columns to different 
        genes (n_cells x n_genes). 
    :return: A genome rankings dataframe (n_cells x n_genes).
    """
    # Do a shuffle would be nice for exactly similar behaviour as R implementation. 
    # 1. Ranks are assigned in the range of 1 to n, therefore we need to subtract 1. 
    # 2. In case of a tie the 'first' method is used, i.e. we keep the order in the original   array. The remove any 
    #    bias we shuffle the dataframe before ranking it. This introduces a performance penalty!  
    # 3. Genes are ranked according to gene expression in descending order, i.e. from highly expressed (0) to low expression (n). 
    # 3. NAs should be given the highest rank numbers. Documentation is bad, so tested   implementation via code snippet: 
    # 
    #    import pandas as pd 
    #    import numpy as np 
    #    df = pd.DataFrame(data=[4, 1, 3, np.nan, 2, 3], columns=['values']) 
    #    # Run below statement multiple times to see effect of shuffling in case of a tie. 
    #    df.sample(frac=1.0, replace=False).rank(ascending=False, method='first' na_option='bottom').sort_index() - 1 
    # 
    return ( 
        ex_mtx.sample(frac=1.0, replace=False, axis=1, random_state=seed) 
        .rank(axis=1, ascending=False, method='first', na_option='bottom') 
        .astype(DTYPE) 
        - 1 
    ) 
```
b.	Filter out some low-threshold cell
c.	```https://github.com/aertslab/pySCENIC/blob/c120979378f2e6a1415b1d9e973b19a718531484/src/pyscenic/aucell.py```
d.	

* Clustering based on the regulon and get the cell  
* For this step, the results from the previous step looks like some cell filled with aucs with different modules. Then we could do clustering for those cells, we could do UMAP(which is deterministic(result don't change when seed change), TSNE(non-deterministic, result do change when seed change), both of those are non-linear method, or we could do k-means, which is linear method, which might not capture the necessary important information for the results from LightGBM(Genius 3))


Read the code through but have a little bit confusion about the AUC and how to get access to the predefined ranking database?
