!pip install gseapy
import pandas as pd
from gseapy.plot import gseaplot 
import gseapy as gp
import numpy as np

df = pd.read_csv('input.csv')
df = df[df.pval <0.05]

df['Rank'] = -np.log10(df.pval)*df.log2FoldChange
df = df.sort_values('Rank', ascending = False).reset_index(drop = True)
df

gp.get_library_name()

ranking = df[['Gene','Rank']]
ranking

## convert the gene name as uppercase
ranking['Gene'] = ranking['Gene'].str.upper()
ranking

# select the dataset that you are interested in
pre_res =gp.prerank(rnk = ranking, gene_sets = "WikiPathways_2019_Mouse")

pre_res.results

# put result in a dataframe 
out = []

for term in list(pre_res.results):
    out.append([term,
               pre_res.results[term]['fdr'],
               pre_res.results[term]['es'],
               pre_res.results[term]['nes']])
out_df = pd.DataFrame(out, columns = ['Term','fdr','es','nes']).sort_values('fdr').reset_index(drop = True)
out_df

## save the result in local 
out_df.to_csv('WikiPathways_2019_Mouse.csv', index = False)

term_to_graph = out_df.iloc[0].Term
term_to_graph


# plot gsea
gseaplot(pre_res.ranking, term= term_to_graph, **pre_res.results[term_to_graph])

## check the gene list in specific item
print(pre_res.results['Myc Targets V1'])


## end 
### from youtube Sanbomics