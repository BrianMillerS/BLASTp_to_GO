# BLASTp_to_GO_annotaions
For some model organisms, their gene ontology (GO) annotations are not as complete as more closely studied organisms. By finding proteins with similar sequences one can predict new GO anntotations based on assumed homology. This set of scripts can be used to predict novel GO annotations based on protein-protein sequence similarity using BLASTp, the outputs are parsed using R to organize the new GO annotations into a tsv.

For every query protein, if there was a protein sequence match with an evalue less than 1e-5 then it was kept and its GO terms were then associated with the matched protein. Each GO annotation produced has three fields: a GO term (from the organism with more complete anntoations), an evidence code, a gene name (from the organism with less complete annotations). Here are the first few lines from the output tsv when extending the incomplete GO annotations from T. thermophila using annotations  from A. thaliana.

![alt text](https://raw.githubusercontent.com/BrianMillerS/BLASTp_to_GO_annotaions/master/example_output.png)
