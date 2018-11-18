# BLASTp_to_GO_annotaions
For some model organisms, their gene ontology (GO) annotations are not as complete as more closely studied organisms. By finding proteins with similar sequences one can predict new GO anntotations based on assumed homology. This set of scripts can be used to predict novel GO annotations based on protein-protein sequence similarity between multiple organisms, and then organize this information using R to generate a tsv file with the new annotations.

For every query protein, if there was protein sequence match with an evalue less than 1e-5 then it was kept and its GO terms were then associated with the matched protein. Each GO annotation produced has three fields: a gene name (from the organism with less complete annotations), an evidence code, and an GO term (from the organism with more complete anntoations). Here are the first few lines from the output tsv when extending the incomplete GO annotations from T. thermophila using annotations  from A. thaliana.

![alt text](https://raw.githubusercontent.com/BrianMillerS/BLASTp_to_GO_annotaions/master/example_output.png)
