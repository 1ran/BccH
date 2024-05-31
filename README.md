# BccH
This repository is for the following article:

This repository contains basic and advanced algorithms proposed in this article. 
## datasets
4 datasets are used in this paper:
* Movies：downloadable from https://www.aminer.cn/data-sna
* Yelp：downloadable from https://www.yelp.com/dataset/download
* DBLP：downloadable from https://www.aminer.org/billboard/citation
* IMDb：downloadable from https://www.imdb.com/interfaces/
* KEGG: downloadable from https://www.genome.jp/kegg/pathway.html<br />

The preprocessed data used in our experiments can be found in the "data" folder. <br />
## Rreproduce major results
### experimental environment
Compiled by Microsoft Visual Studio 2019 with MSVC++ 14.2
### algorithms
Please see c_base.cpp for the basic algorithm to compute cBC. <br />
Please see f_base.cpp for the basic algorithm to compute fBC. <br />
Please see c_ba_s2_i1.cpp for the algorithm with BA, SD2 and ID1 optimization strategies to compute cBC. <br />
Please see f_ba_i.cpp for the algorithm with BA and ID optimization strategies to compute fBC.<br />
You can run the .cpp file directly to get the experimental results.
### parameter
Before running the .cpp files, please refine the parameter.<br/>
Please refine the parameter "file" with the local dataset folder. <br />
Please refine the parameter "P" with the meta-path seleced for this dataset.<br />
