关键基因集处理：
原始数据在The database of Online GEne Essentiality (OGEE) (downloaded at 20/10/2020)  http://ogee.medgenius.info/browse/中下载，共131823列，21556种基因
对原始数据去重，多个相同的locus中只要存在一个以上为E，则essential置为E，全为NE则置为NE
将locus放进uniprot中映射出对应gene ，最终得到18900个gene，7123个E，11777个NE

PPI网络处理：
BioGRID (Version 3.5.182)  http://thebiogrid.org/  包含1748436条边
1）选取两边ORGANISM_ID均为9606后得到humo的bio网络含504848条边，18900个点
2）与关键基因集做交集，去除自环和重复边后为322406条边和15721个点

基因表达数据：
Gene Expression Omnibus (GEO) database   https://www.ncbi.nlm.nih.gov/geo/
GSE86354提供了基因型-组织表达(GTEx)项目产生的8个组织位点的1,558份样本的表达谱
其中Bladder有11列特征