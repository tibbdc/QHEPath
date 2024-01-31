# Unveiling Metabolic Engineering Strategies by Quantitative Heterologous Pathway Design
 
 ### Description
 we develop a high-quality cross-species metabolic network (CSMN) model and a quantitative heterologous pathway design algorithm (QHEPath) to predict  the potential for improvement in product pathway yield in different species. This algorithm can compute calculate multiple suboptimal and optimal pathways for the target product.

### Installation

To create a stand-alone environment named qhepath with Python and all the required package versions, run the following code:

```shell
$ conda create -n qhepath python=3.6.13
```
```shell
$ conda activate qhepath
```
```shell
$ pip install ipykernel
```
```shell
$ python -m ipykernel install --user --name qhepath --display-name "qhepath"
```
```shell
$ pip install -r requirements.txt
```

### Usage

### 1. Reconstructing a high-quality cross-species metabolic network model (CSMN)
+ quality_control_workflow.ipynb

### 2. Calculate multiple suboptimal and optimal pathways for the target product
+ example: 
    + chassis：iML1515 (E. coli)
    + substrate: glc__D_c (Glucose)
    + product: ac_c (Acetate)
    + number of optimization: 10
```shell
$ python pathway_calculation.py iML1515 glc__D_c ac_c 10
```


