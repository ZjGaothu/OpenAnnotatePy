# OpenAnnotatePy

#### A python package for efficiently annotating the chromatin accessibility of genomic regions

Chromatin accessibility is a measure of the ability of nuclear macromolecules to physically contact DNA, and is essential for understanding regulatory mechanisms.

OpenAnnotate facilitates the chromatin accessibility annotation for massive genomic regions by allowing ultra-efficient annotation across various biosample types based on chromatin accessibility profiles accumulated in public repositories (1236 samples from ENCODE and 1493 samples from ATACdb).

For more information, please refer to the web: http://health.tsinghua.edu.cn/openannotate/

### Install OpenAnnotate via Pypi
Anaconda users can first create a new Python environment and activate it via(this is unnecessary if your Python environment is managed in other ways)

```Python
conda create python=3.6 -n OpenAnnotatePy
conda activate OpenAnnotatePy
```

OpenAnnotate is available on pypi here and can be installed via
```Python
pip install OpenAnnotatePy
```

### Functions of an Annotate() object

| Code | Function |
| ------ | ------  |
| testWebserver() | test whether the web server is working normally |
| setAddress(IP, port) | set the address of the web server |
| help() | get a list of the various functions and arguments that the package contains. |
| getParams() | get params list |
| getCelltypeList(protocol, species) | get cell type list |
| getTissueList(protocol, species) | get cell type list |
| getSystemList(protocol, species) | get cell type list |
| searchCelltype(protocol, species, keyword) | search for cell types that contain keyword |
| searchTissue(protocol, species, keyword) | search for cell types that contain keyword |
| searchSystem(protocol, species, keyword) | search for cell types that contain keyword |
| setParams(assay, species, cell_type, perbase)| set parameters |
| runAnnotate(input) | upload file to server |
| getProgress(task_id)| you can view the annotation progress |
| getAnnoResult(result_type,task_id,cell_type) | download the annotation result |
| getInputFile(save_path, task_id) | get your input file from server |
| viewParams(task_id) | view parameters|
| getExampleTaskID() | get example task id|
| getExampleInputFile(save_path) | get example input file to the save_path|
| fromOpen2EpiScanpy(data, head) | generate anndata from annotation result |

### A simple example

Upload a region file to the web server and download the head and the readopen of the annotation result to the local path, then initialize an anndata for downstream analysis.

```python
from OpenAnnotatePy import OpenAnnotateApi
oaa=OpenAnnotateApi.Annotate()

oaa.setParams(species=1, protocol=1, cell_type=1, perbase=1)

task_id=oaa.runAnnotate(input='./EXAMPLE.bed.gz')

anno_data = oaa.getAnnoResult(result_type = 2,task_id = task_id ,cell_type = 1)

anno_head = oaa.getAnnoResult(result_type = 1,task_id = task_id ,cell_type = 1)

ann_data = fromOpen2EpiScanpy(data = anno_data, head = anno_head)

```


### Usage  

**Import**

The package inclues a class named `OpenAnnotatePy`, All functions are implemented by instantiating objects of this class.
```Python
from OpenAnnotatePy import OpenAnnotateApi
```

**Instantiate object**

Instantiate an object with the data path.

```Python
oaa=OpenAnnotateApi.Annotate()
```


**Help**

Get a list of the various functions and arguments that the package contains.
```Python
oaa.help()
'''
testWebserver() : test whether the web server is working normally
setAddress(IP, port) : set the address of the web server
getParams() : get params list
getCelltypeList(protocol,species) : get cell type list
getTissueList(protocol,species) : get tissue list
getSystemList(protocol,species) : get system list

searchCelltype(protocol, species, keyword) : search for cell types that contain keyword
searchTissue(protocol, species, keyword) : search for tissues that contain keyword and the corresponding cell types
searchSystem(protocol, species, keyword) : search for systems that contain keyword and the corresponding cell types
setParams(assay,species,cell_type,perbase) : set params list

runAnnotate(input) : Upload file to server
getProgress(task_id) : query the annotation progress
getAnnoResult(result_type,task_id,cell_type) : download annotation result to local path
getInputFile(save_path, task_id) : get your input file from server
viewParams(task_id) : view parameters
getExampleTaskID() : get example task id
getExampleInputFile(save_path) : get example input file to the save_path
fromOpen2EpiScanpy(data, head) : generate anndata from annotation result
'''
```

**Get parameters**

Get the parameters to be set.
```python
# get basic parameters you need to set
oaa.getParams()

# get the corresponding cell type list
oaa.getCelltypeList(protocol, species)

# get the corresponding tissues list
oaa.getTissueList(protocol, species)

# get the corresponding systems list
oaa.getSystemList(protocol, species)

# search cell type
oaa.searchCelltype(protocol, species, keyword)

# search tissue and corresponding cell types
oaa.searchTissue(protocol, species, keyword)

# search system and corresponding cell types
oaa.searchSystem(protocol, species, keyword)
```
- `getParams()`: Return the parameter list of `species`,`protocol` and `Annotate method`.
- `getCelltypeList(protocol,species)` : Return the cell type list of the corresponding `protocol` and `species`.
- `species` : 
  - 1 : GRCh37/hg19 
  - 2 : GRCh38/hg38 
  - 3 : GRCm37/mm9 
  - 4 : GRCm38/mm10
- `protocol`: 
  - 1 : DNase-seq(ENCODE)
  - 2 : ATAC-seq(ENCODE) 
  - 3 : ATAC-seq(ATACdb)
- `keyword`: Key word for search. Such as `K562` and `Blood`.


**Set parameters**

Set parameters for your object.

```python
oaa.setParams(species, protocol, cell_type, perbase)
```
- `species` : 
  - 1 : GRCh37/hg19 
  - 2 : GRCh38/hg38 
  - 3 : GRCm37/mm9 
  - 4 : GRCm38/mm10
- `protocol`: 
  - 1 : DNase-seq(ENCODE)
  - 2 : ATAC-seq(ENCODE) 
  - 3 : ATAC-seq(ATACdb)
- `cell_type`: refer to the function `getCelltypeList()`.
- `perbase`: 1 : Region based,2 : Per-base based.


**Example file**

The format of the chromatin regions in the input file.

```python
chr1	10732070	10733118	.	.	.
chr1	10781239	10781744	.	.	.
chr1	10795106	10799241	.	.	.
chr1	10851570	10852173	.	.	.
chr1	10965129	10966144	.	.	.
chr1	11906876	11908666	.	.	.
```


Example `task_id` and `EXAMPLE.bed` file.
```python
oaa.getExampleInputFile(save_path)

task_id=oaa.getExampleTaskID()
```
- `task_id`: The 16-bit identity of the submitted task.


**Submit**

Submit your file to server and return a `task_id` for query progress and download results.

```python
task_id=oaa.runAnnotate(input)
```

- `input`: The path of the '.bed' or '.bed.gz' file or a `list/pandas.DataFrame` format variable to be uploaded, such as `'/Users/example/example.bed'`.

**Get Result**

Get the current progress according to the `task_id`, download the result file to the local path.

```python
# You can view the annotation progress
oaa.getProgress(task_id)

# You can view the parameters you set before
oaa.viewParams(task_id)

oaa.getResultType()
'''
1 - head
2 - readopen
3 - peakopen
4 - spotopen
5 - foreread
'''

# download the annotate result
oaa.getAnnoResult(result_type, task_id ,cell_type )

# download the bed file from web server
oaa.getInputFile(save_path, task_id)
```

- `result_type`: The file type of the result, 1 - head, 2 - readopen, 3 - peakopen, 4 - spotopen, 5 - foreread.
- `save_path`: Path to save download file.
- `task_id`: The 16-bit identity of the submitted task.
- `cell_type`: You can choose one specific or more cell types in the form of `list`


Then we provide an interface `anndata`, which can embed openness data into anndata structure for downstream analysis

```python
# build ann data matrix from openness annotation result 
fromOpen2EpiScanpy(self, data, head)
```

- `data`: path to the openness result file or the output from the function `getAnnoResult()`
- `head`: path to the openness head file or the output from the function `getAnnoResult(result_type = 1)`







### Example

```python
# initial and get parameters
from OpenAnnotate import OpenAnnotateApi
oaa=OpenAnnotateApi.Annotate()
oaa.help()
oaa.getParams()
```
output:
```
Species list :
1 - GRCh37/hg19
1 - GRCh38/hg38
3 - GRCm37/mm9
4 - GRCm38/mm10
Protocol list :
1 - DNase-seq(ENCODE)
2 - ATAC-seq(ENCODE)
3 - ATAC-seq(ATACdb)
Annotate mode :
1 - Region based
2 - Per-base based
```
```python
# get example bed and task id.
# download bed file from server
task_id=oaa.getExampleTaskID()

oaa.getExampleInputFile(save_path='.')

oaa.getInputFile(save_path='.', task_id=2021061544690865)
```
output:
```
Example task id: 2020121013091517
get the result to ./EXAMPLE.bed.gz
get the result to ./2021061544690865.bed
```
Then search for the system, tissue and cell type. After setting parameters, you can submit your job to the server.

```python
oaa.getCelltypeList(protocol=1, species=1)

oaa.getTissueList(protocol=1, species=1)

oaa.getSystemList(protocol=1, species=1)

oaa.searchCelltype(protocol=1, species=1, keyword='K562')

oaa.searchTissue(protocol=1, species=1, keyword='blood')

oaa.searchSystem(protocol=1, species=1, keyword='Stem')

oaa.setParams(species=1, protocol=1, cell_type=1, perbase=1)

task_id=oaa.runAnnotate(input='./EXAMPLE.bed.gz')

# view parameters
oaa.viewParams(task_id=2021061817196919)
```
Or you can submit a bed file in list or pd.Dataframe format

```python
import pandas as pd
regions = []
with open("./EXAMPLE.bed", "r") as file:
  lines = file.readlines()
for line in lines:
  regions.append(line.split('\t'))
task_id=oaa.runAnnotate(input=regions)


pd_regions = pd.Dataframe(regions)
task_id=oaa.runAnnotate(input=pd_regions)
```



output (Omit cell type):
```
Your task id is: 2021061915336302
You can get the progress of your task through getProgress(task_id=2021061915336302)

Your task's parameters:
Protocol: DNase-seq(ENCODE)
Species: GRCh37/hg19
Cell type: All biosample types
Annotate mode: perbase based
```


```python
# download the result
oaa.getProgress(task_id=2021061817196919)
head = oaa.getAnnoResult(result_type=1, task_id=2021061817196919,cell_type=1)
```
output：
```
Your task has been completed!
You can get the result file type first through getResultType()
You can download result file through getAnnoResult(result_type, 2021061817196919)

get the result to ./head.txt.gz
```

```python
# download the result
anndata = oaa.fromOpen2EpiScanpy('./results/readopen_2021061817196919.txt', './results/head_2021061817196919.txt')
```




