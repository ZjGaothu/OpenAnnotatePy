# OpenAnnotate

#### A python package for efficiently annotating the chromatin accessibility of genomic regions

Chromatin accessibility is a measure of the ability of nuclear macromolecules to physically contact DNA, and is essential for understanding regulatory mechanisms.

OpenAnnotate facilitates the chromatin accessibility annotation for massive genomic regions by allowing ultra-efficient annotation across various biosample types based on chromatin accessibility profiles accumulated in public repositories (1236 samples from ENCODE and 1493 samples from ATACdb).

For more information, please refer to the web: http://health.tsinghua.edu.cn/openannotate/

### Install OpenAnnotate via Pypi
Anaconda users can first create a new Python environment and activate it via(this is unnecessary if your Python environment is managed in other ways)

```Python
conda create python=3.6 -n OpenAnnotate
conda activate OpenAnnotate
```

OpenAnnotate is available on pypi here and can be installed via
```Python
pip install OpenAnnotate
```


### Usage  

**Import**

The package inclues a class named `OpenAnnotate`, All functions are implemented by instantiating objects of this class.
```Python
from OpenAnnotate import OpenAnnotateApi
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
getParams() : get params list
getCelltypeList(protocol,species) : get cell type list
searchCelltype(protocol, species, keyword) : search for cell types that contain keyword
setParams(assay,species,cell_type,perbase) : set params list
runAnnotate(file_path) : Upload file to server
getProgress(task_id)
getAnnoResult(result_type,save_path,task_id)
getInputFile(save_path, task_id) : get your input file from server
viewParams(task_id) : view parameters
exampleTaskID() : get example task id
exampleInputFile(save_path) : get example input file to the save_path
'''
```

**Get parameters**

Get the parameters to be set.
```python
# get basic parameters you need to set
oaa.getParams()

# get the corresponding cell type list
oaa.getCelltypeList(protocol, species)

# search cell type
oaa.searchCelltype(protocol, species, keyword)
```
- `getParams()`: Return the parameter list of `species`,`protocol` and `Annotate method`.
- `getCelltypeList(protocol,species)` : Return the cell type list of the corresponding `protocol` and `species`.
- `species` : 
  - 11 : GRCh37/hg19 
  - 12 : GRCh38/hg38 
  - 21 : GRCm37/mm9 
  - 22 : GRCm38/mm10
- `protocol`: 
  - 1 : DNase-seq(ENCODE)
  - 2 : ATAC-seq(ENCODE) 
  - 3 : ATAC-seq(ATACdb)
- `keyword`: Key word for search. Such as `K562`.


**Set parameters**

Set parameters for your object.

```python
oaa.setParams(species, protocol, cell_type, perbase)
```
- `species` : 
  - 11 : GRCh37/hg19 
  - 12 : GRCh38/hg38 
  - 21 : GRCm37/mm9 
  - 22 : GRCm38/mm10
- `protocol`: 
  - 1 : DNase-seq(ENCODE)
  - 2 : ATAC-seq(ENCODE) 
  - 3 : ATAC-seq(ATACdb)
- `cell_type`: refer to the function `getCelltypeList()`.
- `perbase`: 0 : Region based,1 : Per-base based.


**Example file**

Example `task_id` and `EXAMPLE.bed` file.
```python
oaa.exampleInputFile(save_path)

task_id=oaa.exampleTaskID()
```
- `task_id`: The 16-bit identity of the submitted task.


**Submit**

Submit your file to server and return a `task_id` for query progress and download results.

```python
task_id=oaa.runAnnotate(file_path)
```

- `file_path`: The path of the '.bed' or '.bed.gz' file to be uploaded, such as `'/Users/example/example.bed'`.

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
oaa.getAnnoResult(result_type, save_path, task_id)

# download the bed file from web server
oaa.getInputFile(save_path, task_id)
```

- `result_type`: The file type of the result, 1 - head, 2 - readopen, 3 - peakopen, 4 - spotopen, 5 - foreread.
- `save_path`: Path to save download file.
- `task_id`: The 16-bit identity of the submitted task.



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
11 - GRCh37/hg19
12 - GRCh38/hg38
21 - GRCm37/mm9
22 - GRCm38/mm10
Protocol list :
1 - DNase-seq(ENCODE)
2 - ATAC-seq(ENCODE)
3 - ATAC-seq(ATACdb)
Annotate mode :
0 - Region based
1 - Per-base based
```
```python
# get example bed and task id.
# download bed file from server
task_id=oaa.exampleTaskID()

oaa.exampleInputFile(save_path='.')

oaa.getInputFile(save_path='.', task_id=2021061544690865)
```
output:
```
Example task id: 2020121013091517
get the result to ./EXAMPLE.bed.gz
get the result to ./2021061544690865.bed
```


```python
oaa.getCelltypeList(protocol=1, species=11)

oaa.searchCelltype(protocol=1, species=11, keyword='K562')

oaa.setParams(species=11, protocol=1, cell_type=1, perbase=1)

task_id=oaa.runAnnotate(file_path='./EXAMPLE.bed.gz')

# view parameters
oaa.viewParams(task_id=2021061817196919)
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
oaa.getAnnoResult(result_type=1, save_path='.', task_id=2021061817196919)
```
output：
```
Your task has been completed!
You can get the result file type first through getResultType()
You can download result file through getAnnoResult(result_type, save_path, 2021061817196919)

get the result to ./head.txt.gz
```



