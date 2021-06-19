import requests
import time
import datetime
import gzip
import os


protocolDict = {1:'dseq',2:'aseq',3:'atbd'}
speciesDict = {11:'hg19',12:'hg38',21:'mm09',22:'mm10'}
resultDict = {1:'head',2:'readopen',3:'peakopen',4:'spotopen',5:'foreread'}
perbaseDict = {1: 'perbase based', 0:'region based'}
protocolDict_w = {1:'DNase-seq(ENCODE)',2:'ATAC-seq(ENCODE)',3:'ATAC-seq(ATACdb)'}
speciesDict_w = {11:'GRCh37/hg19',12:'GRCh38/hg38',21:'GRCm37/mm9',22:'GRCm38/mm10'}
class Annotate(object):
    def __init__(self):
        super().__init__()
        self.protocol = -1
        self.species = -1
        self.celltype = -1
        self.perbase = -1
        self.data_path = -1
        self.task_id = -1
        print('Use help() to get basic functions and arguments')
    
    def help(self):
        print('getParams() : get params list')
        print('getCelltypeList(protocol, species) : get cell type list')
        print('searchCelltype(protocol, species, keyword) : search for cell types that contain keyword')
        print('setParams(assay, species, cell_type, perbase) : set params list')
        print('runAnnotate(file_path) : upload file to server')
        print('getProgress(task_id)')
        print('getAnnoResult(result_type, save_path, task_id)')
        print('getInputFile(save_path, task_id) : get your input file from server')
        print('viewParams(task_id) : view parameters')
        print('exampleTaskID() : get example task id')
        print('exampleInputFile(save_path) : get example input file to the save_path')
        

    def getParams(self):
        print('Species list : ')
        self.getSpeciesList()
        print('Protocol list : ')
        self.getProtocolList()
        print('Annotate mode : ')
        print('0 - Region based')
        print('1 - Per-base based')
    
    def setParams(self,species,protocol,cell_type,perbase):
        if protocol not in [1,2,3]:
            print('Wrong parameter! Please reset protocol')
            return
        if species not in [11,12,21,22]:
            print('Wrong parameter! Please reset species')
            return
        self.protocol = protocol
        self.species = species
        url = 'http://health.tsinghua.edu.cn/openness/anno/info/stat/celltp_%s_%s.txt'%(speciesDict[self.species],protocolDict[self.protocol])
        result = requests.get(url,stream=True)
        result = result.text
        result = result.split('\n')
        if cell_type not in list(range(1,len(result) + 1)):
           print('Wrong parameter! Please reset cell type') 
           return
        if perbase not in [0,1]:
            print('Wrong parameter! Please reset perbase')
            return 
        self.celltype = cell_type
        self.perbase = perbase
    
    def checkParams(self):
        if self.protocol == -1 or self.species == -1 or self.celltype == -1 or self.perbase == -1:
            return 0
        else:
            return 1

    def checkBED(self,data_path):
        if data_path.split('.')[-1] == 'bed':
            with open(data_path, "r") as f:
                lines = f.readlines()
        elif data_path.split('.')[-1] == 'gz':
            with gzip.open(data_path, 'rb') as f:
                lines = f.read()
                lines = lines.decode('utf-8')
                lines = lines.split('\n')
                lines = lines[:len(lines) - 1]
        else:
            return 0
        for i in range(len(lines)):
            if len(lines[i].strip().split()) != 6:
                return 0
            # 判断是否为整数
            if not lines[i].strip().split()[1].isdigit():
                return 0
            if not lines[i].strip().split()[2].isdigit():
                return 0
            if lines[i].strip().split()[0][:3] != 'chr':
                return 0
        return 1

    def getInputFile(self, save_path, task_id):
        print('get the result to %s/%s.bed'%(save_path,task_id))
        task_id = str(task_id)
        url = 'http://health.tsinghua.edu.cn/openness/anno/task/task/%s/%s.bed'%(task_id,task_id[8:])
        r = requests.get(url,stream=True)
        f = open("%s/%s.bed"%(save_path,task_id), "wb")
        for chunk in r.iter_content(chunk_size=512):
            if chunk:
                f.write(chunk)

    def runAnnotate(self,file_path):
        if not os.path.exists(file_path):
            print('ERROR! No such file or directory!')
            return
        temp_path = file_path.split('/')[-1]
        if temp_path.split('.')[1] != 'bed':
            print('ERROR! Please upload \'.bed\' or \'.bed.gz\' file')
            return
        if not self.checkBED(file_path):
            print('Wrong file format! Please check the format of your bed file.')
            return
        task_id = self.generateID()
        self.data_path = file_path
        species = self.species
        protocol = self.protocol
        celltype = self.celltype
        perbase = self.perbase

        if self.checkParams() == 0:
            print('Please set parameters first!')
            return 0

        files = {
            "species" : species,
            "protocol" : protocol,
            "celltype":celltype,
            "perbasepair":perbase,
            "taskname" : task_id
            }
        file = {"file" :open(file_path,"rb")}
        url = 'http://health.tsinghua.edu.cn/openness/anno/phpa/stepu_api.php'
        r = requests.post(url,files= file,data=files)
        self.task_id = task_id
        print('Your task id is: ' + str(task_id))
        print('You can get the progress of your task through getProgress(task_id=%s)'%(task_id))
        return task_id

    def getCelltypeList(self,protocol,species):
        if species == 12:
            species = 11
        if species == 22:
            species = 21
        if species == 21 and protocol == 3:
            print('The corresponding cell type was not found. Please reselect the parameters.')
            return
        url = 'http://health.tsinghua.edu.cn/openness/anno/info/stat/celltp_%s_%s.txt'%(speciesDict[species],protocolDict[protocol])
        result = requests.get(url,stream=True)
        result = result.text
        result = result.split('\n')
        print('1 - All biosample types')
        for i in range(len(result) - 1):
            print(str(i + 2) + " - " + result[i][8:])
    
    def searchCelltype(self,protocol,species,keyword):
        if species == 12:
            species = 11
        if species == 22:
            species = 21
        if species == 21 and protocol == 3:
            print('The corresponding cell type was not found. Please reselect the parameters.')
            return
        url = 'http://health.tsinghua.edu.cn/openness/anno/info/stat/celltp_%s_%s.txt'%(speciesDict[species],protocolDict[protocol])
        result = requests.get(url,stream=True)
        result = result.text
        result = result.split('\n')
        count = 0
        if keyword in '1 - All biosample types':
            print('1 - All biosample types')
            count += 1
        for i in range(len(result) - 1):
            if keyword in result[i][8:]:
                print(str(i + 2) + " - " + result[i][8:])
                count += 1
        if count == 0:
            print('Your keyword can not be found.')

    
    def getProtocolList(self):
        print('1 - DNase-seq(ENCODE)')
        print('2 - ATAC-seq(ENCODE)')
        print('3 - ATAC-seq(ATACdb)')

    def getSpeciesList(self):
        print('11 - GRCh37/hg19')
        print('12 - GRCh38/hg38')
        print('21 - GRCm37/mm9')
        print('22 - GRCm38/mm10')

    def getResultType(self):
        print('1 - head')
        print('2 - readopen')
        print('3 - peakopen')
        print('4 - spotopen')
        print('5 - foreread')

    # random generate task id
    def generateID(self):
        date_id = str(datetime.datetime.now().year) + str(datetime.datetime.now().month).zfill(2) + str(datetime.datetime.now().day).zfill(2)
        random_id = time.strftime("%H%M", time.localtime()) + str(int(time.time() * 10000 %10000))
        # random_id = "".join(map(lambda x:random.choice(string.digits), range(8)))
        task_id = date_id + random_id
        return task_id
    
    def viewParams(self,task_id):
        print('Your task\'s parameters:')
        url = 'http://health.tsinghua.edu.cn/openness/anno/task/task/%s/openanno.ret'%(task_id)
        r = requests.get(url,stream=True)
        task_info = r.text.split('\t')[-4:]
        task_info[3] = task_info[3].split('\n')[0]
        protocol = task_info[0]
        species = task_info[1]
        cell_type = task_info[2]
        perbase = task_info[3]
        print('Protocol: ' + protocolDict_w[int(protocol)])
        print('Species: ' + speciesDict_w[int(species)])
        if species == 12:
            species = 11
        if species == 22:
            species = 21
        if species == 21 and protocol == 3:
            print('The corresponding cell type was not found. Please reselect the parameters.')
            return
        url = 'http://health.tsinghua.edu.cn/openness/anno/info/stat/celltp_%s_%s.txt'%(speciesDict[int(species)],protocolDict[int(protocol)])
        result = requests.get(url,stream=True)
        result = result.text
        result = result.split('\n')
        if int(cell_type) ==  1:
            print('Cell type: '+ 'All biosample types')
        else:
            print('Cell type: ' + result[int(cell_type) - 2][8:])
        print('Annotate mode: '+ perbaseDict[int(perbase)] )

        

    
    
    def exampleTaskID(self):
        task_id = 2020121013091517
        print('Example task id: ' + str(task_id))
        return task_id

    def exampleInputFile(self,save_path):
        print('Get the result to %s/EXAMPLE.bed.gz'%(save_path))
        url = 'http://health.tsinghua.edu.cn/openness/anno/phpa/help/EXAMPLE.bed.gz'
        r = requests.get(url,stream=True)
        f = open("%s/EXAMPLE.bed.gz"%(save_path), "wb")
        for chunk in r.iter_content(chunk_size=512):
            if chunk:
                f.write(chunk)
        
    def getProgress(self,task_id = -1):
        if task_id == -1:
            task_id = self.task_id
        task_id = str(task_id)
        url_finish = 'http://health.tsinghua.edu.cn/openness/anno/task/%s/%s/%s/logs/openanno.ret'%(task_id[:4],task_id[4:8],task_id[8:])
        result = requests.get(url_finish,stream=True)
        if result.status_code == 404:
            url = 'http://health.tsinghua.edu.cn/openness/anno/task/%s/%s/%s/logs/openanno.sta'%(task_id[:4],task_id[4:8],task_id[8:])
            result = requests.get(url,stream=True)
            result = result.text
            print(result)
        else:
            print('Your task has been completed!')
            print('You can get the result file type first through getResultType()')
            print('You can download result file through getAnnoResult(result_type, save_path, %s)'%(task_id))


    def getAnnoResult(self,result_type,save_path,task_id = -1):
        if task_id == -1:
            task_id = self.task_id
        task_id = str(task_id)
        result_type = resultDict[result_type]
        print('Get the result to %s/%s.txt.gz'%(save_path,result_type))
        url = 'http://health.tsinghua.edu.cn/openness/anno/task/%s/%s/%s/anno/%s.txt.gz'%(task_id[:4],task_id[4:8],task_id[8:],result_type)
        r = requests.get(url,stream=True)
        f = open("%s/%s.txt.gz"%(save_path,result_type), "wb")
        for chunk in r.iter_content(chunk_size=512):
            if chunk:
                f.write(chunk)
