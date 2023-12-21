from itertools import count
import requests
import time
import datetime
import gzip
import os
import pandas as pd
import numpy as np
import anndata as ad
import platform
import scanpy as sc
protocolDict = {1:'dseq',2:'aseq',3:'atbd'}
speciesDict = {11:'hg19',12:'hg38',21:'mm09',22:'mm10'}
speciesDict = {1:'hg19',2:'hg38',3:'mm09',4:'mm10'}
resultDict = {1:'head',2:'readopen',3:'peakopen',4:'spotopen',5:'foreread'}
perbaseDict = {2: 'perbase based', 1:'region based'}
protocolDict_w = {1:'DNase-seq(ENCODE)',2:'ATAC-seq(ENCODE)',3:'ATAC-seq(ATACdb)'}
speciesDict_w = {1:'GRCh37/hg19',2:'GRCh38/hg38',3:'GRCm37/mm9',4:'GRCm38/mm10'}


class Annotate(object):
    def __init__(self):
        super().__init__()
        self.protocol = -1
        self.species = -1
        self.celltype = -1
        self.perbase = -1
        self.data_path = -1
        self.task_id = -1
        self.IP_addr = '159.226.47.242'
        self.port = '65533'   
        self.IP_addr_ = '166.111.5.185'
        self.port_ = '80'     
        print('Use object.help() to get basic functions and arguments')
    
    def help(self):
        print('testWebserver() : test whether the web server is working normally')
        print('setAddress(IP, port) : set the address of the web server')
        print('getParams() : get params list')
        print('getCelltypeList(protocol, species) : get cell type list')
        print('getTissueList(protocol,species) : get tissue list')
        print('getSystemList(protocol,species) : get system list')
        print('getExampleTaskID() : get example task id')
        print('getExampleInputFile(save_path) : get example input file to the save_path')
        print('searchCelltype(protocol, species, keyword) : search for cell types that contain keyword')
        print('searchTissue(protocol, species, keyword) : search for tissues that contain keyword and the corresponding cell types')
        print('searchSystem(protocol, species, keyword) : search for systems that contain keyword and the corresponding cell types')
        print('setParams(assay, species, cell_type, perbase) : set params list')
        print('runAnnotate(input) : upload file to server')
        print('getProgress(task_id) : query the annotation progress')
        print('getAnnoResult(result_type,task_id,cell_type) : download annotation result to local path')
        print('getInputFile(save_path, task_id) : get your input file from server')
        print('viewParams(task_id) : view parameters')
        print('fromOpen2EpiScanpy(data, head) : generate anndata from annotation result')
        

    def getParams(self):
        print('Species list : ')
        self.getSpeciesList()
        print('Protocol list : ')
        self.getProtocolList()
        print('Annotate mode : ')
        print('1 - Region based')
        print('2 - Per-base based')

    def SetAddress(self,IP,port):
        self.IP_addr = IP
        self.port = port

    def setParams(self,species,protocol,cell_type,perbase):
        if protocol not in [1,2,3]:
            print('Wrong parameter! Please reset protocol')
            print('Protocol list : ')
            self.getProtocolList()
            return
        if species not in [1,2,3,4]:
            print('Wrong parameter! Please reset species')
            print('Species list : ')
            self.getSpeciesList()
            return
        if species == 3 and protocol == 3:
            print('The corresponding cell type was not found. Please reselect the parameters.')
            return
        if species == 4 and protocol == 3:
            print('The corresponding cell type was not found. Please reselect the parameters.')
            return
        self.protocol = protocol
        self.species = (species//2 +1)*10 + (species%2)
        if species == 2:
            species = 1
        if species == 4:
            species = 2
        url = 'http://%s:%s/openness/anno/info/stat/celltp_%s_%s.txt'%(self.IP_addr,self.port,speciesDict[species],protocolDict[self.protocol])
        result = requests.get(url,stream=True)
        result = result.text
        result = result.split('\n')
        if isinstance(cell_type,list):
            for cell in cell_type:
                if cell not in list(range(1,len(result) + 1)):
                    print('Wrong parameter! Please reset cell type') 
                    print('Please set cell types in range 1 - %s'%(len(result)))
                return
        else:    
            if cell_type not in list(range(1,len(result) + 1)):
                print('Wrong parameter! Please reset cell type') 
                print('Please set cell types in range 1 - %s'%(len(result)))
                return
        if perbase not in [2,1]:
            print('Wrong parameter! Please reset perbase in {1, 2}')
            print('Annotate mode : ')
            print('1 - Region based')
            print('2 - Per-base based')
            return 
        self.celltype = cell_type
        self.perbase = perbase - 1
        print('Species: ', speciesDict_w[species])
        print('Protocol: ', protocolDict_w[self.protocol])
        print('Cell type: ')
        if cell_type == 1:
            print('All biosamples')
        else:
            if isinstance(cell_type,list):
                for cell in cell_type:
                    print(result[cell][8:])
            else:
                print(result[cell_type][8:])
        print('Annotate mode: ', perbaseDict[self.perbase+1])
    
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
        if os.path.exists(save_path):
            print('get the result to %s/%s.bed'%(save_path,task_id))
            task_id = str(task_id)
            url = 'http://%s:%s/openness/anno/task/task/%s/%s.bed'%(self.IP_addr,self.port, task_id,task_id[8:])
            r = requests.get(url,stream=True)
            f = open("%s/%s.bed"%(save_path,task_id), "wb")
            for chunk in r.iter_content(chunk_size=512):
                if chunk:
                    f.write(chunk)
        else:
            print('Path does not exist， please resubmit the parameters.')
    def runAnnotate(self,input):
        """ Returns the task id of annotation task

        The file_path can either be the path to the bed file, or value in 'list' format
        
        """
        file_path = input
        if isinstance(file_path,str):
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
            celltype = 1
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
            url = 'http://%s:%s/openness/anno/phpa/stepu.php'%(self.IP_addr,self.port)
            print('Uploading...')
            r = requests.post(url,files= file,data=files)
            self.task_id = task_id
            url = 'http://%s:%s/openness/anno/task/%s/%s/%s/logs/openanno.sta'%(self.IP_addr,self.port,task_id[:4],task_id[4:8],task_id[8:])
            result = requests.get(url,stream=True)
            print('Initializing...')
            while result.status_code != 200:
                time.sleep(1)
                result = requests.get(url,stream=True)
            print('Your task id is: ' + str(task_id))
            print('You can get the progress of your task through object.getProgress(task_id=%s)'%(task_id))
            return task_id
        elif isinstance(file_path ,list) or isinstance(file_path,pd.DataFrame):
            if isinstance(file_path,pd.DataFrame):
                file_path = np.array(file_path).tolist()
            if not os.path.exists('tmp'):
                os.mkdir('tmp')
            f = open('./tmp/annotatefile.bed', 'w')
            for i in range(len(file_path)):
                f.write('\t'.join(file_path[i])+'\n')
            f.close()
            file_path = './tmp/annotatefile.bed'
            task_id = self.runAnnotate(file_path)
            return(task_id)

    def getCelltypeList(self,protocol,species):
        if protocol not in [1,2,3]:
            print('Wrong parameter! Please reset protocol')
            return
        if species not in [1,2,3,4]:
            print('Wrong parameter! Please reset species')
            return
        if species == 3 and protocol == 3:
            print('The corresponding cell type was not found. Please reselect the parameters.')
            return
        if species == 4 and protocol == 3:
            print('The corresponding cell type was not found. Please reselect the parameters.')
            return
        if species == 2:
            species = 1
        if species == 4:
            species = 3
        if species == 3 and protocol == 3:
            print('The corresponding cell type was not found. Please reselect the parameters.')
            return
        url = 'http://%s:%s/openness/anno/info/stat/celltp_%s_%s.txt'%(self.IP_addr,self.port,speciesDict[species],protocolDict[protocol])
        result = requests.get(url,stream=True)
        result = result.text
        result = result.split('\n')
        print('1 - All biosample types')
        for i in range(len(result) - 1):
            print(str(i + 2) + " - " + result[i][8:])

    def getTissueList(self,protocol,species):
        if protocol not in [1,2,3]:
            print('Wrong parameter! Please reset protocol')
            return
        if species not in [1,2,3,4]:
            print('Wrong parameter! Please reset species')
            return
        if species == 3 and protocol == 3:
            print('The corresponding cell type was not found. Please reselect the parameters.')
            return
        if species == 4 and protocol == 3:
            print('The corresponding cell type was not found. Please reselect the parameters.')
            return
        if species == 2:
            species = 1
        if species == 4:
            species = 3
        if species == 3 and protocol == 3:
            print('The corresponding cell type was not found. Please reselect the parameters.')
            return
        url = 'http://%s:%s/openness/anno/info/stat/tissue_%s_%s.txt'%(self.IP_addr,self.port,speciesDict[species],protocolDict[protocol])
        result = requests.get(url,stream=True)
        result = result.text
        result = result.split('\n')
        for i in range(len(result) - 1):
            print(i+1 , '-' ,result[i][8:])

    def getSystemList(self,protocol,species):
        if species == 2:
            species = 1
        if species == 4:
            species = 3
        if species == 3 and protocol == 3:
            print('The corresponding cell type was not found. Please reselect the parameters.')
            return
        url = 'http://%s:%s/openness/anno/info/stat/system_%s_%s.txt'%(self.IP_addr,self.port,speciesDict[species],protocolDict[protocol])
        result = requests.get(url,stream=True)
        result = result.text
        result = result.split('\n')
        for i in range(len(result) - 1):
            print(i+1, '-', result[i][8:])
    
    def testWebserver(self):
        url = 'http://%s:%s/openness/anno/phpa/help/EXAMPLE.bed.gz'%(self.IP_addr,self.port)
        r = requests.get(url,stream=True)
        if r.status_code == 200:
            print('The web server is working properly!')
        else:
            print('The web server is not working')

    def searchTissue(self,protocol,species,keyword):
        if protocol not in [1,2,3]:
            print('Wrong parameter! Please reset protocol')
            return
        if species not in [1,2,3,4]:
            print('Wrong parameter! Please reset species')
            return
        if species == 3 and protocol == 3:
            print('The corresponding cell type was not found. Please reselect the parameters.')
            return
        if species == 4 and protocol == 3:
            print('The corresponding cell type was not found. Please reselect the parameters.')
            return
        if species == 2:
            species = 1
        if species == 4:
            species = 3
        if species == 3 and protocol == 3:
            print('The corresponding cell type was not found. Please reselect the parameters.')
            return
        url = 'http://%s:%s/openness/anno/info/meta/headFiles/%s_%s.txt'%(self.IP_addr,self.port,speciesDict[species],protocolDict[protocol])
        result = requests.get(url,stream=True)
        result = result.text
        result = result.split('\n')
        results = []
        for i in range(len(result)-1):
            line = result[i]
            results.append(line.split('\t')[:3])
        results = pd.DataFrame(results)
        results.columns = ['System','Tissue','Celltype']
        match_tissues = []
        for tissue in np.unique(results['Tissue'].values).tolist():
            if keyword.lower() in tissue.lower():
                match_tissues.append(tissue)
        counts = []
        for match_tissue in match_tissues:
            print('Tissue:  ' + match_tissue)
            idx = np.where(results['Tissue'].values == match_tissue)[0]
            print('Cell types: ')
            print('Index - Cell type')
            count = 0
            for cell_type in np.unique(results['Celltype'].values[idx]).tolist():
                count += 1
                idx_cell = self.searchCelltypeIndex(species,protocol,cell_type)
                print(idx_cell ,' - ' ,cell_type)
            counts.append(count)
        if np.sum(counts) == 0:
            print('Your keyword can not be found.')

    def searchSystem(self,protocol,species,keyword):
        if protocol not in [1,2,3]:
            print('Wrong parameter! Please reset protocol')
            return
        if species not in [1,2,3,4]:
            print('Wrong parameter! Please reset species')
            return
        if species == 3 and protocol == 3:
            print('The corresponding cell type was not found. Please reselect the parameters.')
            return
        if species == 4 and protocol == 3:
            print('The corresponding cell type was not found. Please reselect the parameters.')
            return
        if species == 2:
            species = 1
        if species == 4:
            species = 3
        if species == 3 and protocol == 3:
            print('The corresponding cell type was not found. Please reselect the parameters.')
            return
        url = 'http://%s:%s/openness/anno/info/meta/headFiles/%s_%s.txt'%(self.IP_addr,self.port,speciesDict[species],protocolDict[protocol])
        result = requests.get(url,stream=True)
        result = result.text
        result = result.split('\n')
        results = []
        for i in range(len(result)-1):
            line = result[i]
            results.append(line.split('\t')[:3])
        results = pd.DataFrame(results)
        results.columns = ['System','Tissue','Celltype']
        match_systems = []
        for system in np.unique(results['System'].values).tolist():
            if keyword.lower() in system.lower():
                match_systems.append(system)
        counts = []
        for match_system in match_systems:
            print('System:  ' + match_system)
            idx = np.where(results['System'].values == match_system)[0]
            print('Cell types: ')
            print('Index - Cell type')
            count = 0
            for cell_type in np.unique(results['Celltype'].values[idx]).tolist():
                count += 1
                idx_cell = self.searchCelltypeIndex(species,protocol,cell_type)
                print(idx_cell ,' - ' ,cell_type)
            counts.append(count)
        if np.sum(counts) == 0:
            print('Your keyword can not be found.')

    def searchCelltypeIndex(self,protocol,species,keyword):
        if species == 2:
            species = 1
        if species == 4:
            species = 3
        if species == 3 and protocol == 3:
            print('The corresponding cell type was not found. Please reselect the parameters.')
            return
        url = 'http://%s:%s/openness/anno/info/stat/celltp_%s_%s.txt'%(self.IP_addr,self.port,speciesDict[species],protocolDict[protocol])
        result = requests.get(url,stream=True)
        result = result.text
        result = result.split('\n')
        if keyword in '1 - All biosample types':
            return 1
        for i in range(len(result) - 1):
            if keyword.lower() in (result[i][8:]).lower():
                return(i+2)

    def searchCelltype(self,protocol,species,keyword):
        if protocol not in [1,2,3]:
            print('Wrong parameter! Please reset protocol')
            return
        if species not in [1,2,3,4]:
            print('Wrong parameter! Please reset species')
            return
        if species == 3 and protocol == 3:
            print('The corresponding cell type was not found. Please reselect the parameters.')
            return
        if species == 4 and protocol == 3:
            print('The corresponding cell type was not found. Please reselect the parameters.')
            return
        if species == 2:
            species = 1
        if species == 4:
            species = 3
        if species == 3 and protocol == 3:
            print('The corresponding cell type was not found. Please reselect the parameters.')
            return
        url = 'http://%s:%s/openness/anno/info/stat/celltp_%s_%s.txt'%(self.IP_addr,self.port,speciesDict[species],protocolDict[protocol])
        result = requests.get(url,stream=True)
        result = result.text
        result = result.split('\n')
        count = 0
        print('Cell types:')
        print('Index - Cell type')
        if keyword in '1 - All biosample types':
            print('1 - All biosample types')
            count += 1
        for i in range(len(result) - 1):
            if keyword.lower() in (result[i][8:]).lower():
                print(str(i + 2) + " - " + result[i][8:])
                count += 1
        if count == 0:
            print('Your keyword can not be found.')

    
    def getProtocolList(self):
        print('1 - DNase-seq(ENCODE)')
        print('2 - ATAC-seq(ENCODE)')
        print('3 - ATAC-seq(ATACdb)')

    def getSpeciesList(self):
        print('1 - GRCh37/hg19')
        print('2 - GRCh38/hg38')
        print('3 - GRCm37/mm9')
        print('4 - GRCm38/mm10')

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
        url = 'http://%s:%s/openness/anno/task/task/%s/openanno.ret'%(self.IP_addr,self.port,task_id)
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
        url = 'http://%s:%s/openness/anno/info/stat/celltp_%s_%s.txt'%(self.IP_addr,self.port,speciesDict[int(species)],protocolDict[int(protocol)])
        result = requests.get(url,stream=True)
        result = result.text
        result = result.split('\n')
        if int(cell_type) ==  1:
            print('Cell type: '+ 'All biosample types')
        else:
            print('Cell type: ' + result[int(cell_type) - 2][8:])
        print('Annotate Mode: '+ perbaseDict[int(perbase)+1])

    
    
    def getExampleTaskID(self):
        task_id = 2020121013091517
        print('Example task id: ' + str(task_id))
        return task_id

    def getExampleInputFile(self,save_path):
        print('Get the result to %s/EXAMPLE.bed.gz'%(save_path))
        url = 'http://%s:%s/openness/anno/phpa/help/EXAMPLE.bed.gz'%(self.IP_addr,self.port)
        r = requests.get(url,stream=True)
        f = open("%s/EXAMPLE.bed.gz"%(save_path), "wb")
        for chunk in r.iter_content(chunk_size=512):
            if chunk:
                f.write(chunk)
    # 需要修改
    def getProgress(self,task_id = -1):
        if task_id == -1:
            task_id = self.task_id
        task_id = str(task_id)
        url = 'http://%s:%s/openness/anno/task/%s/%s/%s/logs/openanno.ret'%(self.IP_addr,self.port,task_id[:4],task_id[4:8],task_id[8:])
        result = requests.get(url,stream=True)
        if result.status_code == 404:
            url = 'http://%s:%s/openness/anno/task/%s/%s/%s/logs/openanno.sta'%(self.IP_addr,self.port,task_id[:4],task_id[4:8],task_id[8:])
            result = requests.get(url,stream=True)
            result = result.text
            print(result)
        else:
            print('Your task has been completed!')
            print('You can get the result file type first through object.getResultType()')
            print('You can download result file through data = object.getAnnoResult(result_type, %s)'%(task_id))


    
    def getAnnoResult(self,result_type,task_id,cell_type):
        save_path = './results'
        if not os.path.exists('./results'):
            os.mkdir('results')
        task_id = str(task_id)
        result_type = resultDict[result_type]
        url = 'http://%s:%s/openness/anno/task/%s/%s/%s/anno/%s.txt.gz'%(self.IP_addr,self.port,task_id[:4],task_id[4:8],task_id[8:],result_type)
        r = requests.get(url,stream=True)
        f = open("%s/%s_%s_origin.txt.gz"%(save_path,result_type,task_id), "wb")
        for chunk in r.iter_content(chunk_size=1024):
            if chunk:
                f.write(chunk)
        f.close()
        print('Get the result to %s/%s_%s.txt'%(save_path,result_type,task_id))
        if platform.system() == 'Linux' or platform.system() == 'Darwin':
            os.system('gunzip %s/%s_%s_origin.txt.gz'%(save_path,result_type,task_id))
            if cell_type != 1:
                if result_type == 'head':
                    command = 'cut -f 6 %s/%s_%s_origin.txt'%(save_path,result_type,task_id)
                    command += ' > %s/%s_%s_temp.txt'%(save_path,result_type,task_id)
                    os.system(command)
                    with open("%s/%s_%s_temp.txt"%(save_path,result_type,task_id), "r") as file:
                        lines = file.readlines()
                    heads = []
                    for line in lines:
                        heads.append(line.split('\n')[0])
                    species = (self.species // 10 -1)*2  + self.species % 10
                    protocol = self.protocol
                    if species == 2:
                        species = 1
                    if species == 4:
                        species = 3
                    if species == 3 and protocol == 3:
                        print('The corresponding cell type was not found. Please reselect the parameters.')
                        return
                    url = 'http://%s:%s/openness/anno/info/stat/celltp_%s_%s.txt'%(self.IP_addr,self.port,speciesDict[species],protocolDict[protocol])
                    result = requests.get(url,stream=True)
                    result = result.text
                    result = result.split('\n')
                    cells = []
                    for i in range(len(result) - 1):
                        cells.append(result[i][8:])
                    cells = np.array(cells)
                    # filter if not all bio sample selected
                    if isinstance(cell_type,list):
                        cells = cells[[int(tmp) - 1 for tmp in cell_type]]
                    idxs = []
                    heads = np.array(heads)
                    for j in range(len(cells)):
                        idxs.extend(np.where(heads == cells[j])[0].tolist())
                    command = 'sed -n \"'
                    for idx in idxs:
                        command = command  + str(idx+1) + 'p;' 
                    command = command[:-1]
                    command += "\" %s/%s_%s_origin.txt"%(save_path,result_type,task_id)
                    command += ' > %s/%s_%s.txt'%(save_path,result_type,task_id)
                    os.system(command)
                    os.remove("%s/%s_%s_origin.txt"%(save_path,result_type,task_id))
                    os.remove("%s/%s_%s_temp.txt"%(save_path,result_type,task_id))
                    with open("%s/%s_%s.txt"%(save_path,result_type,task_id), "r") as file:
                        lines = file.readlines()
                    datas = []
                    for line in lines:
                        datas.append(line[:-1].split('\t'))
                    return datas
                else:
                    url = 'http://%s:%s/openness/anno/task/%s/%s/%s/anno/%s.txt.gz'%(self.IP_addr,self.port,task_id[:4],task_id[4:8],task_id[8:],'head')
                    r = requests.get(url,stream=True)
                    f = open("%s/%s_%s.txt.gz"%(save_path,'head',task_id), "wb")
                    for chunk in r.iter_content(chunk_size=1024):
                        if chunk:
                            f.write(chunk)
                    f.close()
                    os.system('gunzip %s/%s_%s.txt.gz'%(save_path,'head',task_id))
                    command = 'cut -f 6 %s/%s_%s.txt'%(save_path,'head',task_id)
                    command += ' > %s/%s_%s_temp.txt'%(save_path,'head',task_id)
                    os.system(command)
                    with open("%s/%s_%s_temp.txt"%(save_path,'head',task_id), "r") as file:
                        lines = file.readlines()
                    heads = []
                    for line in lines:
                        heads.append(line.split('\n')[0])
                    print('Waiting ......')
                    # get cell types
                    species = (self.species // 10 -1)*2  + self.species % 10
                    protocol = self.protocol
                    if species == 2:
                        species = 1
                    if species == 4:
                        species = 3
                    if species == 3 and protocol == 3:
                        print('The corresponding cell type was not found. Please reselect the parameters.')
                        return
                    url = 'http://%s:%s/openness/anno/info/stat/celltp_%s_%s.txt'%(self.IP_addr,self.port,speciesDict[species],protocolDict[protocol])
                    result = requests.get(url,stream=True)
                    result = result.text
                    result = result.split('\n')
                    cells = []
                    for i in range(len(result) - 1):
                        cells.append(result[i][8:])
                    cells = np.array(cells)
                    # filter if not all bio sample selected
                    if isinstance(cell_type,list):
                        cells = cells[[int(tmp) - 1 for tmp in cell_type]]
                    idxs = [1,2,3,4]
                    heads = np.array(heads)
                    for j in range(len(cells)):
                        idxs.extend((np.where(heads == cells[j])[0]+5).tolist())
                    command = 'cut -f '
                    for idx in idxs:
                        command = command  + str(idx) + ',' 
                    command = command[:-1]
                    command += " %s/%s_%s_origin.txt"%(save_path,result_type,task_id)
                    command += ' > %s/%s_%s.txt'%(save_path,result_type,task_id)
                    os.system(command)
                    os.remove("%s/%s_%s_origin.txt"%(save_path,result_type,task_id))
                    os.remove("%s/%s_%s_temp.txt"%(save_path,'head',task_id))
                    os.remove("%s/%s_%s.txt"%(save_path,'head',task_id))
                    with open("%s/%s_%s.txt"%(save_path,result_type,task_id), "r") as file:
                        lines = file.readlines()
                    datas = []
                    for line in lines:
                        datas.append(line[:-1].split('\t'))
                    return datas
            else:
                os.rename("%s/%s_%s_origin.txt"%(save_path,result_type,task_id),"%s/%s_%s.txt"%(save_path,result_type,task_id))
        else:
            print('Error! Please use Linux or MacOS')
            if result_type == 'head':
                with gzip.open("%s/%s_%s.txt.gz"%(save_path,result_type,task_id), "rt") as file:
                    lines = file.readlines()
                heads = []
                for line in lines:
                    heads.append(line.split('\t')[5])
                datas = []
                for line in lines:
                    datas.append(line.split('\t'))
                species = (self.species // 10 -1)*2  + self.species % 10
                protocol = self.protocol
                if species == 12:
                    species = 11
                if species == 22:
                    species = 21
                if species == 21 and protocol == 3:
                    print('The corresponding cell type was not found. Please reselect the parameters.')
                    return
                url = 'http://%s:%s/openness/anno/info/stat/celltp_%s_%s.txt'%(self.IP_addr,self.port,speciesDict[species],protocolDict[protocol])
                result = requests.get(url,stream=True)
                result = result.text
                result = result.split('\n')
                cells = []
                for i in range(len(result) - 1):
                    cells.append(result[i][8:])
                cells = np.array(cells)
                # filter if not all bio sample selected
                if isinstance(cell_type,list):
                    cells = cells[[int(tmp) - 1 for tmp in cell_type]]
                idxs = []
                heads = np.array(heads)
                for j in range(len(cells)):
                    idxs.extend(np.where(heads == cells[j])[0].tolist())
                datas = np.array(datas)
                datas = datas[idxs]
                datas = datas.tolist()
                f = open("%s/%s_%s.txt"%(save_path,result_type,task_id), 'w')
                os.remove("%s/%s_%s.txt.gz"%(save_path,result_type,task_id))
                for i in range(len(datas)):
                    f.write('\t'.join(datas[i]))
                f.close()
                return datas
            else:
                r = requests.get(url,stream=True)
                f = open("%s/%s_%s.txt.gz"%(save_path,'head',task_id), "wb")
                for chunk in r.iter_content(chunk_size=512):
                    if chunk:
                        f.write(chunk)
                f.close()
                with gzip.open("%s/%s_%s.txt.gz"%(save_path,'head',task_id), "rt") as file:
                    lines = file.readlines()
                heads = []
                for line in lines:
                    heads.append(line.split('\t')[5])
                with gzip.open("%s/%s_%s.txt.gz"%(save_path,result_type,task_id), "rt") as file:
                    lines = file.readlines()
                datas = []
                print('Waiting ......')
                for line in lines:
                    datas.append(line.split('\t'))
                
                # get cell types
                species = (self.species // 10 -1)*2  + self.species % 10
                protocol = self.protocol
                if species == 12:
                    species = 11
                if species == 22:
                    species = 21
                if species == 21 and protocol == 3:
                    print('The corresponding cell type was not found. Please reselect the parameters.')
                    return
                url = 'http://%s:%s/openness/anno/info/stat/celltp_%s_%s.txt'%(self.IP_addr,self.port,speciesDict[species],protocolDict[protocol])
                result = requests.get(url,stream=True)
                result = result.text
                result = result.split('\n')
                cells = []
                for i in range(len(result) - 1):
                    cells.append(result[i][8:])
                cells = np.array(cells)
                # filter if not all bio sample selected
                if isinstance(cell_type,list):
                    cells = cells[[int(tmp) - 1 for tmp in cell_type]]
                idxs = [0,1,2,3]
                heads = np.array(heads)
                for j in range(len(cells)):
                    idxs.extend((np.where(heads == cells[j])[0]+4).tolist())
                datas = np.array(datas)
                datas = datas[:,idxs]
                datas = datas.tolist()
                f = open("%s/%s_%s.txt"%(save_path,result_type,task_id), 'w')
                os.remove("%s/%s_%s.txt.gz"%(save_path,result_type,task_id))
                for i in range(len(datas)):
                    f.write('\t'.join(datas[i])+'\n')
                f.close()
                return datas
            
    # 修改三
    def fromOpen2EpiScanpy(self, data_path, head_path):
        """
        build ann data matrix from openness annotation result    
        """
        # load openness data
        datas = []
        heads = []
        if isinstance(data_path,str):
            with open(data_path, "r") as file:
                lines = file.readlines()
            for line in lines:
                ls_line = line.split('\t')
                ls_line[-1] = ls_line[-1].split('\n')[0]
                datas.append(ls_line[4:])
                heads.append('-'.join(ls_line[:3]))
        elif isinstance(data_path,list):
            for ls_line in data_path:
                ls_line[-1] = ls_line[-1].split('\n')[0]
                datas.append(ls_line[4:])
                heads.append('-'.join(ls_line[:3]))
        else:
            print('Error: Please submit parameters in the correct format')
            return 0
        datas = np.array(datas).T
        pd_data = pd.DataFrame(data=datas,columns = heads)
        ls_path = data_path.split('.')
        save_path = '.'.join(ls_path[:len(ls_path)-1]) + '.csv'
        pd_data.to_csv(save_path)
        anndata = sc.read_csv(save_path)
        # header
        headers = []
        if isinstance(head_path,str):
            with open(head_path, "r") as file:
                lines = file.readlines()
            for line in lines:
                ls_line = line.split('\t')
                ls_line[-1] = ls_line[-1].split('\n')[0]
                headers.append(ls_line)    
        elif isinstance(head_path,list):
            for ls_line in lines:
                ls_line[-1] = ls_line[-1].split('\n')[0]
                headers.append(ls_line)    
        else:
            print('Error: Please submit parameters in the correct format')
            return 0
        anndata.obs['biosample'] = np.array(headers)[:,5]
        return anndata
        
        
        
        
        
            
            
            


