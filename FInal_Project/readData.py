import csv
import numpy as np
class CSVData:

    def __init__(self,fileName):
        self.fileName = fileName
        with open(self.fileName, mode = 'r') as file:
            self.csvFile = csv.reader(file)
            self.data = np.loadtxt(self.csvFile, delimter=",",skip_header = 1)

    def getRows(self):
        return [row for row in self.data]
    def getColumns(self):
        return [[self.data[:, i] for i in range(self.data.shape[1])]]
    
    def updateCSV(self,fileName,data):
        np.savetxt(fileName, data, delimiter=",")
        
