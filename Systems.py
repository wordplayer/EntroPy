# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 01:04:42 2016

@author: Shashwat
"""
import openpyxl as xl

class BinarySystem(object):
    """
    Creates a Binary System
    
    """
    
    def __init__(self):
        
        self.initialized = False
        
    def create(self, data):
        """Assign various attributes required 
        for defining a binary system 
        """
        
        self.compoundA = data[0]
        self.compoundB = data[1]
        self.etaA = data[2]
        self.etaB = data[3]
        self.rhoA = data[4]
        self.rhoB = data[5]
        self.massA = data[6]
        self.massB = data[7]
        self.x1 = data[8]
        self.etaSystem = data[9]
        self.rhoSystem = data[10]
        self.temperature = data[11]
        self.reference = data[12]
        
        self.initialized = True
        print "Binary System created successfully."
    
    def loadViscosityDataFromExcel(self, path):
        """
        Loads Viscosity data from an Excel file
        Important to use the format of data as 
        given in the sample excel file
        
        Params
        path: path to the excel file
        
        returns
        data: an array of data
        """
        wb = xl.load_workbook(path)
        sheet = wb.get_sheet_by_name('Sheet1')
        compoundA = sheet['B1'].value
        compoundB = sheet['B2'].value
        etaA = sheet['B3'].value
        etaB = sheet['B4'].value
        rhoA = sheet['B5'].value
        rhoB = sheet['B6'].value
        temperature = sheet['D1'].value
        reference = sheet['D2'].value
        massA = sheet['D3']
        massB = sheet['D4']
        x1 = []
        etaSystem = []
        rhoSystem = []
        for i in range(1, sheet.max_row):
            x1.append(sheet.columns[4][i].value)
            etaSystem.append(sheet.columns[5][i].value)
            rhoSystem.append(sheet.columns[6][i].value)
        data = [compoundA, compoundB, etaA, etaB, rhoA, rhoB, massA, massB, x1, etaSystem, rhoSystem, temperature, reference]
        return data
    
    def __str__(self):

        """ 
        Prints the components of binary system and the temperature
        """
        
        return self.compoundA + " + " + self.compoundB + " at " + str(self.temperature) + " K"
        
"""testing code below"""

B = BinarySystem()
data = B.loadViscosityDataFromExcel('Data.xlsx')
B.create(data)