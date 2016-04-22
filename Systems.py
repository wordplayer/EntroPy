# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 01:04:42 2016

@author: Shashwat
"""
import openpyxl as xl
import numpy as np
import h5py as hdf

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
        self.x1 = np.array(data[8])
        self.x2 = 1-self.x1
        self.etaSystem = np.array(data[9])
        self.rhoSystem = np.array(data[10])
        self.temperature = data[11]
        self.reference = data[12]
        self.initialized = True #for providing a flag that all data is available
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
        massA = sheet['D3'].value
        massB = sheet['D4'].value
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
        
    
    def save(self):
        """Saves data in the HDF5 format
        
        """
        
        self.f = hdf.File("dbtest3.hdf5","a")
        dt = hdf.special_dtype(vlen=str)
        name = self.compoundA + " + " + self.compoundB + " @ " + str(self.temperature)
        subgrp = self.f.create_group(name)
        dset = subgrp.create_dataset("Specifics",(10,),dtype=dt)
        dset[0] = self.compoundA
        dset[1] = self.compoundB
        dset[2] = self.reference
        dset[3] = self.etaA
        dset[4] = self.etaB
        dset[5] = self.rhoA
        dset[6] = self.rhoB
        dset[7] = self.massA
        dset[8] = self.massB
        dset[9] = self.temperature
        
        dset = subgrp.create_dataset("Data",(3,len(self.x1)),dtype="f")
        dset[0] = self.x1
        dset[1] = self.etaSystem
        dset[2] = self.rhoSystem
        print "Saving data for " + name
        
    def showAll(self):
        """
        Prints all systems 
        """
        f = hdf.File("dbtest3.hdf5","r")
        count = 0
        for groupname in f:
            print "Name of Compound A: " + f[groupname]["Specifics"][0]
            print "Name of Compound B: " + f[groupname]["Specifics"][1]
            print "Ref: " + str(f[groupname]["Specifics"][2])
            print "Viscosity of Compound A: " + f[groupname]["Specifics"][3]
            print "Viscosity of Compound B: " + f[groupname]["Specifics"][4]
            print "Density of Compound A: " + f[groupname]["Specifics"][5]
            print "Density of Compound B: " + f[groupname]["Specifics"][6]
            print "Mass of Compound A: " + f[groupname]["Specifics"][7]
            print "Mass of Compound B: " + f[groupname]["Specifics"][8]
            print "Temperature: " + f[groupname]["Specifics"][9]
            print " "
            print "Mole Fractions: " + str(f[groupname]["Data"][0])
            print "Eta System: " + str(f[groupname]["Data"][1])
            print "Rho System: " + str(f[groupname]["Data"][2])
            print " "
            print "---------------------------------------------------------"
            print " "
            count+=1
            print " "
        print count
        
    def loadSystem(self,name):
        """
        Loads data from the database
        
        Params: Name of the system
        stored as a std way. Eg: Pentane + Heptane @ 298.15
        """
        f = hdf.File("dbtest3.hdf5","r")
        if name in f:
            self.compoundA = f[name]["Specifics"][0]
            self.compoundB = f[name]["Specifics"][1]
            self.reference = f[name]["Specifics"][2]
            self.etaA = float(f[name]["Specifics"][3])
            self.etaB = float(f[name]["Specifics"][4])
            self.rhoA = float(f[name]["Specifics"][5])
            self.rhoB = float(f[name]["Specifics"][6])
            self.massA = float(f[name]["Specifics"][7])
            self.massB = float(f[name]["Specifics"][8])
            self.temperature = float(f[name]["Specifics"][9])
            self.x1 = f[name]["Data"][0]
            self.x2 = 1-self.x1
            self.etaSystem = f[name]["Data"][1]
            self.rhoSystem = f[name]["Data"][2]
            
    def listAllSystems(self):
        """
        Lists all the systems
        """
        
        f = hdf.File("dbtest3.hdf5","r")
        for name in f:
            print name
            
    def getDBFile(self):
        """returns the database file pointer
        """
        return hdf.File("dbtest3.hdf5","a")
    
    def trial(self):
        computedEta = (self.x1**2) * (self.etaA) + (self.x2**2) * (self.etaB) + self.x1 * self.x2 * (self.etaA *(self.massB/self.massA) + self.etaB*(self.massA/self.massB))
        aapd = self.getAAPD(computedEta)
        return aapd 
        
    def trial2(self):
        computedEta= (np.exp((np.log(self.etaA/self.massA) - np.log(self.etaB/self.massB))*self.x1 + np.log(self.etaB/self.massB)) ) * (self.x1*self.massA + self.x2 * self.massB)
        return self.getAAPD(computedEta)
    
    def doBingham(self):
        """ 
        Params: None
        
        returns 
        
        """
        if self.initialized:
            computedEta = (self.x1)*(self.etaA)  +(self.x2)*(self.etaB)
            aapd=self.getAAPD(computedEta)
            print aapd
            
        else:
            print "No data available for the system."
            
            
    def doKendallMunroe(self):
    
        """
        Params: None
        
        returns void
        """
    
        if self.initialized:
            computedEta = (self.x1) * np.log(self.etaA) + (self.x2) * np.log(self.etaB)
            computedEta = np.exp(computedEta)
            aapd=self.getAAPD(computedEta)
            print aapd
    
    def doFrenkel(self):
        
        """
        For non-ideal binary mixtures
        
        Params: None
        
        returns void
        """
        if self.initialized:
            computedEta=np.exp((self.x1*self.x1*np.log(self.etaA))+(self.x2*self.x2*np.log(self.etaB))+(2*self.x1*self.x2*np.log((self.etaA+self.etaB)/2)))
            aapd=self.getAAPD(computedEta)
            print aapd
    
    def doHind(self):

       """
       Params: None

       returns void
       """
       if self.initialized:
           computedEta=(self.x1*self.x1*self.etaA)+(self.x2*self.x2*self.etaB)+(2*self.x1*self.x2*((self.etaA+self.etaB)/2))
           aapd=self.getAAPD(computedEta)
           print aapd
    
    def doEyring(self):
        
        """
        Params: None
        
        returns void
        """
        if self.initialized:
            V=((((self.x1)*(self.massA))+((self.x2)*(self.massB)))/self.rhoSystem)
            V1=self.massA/self.rhoA
            V2=self.massB/self.rhoB
            computedEta=np.exp((self.x1*np.log(self.etaA*V1))+(self.x2*np.log(self.etaB*V2)))/V
            aapd=self.getAAPD(computedEta)
            print aapd
    
    def doSW(self):
        """
        Sutherland-Wassiljewa correlation
        Params: None
        returns void
        """
        
        A11=1
        A12=np.power((1+(np.sqrt(self.etaA/self.etaB)*np.power((self.massB/self.massA),0.375))),2)/4
        A21=np.power((1+(np.sqrt(self.etaB/self.etaA)*np.power((self.massA/self.massB),0.375))),2)/4
        A22=1
        
        """
        Aij s are the Wassiljewa coefficients
        """
        
        computedEta=((self.x1*self.etaA)/((A11*self.x1)+(A21*self.x2)))+((self.x2*self.etaB)/((A12*self.x1)+(A22*self.x2)))
        aapd=self.getAAPD(computedEta)
        print aapd
    
    def doMc3b(self):
        """
        McAllister 3-body correlation
        Params: None
        returns void
        """
        
        t1=np.log((self.etaSystem/self.rhoSystem)/1000000)
        t2=self.x1*self.x1*self.x1*np.log((self.etaA/self.rhoA)/1000000)
        t3=self.x2*self.x2*self.x2*np.log((self.etaB/self.rhoB)/1000000)
        t6=np.log(self.x1+(self.x2*self.massA/self.massB))
        t7=3*self.x1*self.x1*self.x2*np.log((2+(self.massB/self.massA))/3)
        t8=3*self.x1*self.x2*self.x2*np.log((1+(2*self.massB/self.massA))/3)
        t9=self.x2*self.x2*self.x2*np.log(self.massB/self.massA)
        Y=t1-t2-t3+t6-t7-t8-t9
        X1=3*self.x1*self.x1*self.x2
        X2=3*self.x1*self.x2*self.x2
        X=np.column_stack((X1,X2))
        X=np.matrix(X)
        Y=np.matrix(Y)
        theta=np.array(np.dot(np.dot(np.linalg.inv(np.dot(X.T,X)),X.T),Y.T))
        computedEta=np.exp(t2+t3+(theta.item(0)*X1)+(theta.item(1)*X2)-t6+t7+t8+t9)*self.rhoSystem*1000000
        aapd=self.getAAPD(computedEta)
        print aapd
    
    def doGN(self):
        """
        Grunberg-Nissan correlation
        Params: None
        returns void
        """
        
        t1=np.log(self.etaSystem)
        t2=self.x1*np.log(self.etaA)
        t3=self.x2*np.log(self.etaB)
        Y=np.matrix(t1-t2-t3)
        X=np.matrix(self.x1*self.x2)
        G12=np.matrix(np.dot(np.dot(np.linalg.inv(np.dot(X,X.T)),X),Y.T))
        computedEta=np.exp(t2+t3+np.dot(G12,X))
        aapd=self.getAAPD(computedEta)
        print aapd
    
    def doRefutas(self):
        """
        Refutas correlation
        Params: None
        returns void
        """
        
        kv1 = (self.etaA/self.rhoA)*1000
        kv2 = (self.etaB/self.rhoB)*1000
        kact = (self.etaSystem/self.rhoSystem)*1000
        vbn1 = 14.534*np.log(np.log(kv1+0.8))+10.975
        vbn2 = 14.534*np.log(np.log(kv2+0.8))+10.975
        blend = self.x1*vbn1 + self.x2*vbn2
        kv = np.exp(np.exp((blend-10.975)/14.534))-0.8
        computedEta = (kv*self.rhoSystem)/1000
        aapd = self.getAAPD(computedEta)
        print aapd
        
    def doGambill(self):
        """
        Gambill correlation
        Params: None
        returns void
        """
        
        kv1 = (self.etaA/self.rhoA)*1000
        kv2 = (self.etaB/self.rhoB)*1000
        kv = self.x1*np.power(kv1,1.0/3) + self.x2*np.power(kv2,1.0/3)
        kv = np.power(kv,3)
        computedEta = (kv*self.rhoSystem)/1000
        aapd = self.getAAPD(computedEta)
        print aapd
        
        
    def getAAPD(self,computedEta):
        ad = np.abs(self.etaSystem - computedEta)
        apd = (ad/self.etaSystem)*100
        aapd = np.mean(apd)
        return aapd
        
    def save(self):
        pass
    
            
"""testing code below"""

B = BinarySystem()
data = B.loadViscosityDataFromExcel('Data.xlsx')
print data
B.create(data)
B.doKendallMunroe()
B.doBingham()
B.doFrenkel()
B.doHind()
B.doEyring()
B.doSW()
B.doMc3b()
B.doGN()
B.doRefutas()
B.doGambill()