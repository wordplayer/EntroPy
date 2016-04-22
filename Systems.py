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
        
        computedEta=((self.x1*self.etaA)/((A11*self.x1)+(A12*self.x2)))+((self.x2*self.etaB)/((A21*self.x1)+(A22*self.x2)))
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
        #theta=np.array(np.dot(np.dot(np.linalg.inv(np.dot(X.T,X)),X.T),Y.T))
        theta=self.getTheta(X,Y)
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
        X=np.matrix(self.x1*self.x2).T
        #G12=np.matrix(np.dot(np.dot(np.linalg.inv(np.dot(X,X.T)),X),Y.T))
        G12=self.getTheta(X,Y)
        computedEta=np.exp(t2+t3+X.T*G12.item(0))
        aapd=self.getAAPD(computedEta)
        print aapd
    
    def doWijk(self):
        """
        Wijk correlation
        Params: None
        returns void
        """
        
        t1=np.log(self.etaSystem)
        t2=self.x1*self.x1*np.log(self.etaA)
        t3=self.x2*self.x2*np.log(self.etaB)
        Y=np.matrix(t1-t2-t3)
        X=np.matrix(2*self.x1*self.x2).T
        #theta=np.matrix(np.dot(np.dot(np.linalg.inv(np.dot(X,X.T)),X),Y.T))
        theta=self.getTheta(X,Y)
        computedEta=np.exp(t2+t3+X.T*theta.item(0))
        aapd=self.getAAPD(computedEta)
        print aapd
    
    def doMc4b(self):
        """
        McAllister 4-body correlation
        Params: None
        returns void
        """
        
        t1=np.log((self.etaSystem/self.rhoSystem)/1000000)
        t2=np.power(self.x1,4)*np.log((self.etaA/self.rhoA)/1000000)
        t6=np.power(self.x2,4)*np.log((self.etaB/self.rhoB)/1000000)
        t7=np.log(self.x1+(self.x2*self.massB/self.massA))
        t8=4*np.power(self.x1,3)*self.x2*np.log((3+(self.massB/self.massA))/4)
        t9=6*np.power(self.x1,2)*np.power(self.x2,2)*np.log((1+(self.massB/self.massA))/2)
        t10=4*self.x1*np.power(self.x2,3)*np.log((1+(3*self.massB/self.massA))/4)
        t11=np.power(self.x2,4)*np.log(self.massB/self.massA)
        Y=np.matrix(t1-t2-t6+t7-t8-t9-t10-t11)
        X1=np.matrix(4*(np.power(self.x1,3)*self.x2))
        X2=np.matrix(6*(np.power(self.x1,2)*np.power(self.x2,2)))
        X3=np.matrix(4*(self.x1*np.power(self.x2,3)))
        X=np.matrix(np.hstack((X1.T,X2.T,X3.T)))
        #theta=np.array((np.dot(np.dot(np.linalg.inv(np.dot(X.T,X)),X.T),Y.T)))
        #theta=np.array(np.linalg.inv(X.T*X)*X.T*Y.T)
        theta=self.getTheta(X,Y)
        computedEta=1000000*np.matrix(self.rhoSystem).T*np.matrix(np.exp(t2+(theta.item(0)*X1)+(theta.item(1)*X2)+(theta.item(2)*X3)+t6-t7+t8+t9+t10+t11))
        aapd=self.getAAPD(computedEta)
        print aapd
    
    def doKC(self):
        """
        Katti-Chaudhri correlation
        Params:None
        returns void
        """
        V=((self.x1*self.massA)+(self.x2*self.massB))/self.rhoSystem
        V1=self.massA/self.rhoA
        V2=self.massB/self.rhoB
        t1=np.log(self.etaSystem*V)
        t2=self.x1*np.log(self.etaA*V1)
        t3=self.x2*np.log(self.etaB*V2)
        Y=np.matrix(t1-t2-t3)
        X=np.matrix(self.x1*self.x2).T
        #A12=np.matrix(np.linalg.inv(X.T*X)*X.T*Y.T)
        A12=self.getTheta(X,Y)
        computedEta=np.exp((A12.item(0)*X.T)+t2+t3)/V
        aapd=self.getAAPD(computedEta)
        print aapd
    """
    def doMc4b1(self):
        
        #McAllister 4-body correlation
        #Params:None
        #returns void
        
        t1=np.log((self.etaSystem/self.rhoSystem)/1000000)
        t2=np.power(self.x1,4)*np.log((self.etaA/self.rhoA)/1000000)
        t6=np.power(self.x2,4)*np.log((self.etaB/self.rhoB)/1000000)
        t7=np.log(self.x1+(self.x2*self.massB/self.massA))
        t8=4*np.power(self.x1,3)*self.x2*np.log((3+(self.massB/self.massA))/4)
        t9=6*np.power(self.x1,2)*np.power(self.x2,2)*np.log((1+(self.massB/self.massA))/2)
        t10=4*self.x1*np.power(self.x2,3)*np.log((1+(3*self.massB/self.massA))/4)
        t11=np.power(self.x2,4)*np.log(self.massB/self.massA)
        Y=np.matrix(t1-t2-t6+t7-t8-t9-t10-t11)
        X1=np.matrix(4*np.power(self.x1,3)*self.x2)
        X2=np.matrix(6*np.power(self.x1,2)*np.power(self.x2,2))
        X3=np.matrix(4*self.x1*np.power(self.x2,3))
        X=np.matrix(np.hstack((X1.T,X2.T,X3.T)))
        theta=np.array(np.linalg.inv(X.T*X)*X.T*Y.T)
        computedEta=1000000*np.matrix(self.rhoSystem).T*np.matrix(np.exp(t2+(X1*theta.item(0))+(X2*theta.item(1))+(X3*theta.item(2))+t6-t7+t8+t9+t10+t11))
        aapd=self.getAAPD(computedEta)
        print aapd
    """ 
    
    def doTK(self):
        """
        Tamura-Kurata correlation
        Params: None
        returns void
        """
        #V=((self.x1*self.massA)+(self.x2*self.massB))/self.rhoSystem
        #phi1=(self.x1*self.massA/self.rhoA)/V
        #phi2=(self.x2*self.massB/self.rhoB)/V
        V1=self.x1*self.massA/self.rhoA
        V2=self.x2*self.massB/self.rhoB
        phi1=V1/(V1+V2)
        phi2=1-phi1
        t1=self.etaSystem
        t2=self.x1*phi1*self.etaA
        t3=self.x2*phi2*self.etaB
        Y=np.matrix(t1-t2-t3)
        X=np.matrix(2*np.sqrt(self.x1*self.x2*phi1*phi2)).T
        #T12=(np.linalg.inv(X.T*X)*X.T*Y.T)
        T12=self.getTheta(X,Y)
        computedEta=t2+t3+T12.item(0)*X.T
        aapd=self.getAAPD(computedEta)
        print aapd
    
    def getTheta(self,X,Y):
        theta=np.array(np.linalg.inv(X.T*X)*X.T*Y.T)
        return theta
        
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
B.create(data)
B.doKendallMunroe()
B.doBingham()
B.doFrenkel()
B.doHind()
B.doEyring()
B.doSW()
B.doMc3b()
B.doGN()
B.doWijk()
B.doMc4b()
B.doKC()
#B.doMc4b1()
B.doTK()