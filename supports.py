from __future__ import print_function
import os
import numpy as np
from .log import Log, ERROR, WARNING, NOTICE, DATA
from .parameters import parameters
from .delirium  import DataUtils
from .computing import get_date
from path import fpath
import time

import sys, os
try:
    import pyfits as fits
except:
    from astropy.io import fits

cm_file = os.path.join(os.path.split(__file__)[0], "data", "ControlMatrix.fits")

(CM_V, CM_H, CM_support, CM_opl) = [None]*4
def load_control_matrix():
    global CM_V, CM_H, CM_support, CM_opl
    if CM_V is None:
        cm = fits.open(cm_file)
        (CM_V, CM_H, CM_support, CM_opl) = (cm[k].data for k in ["CM_V", "CM_H", "support", "opl"])

        


class Supports(DataUtils):

    # treshold correction. Correction bellow will not a applyed 
    correction_treshold = 0.007 # [mm]
    # print a warning above this value
    maxcorwarning = 0.05   # [mm]  

    # remove the rail low orders before computing the supports
    # deformations. It appears that removing the low order even if 
    # it has been done before computing attitude helps to remove
    # some correction on the sides
    removeRailLowOrder = True 

    ## the opl range used to remove the rail low orders    
    railLowOrderRange = "good"#
    #railLowOrderRange = "conservative"

    # filter the wobble of rail before computing support 
    # deformation. Normaly it is done before on attitude carriage computing
    # the rail 
    filterRailWobble = False


    params = parameters.restrict("supports", "Hcorrection", "Vcorrection")  
    def __init__(self, rail):        
        self.num = rail.num
        if rail.carriage:
            self.filename = rail.carriage.sensors.delirium.filename
        else:
            self.filename = 'unknown'

        self.parse_period = rail.parse_period
        self.parse_range = rail.parse_range
        self.support2opl = rail.support2opl
        self.opl2support = rail.opl2support
        self.get_id = rail.get_id

        ## open a log 
        wdir, deliriumfile = os.path.split(self.filename)     
        deliriumfile, _ = os.path.splitext(deliriumfile)
        output_file = os.path.join(wdir, deliriumfile+"_CORR.txt")
        # stdout print everything 
        self.log = Log(context=("DELIRIUM", "DL%d"%self.num))
        # print only the data in the file 
        #self.log.add_output(output_file, msgtypes=DATA)
                                                  
        load_control_matrix()

        ##
        # set a correction array with all to zero 
        self.data = np.recarray(CM_support.shape,dtype=self.params.get_dtypes())
        self.data["supports"] = CM_support.astype(int)
        self.data["Hcorrection"] = np.zeros(CM_support.shape, dtype=float)
        self.data["Vcorrection"] = np.zeros(CM_support.shape, dtype=float)
        



        if rail.reconstructed:
            self.get_id = rail.get_id
            self.compute_deformations(rail)
        
    @property
    def raw_data(self):
        """ raw_data is data """
        return self.data
        
    def compute_deformations(self, rail):    
        #####
        # restrict the railarray to what is inside the correction matrix 
        # Assuming that the arrays are just truncated we do not want 
        # to make interpolation

        # clean up the railarray opl so that they are in theoritical position
        self.filename = rail.carriage.sensors.delirium.filename
        
        x = rail.get("x")
        y = rail.get("y",removeLowOrder=self.removeRailLowOrder, filterWobble=self.filterRailWobble, fitrange=self.railLowOrderRange)
        z = rail.get("z",removeLowOrder=self.removeRailLowOrder, filterWobble=self.filterRailWobble, fitrange=self.railLowOrderRange)

        ##
        # save the rail_data has there where to compute the deformation
        # for plot purpose
        self.rail_data = np.rec.fromarrays([x,y,z], names="x,y,z")
        


        self.log.notice("Computing the correction to apply on supports",3)
        delirium = rail.carriage.delirium

        xr = delirium.index2opl(delirium.opl2index(x))
        xm = delirium.index2opl(delirium.opl2index(CM_opl))
        mask = np.in1d(xr, xm)        

        self.data["Hcorrection"] =  -np.asarray(np.matrix(CM_H)*np.matrix(y[mask]).T).squeeze()
        self.data["Vcorrection"] =  -np.asarray(np.matrix(CM_V)*np.matrix(z[mask]).T).squeeze()
                       
        self.data["supports"] = CM_support.astype(int)
        

        self.maxcor = np.max( [np.max(np.abs(self.Hcorrection)), np.max(np.abs(self.Vcorrection)) ])
        l = self.log.notice if self.maxcor< self.maxcorwarning else self.log.warning
        l("Maximum correction %f micron"%(self.maxcor*1000))
    
    def get_data(self):        
        return self.data    

    def get_opl(self):
        return self.support2opl(self.supports)    

    @property
    def H(self):
        return self.data['H']

    @property
    def V(self):
        return self.data['V']

    @property
    def Hcorrection(self):
        return self.data['Hcorrection']

    @property
    def Vcorrection(self):
        return self.data['Vcorrection']

    @property
    def supports(self):
        return self.data['supports']                                        

    def apply_threshold(self):
        """ set the correction to 0 for corrections that are below the treshold 
                           
        """
        self.log.notice("Applying treshold %f to rail data "%self.correction_treshold, 3)
        treshold = self.correction_treshold
        self.Hcorrection[ np.abs(self.Hcorrection)<treshold  ] = 0.0
        self.Vcorrection[ np.abs(self.Vcorrection)<treshold  ] = 0.0
        

    def get_corrections(self, supports=None, unit="micron", treshold=None):
        """ Return the correction recarray 

        Parameters
        ----------
        supports : list of int, optional
            Return the recarray only for supports in the list.
            Otherwhise return only the supports that have at leas one direction 
            above the treshold.

        unit : string, optional    
            The unit for correction, default is micron 

        treshold : float or None
            Values bellow the treshold will be 0 
            treshold unit is always mm 
            if None use the default treshold in  self.correction_treshold   
        """
        ## check is support number exists
        ## in the CM_support 
        treshold = self.correction_treshold if treshold is None else treshold

        scale = self.get_scale(unit)
        if supports is not None:            
            supports = np.asarray(supports)
            test_valid = np.in1d(supports, self.supports)    
            if not np.all(test_valid):
                self.log.warning("Some input supports cannot be corrected because they are out of range ")
                supports = supports[test]
        
            test = np.in1d(self.supports, supports)   
        else:
            ## select all the corrections above the treshold     
            test = (np.abs(self.Hcorrection)>=treshold) + (np.abs(self.Vcorrection)>=treshold)
        if treshold>0.0:    
            Hcorrection = self.Hcorrection.copy()
            Hcorrection[ np.abs(Hcorrection)<treshold ] = 0.0
            Vcorrection = self.Vcorrection.copy()
            Vcorrection[ np.abs(Vcorrection)<treshold ] = 0.0
        else:
            Hcorrection, Vcorrection = self.Hcorrection, self.Vcorrection    
        return np.rec.fromarrays( [self.supports[test], Hcorrection[test]*scale, Vcorrection[test]*scale], names="support,H,V")

    def get_scale(self, unit):
        units_conversion_lookup = {"micron":1000, "mm":1,"m":0.001}
        try:
            scale = units_conversion_lookup[unit]
        except KeyError:
            raise ValueError("unit must be one of 'micron', 'mm', 'm' got '%s'"%unit)
        return scale        

    def get_date(self):
        """ get the date 'yyyymmdd' of the used delirium file"""        
        return get_date(self.filename)


    def txt_corrections(self, supports=None, unit="micron", treshold=None):
        treshold = self.correction_treshold if treshold is None else treshold
        corrections = self.get_corrections(supports, unit, treshold)
        order_lookup = {0:"AVG",1:"LINEAR",2:"PARABOLIC"}
        scale = self.get_scale(unit)

        txt = []
        def log(s, txt=txt):
            txt.append(s)

        log('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
        log('%% Delay Line Number = %s\n'%self.num)
        log('%% Fit along Horizontal axis = %s\n'%(order_lookup.get(parameters.get("horizontal").order,">2")))
        log('%% Fit along Vertical axis = %s\n'%(order_lookup.get(parameters.get("vertical").order, ">2")))
        log('%% Maximum correction amplitude = %f %s\n'%(self.maxcor*scale, unit))
        log('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        log('%% Original File Name = ''%s''\n'%self.filename)
        log('%% Correction Script Executed at = %s\n'%time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime()))
        log('%% all correction in %s\n'%unit);
        for c in corrections:
            log('#{c[support]:2d}  H -> {c[H]:+4.0f}                           #{c[support]:2d}  V -> {c[V]:+4.0f}\n'.format(c=c))  
        
        return "".join(txt)

    def log_corrections(self, supports=None, unit="micron", treshold=None, filename=None):
        """ get and log the correction into the log file 

        Parameters
        ----------
        supports : list of int, optional
            Return the recarray only for supports in the list.
            Otherwhise return only the supports that have at leas one direction 
            above the treshold.

        unit : string, optional    
            The unit for correction, default is micron
        treshold : float or None
            Values bellow the treshold will be 0 
            treshold unit is always mm 
            if None use the default treshold in  self.correction_treshold                  
        """
        treshold = self.correction_treshold if treshold is None else treshold
        corrections = self.get_corrections(supports, unit, treshold)
        order_lookup = {0:"AVG",1:"LINEAR",2:"PARABOLIC"}
        scale = self.get_scale(unit)
        self.log.notice("Writing correction log ")




        if filename:
            flog = fpath(filename).open("w")
            def log(s):
                print(s, end='')
                flog.write(s)
        else:                
            log = self.log.data
        log('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
        log('%% Delay Line Number = %s\n'%self.num)
        log('%% Fit along Horizontal axis = %s\n'%(order_lookup.get(parameters.get("horizontal").order,">2")))
        log('%% Fit along Vertical axis = %s\n'%(order_lookup.get(parameters.get("vertical").order, ">2")))
        log('%% Maximum correction amplitude = %f %s\n'%(self.maxcor*scale, unit))
        log('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n')
        log('%% Original File Name = ''%s''\n'%self.filename)
        log('%% Correction Script Executed at = %s\n'%time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime()))
        log('%% all correction in %s\n'%unit);
        for c in corrections:
            log('#{c[support]:2d}  H -> {c[H]:+4.0f}                           #{c[support]:2d}  V -> {c[V]:+4.0f}\n'.format(c=c))  
        
        if filename:
            flog.close()        
        