
from .log import Log
from .io import _path2dl

from .computing import (rad2arcsec, get_date, remove_low_order, recarray_matrix_convertion, range2mask)
from .wobblefit import WobbleFitPol4, WobbleFitSin
from .parameters import parameters
from path import fpath

import numpy as np 

try:
    unicode
except NameError:
    basestring = (str, bytes)



class DataError(ValueError):
    pass

class DataUtils(object):
    """ This is a buch of functions and setup for all object with data.
        

    """
    ## Set the fitter class for wobble correction 
    #wobblefit = WobbleFitPol4()
    wobblefit = WobbleFitSin()
        
    # The opl range for the fit of low order
    # removal. See range_lookup in dl.py 
    opl_low_order_fit_range = "good"
        
    # The opl range used to remove the low order before
    # fitting the wobble
    opl_wobble_range = "good"



    def get_id(self):
        """generic get_id for unknown delirium overwriten 
           after a delirium has been set. 
        """
        return (getattr(self, "num",0), "UNKNOWN-%s"%id(self), False)


    def filter_wobble(self, x, y, period, order, removeLowOrder=True, fitrange=None):
        """ filter the wobble of a array value

        Parameters
        ----------
        x : array like
            x=opl data
        y : array like,
            y measurement 
        period : float, 
            wobble period in opl 
        order : int,
            the N low order will be removed before wobble fitting.
        
        removeLowOrder : boolean, optional
            default is True.
        fitrange : 2 tuple or string
            min, max tuple (in opl) for the fit                 
            if string must be a defined range (see dl.py)    

        Outputs
        -------
        residuals : array
            y without wobble
        amplitude : float
            wobble amplitude
        fit : WobbleFit object
            the fit object used to performe the fit
        range : 2 tuple
            min, max tuple range used    
        """
        # Position of DELIRUM measurements

        ## convert string range to numerical if needed
        rg = self.parse_range(self.opl_wobble_range if fitrange is None else fitrange)
        mask = range2mask(x,rg)        

        if removeLowOrder:
            residuals, fitpol = remove_low_order(x, y, order, mask, xrange=rg)
        else:
            residuals, fitpol = y, []

        phi   = np.mod(x, period)*2*np.pi/period
        fit = self.wobblefit.new(phi[mask], residuals[mask])                    
        # perform the fit
        fit.fit()

                            
        angle = np.linspace(0, 2*np.pi, 1000)
        model = fit.model(angle)                    
        # Lazy way to find amplitude     
        amplitude = max(model)-min(model)
                        
        return y-fit.model(phi), amplitude, fit, rg  

    def remove_low_order_of_param(self, param, x, value, order=None, fitrange=None):
        """ Remove the low order of a parameter data

        Parameters
        ----------
        param :  string or Parameter
            define one of the delirium parameter 'yctr', 'zctr', etc 
            param is used to get the default fit order from them
            
        """
        if isinstance(param, basestring):
            param = self.params.get(param)
        
        rg = self.opl_low_order_fit_range if fitrange is None else fitrange
        rg = self.parse_range(rg)            
        order = param.order if order is None else order
        if order is None:
            return value, [], rg 

        return remove_low_order(x, value, order, xrange=rg)+(rg,)

    def iopl(self, rge):
        """ Return index of values delimited by min.max range

        Parameters
        ----------
        range : 2 tuple or scalar float
            If a 2 tuple return a list of data index that contain opl between 
            the min,max range.
            If scalar return the index of the closest opl to the given value

        Outputs
        -------
        indexes : integer array or scalar                

        """
        opl = self.get_opl()
        if hasattr(rge, "__iter__"):
            oplmin, oplmax = rge
            if oplmin is None and oplmax is None:
                return np.arange(len(opl))
            if oplmax is None:
                return np.where(opl>=oplmin)[0]
            if oplmin is None:
                return np.where(opl<=oplmax)[0]
            return np.where( (opl<=oplmax)*(opl>=oplmin))[0]            
        ##
        # if a scalr return the closest index in opl
        return np.abs(opl-rge).argmin()

    data_aliases = None    
    def get(self, key, removeLowOrder=None, filterWobble=None, 
            arcsec=False, period=None, order=None, fitrange=None,  
            extras=False, inside="recarray", raw=False, 
            indexes=None, opls=None, wobbleFitrange=None, wobbleOrder=None
            ):
        """ get a measured or computed value 

        Parameters
        ----------
        key : string or  list of string
            the name(s) defining the data parameter

            if key is a list of key the values are returned in a recarray by default but can be change with the *inside* parameter
        
        removeLowOrder: bool, optional
            Remove the lowest order polynome. The order is pre-defined for each parameters (see .params attribute)
            Of course, has no effect if key is 'opl' or 'x'

        filterWobble: bool, optional
            If True, the wobble is filtered. The wobble period is pre-defined for each measurement (see .params attribute)
            This has no effect on 'opl'

        arcsec : bool, optional
            default is False. If True, convert the angle in arcsec, has only effect on  'phi', 'theta' and 'psi'             

        period : string or float or None, optional
            if not None alter the default period of requested params for wobble fit.
            Period can be a string representing a physical thing : 'wheel', 'support'
            or a float in opl [m] unit.
            This has no effect if filterWobble is False     

        order : float or None
            if not None alter the default low order for the requested param
            This has no effect if removeLowOrder is False 
            (Note that the wobble removel if any will still use the default order)

        fitrange : string, 2xtuple or None
            if not None change the default range for the low order fit.
            Can be a (min, max) tuple in opl [m]
            Or a string that define a range for each delay line : 'good', 'full', 'conservative'
            This has no effect if removeLowOrder is False 

        wobbleFitrange :  string, 2xtuple or None
            if not change the default range for the wobble function when 
            removing the low order   

        wobbleOrder : int or None
            if not None, change the default order to fit the low order 
            before fitting the wobble         

        extras : bool, optional
            default is False. If True, return a dictionary with more information 
            about value computation. What is inside the extra dictionary depend 
            on option: 
                if filterWobble is True: 
                    wobble_amplitude : amplitude value 
                    wobble_fit :  the fit object result            

        inside : string, optional
            Define the output type is key is a list of key 
            "recarray" (default), "array", "matrix", "list" or "dict"
            Has no effect if key is a string.

        raw : boolean, optional
            if true return the corresponding raw data if available
            The raw data depend on which object is called.
            For a carriage : raw data is the data without the low order removed
            For a sensors : raw data is the data without WPL roll error
            For other object : raw data is data 

        indexes : array like of int
            Indexes of returned values 
            is indexes is not None, opls must be None

        opls : 2 tuple or float
            if float return the value for the closest opl
            if 2 tuple (min,max) return the values in between or 
            equal min and max.
            if opls is not None, indexes must be None           

        Outputs
        -------
        data : recarray, array, matrix, list or dict 
            see inside parameter.

        """
        #if self.data is None and self.data_aliases is None:
        #    raise ValueError("There is no data, probably not associated to a dlirium file")

        if indexes is not None:
            if opls is not None:
                raise ValueError("cannot set 'opls' when indexes is set")    

        if hasattr(key, "__iter__") and not isinstance(key, basestring):
            if opls is not None:
                indexes = self.iopl(opls)

            datalist = [self.get(k, removeLowOrder=removeLowOrder, 
                                    filterWobble=filterWobble, raw=raw, 
                                    fitrange=fitrange, period=period, 
                                    arcsec=arcsec, order=order,
                                    extras=extras, indexes=indexes) for k in key]
            if extras:
                datalist, extras = zip(*datalist)

            if inside == "recarray":
                array_out =  np.rec.fromarrays(datalist, names=key)
            elif inside == "array":
                array_out =  np.array(datalist)
            elif inside == "matrix":
                array_out =  np.matrix(datalist)
            elif inside == "list":
                array_out =  datalist
            elif inside == "dict":
                array_out =  dict(zip(key,datalist))    
            else:
                raise ValueError("inside should be one of 'recarray', 'array', 'matrix', 'list' or 'dict' got '%s'"%inside)    
            if extras:
                return array_out, extras
            return array_out        
        ############################################
        # scalar case
        ############################################    
        
        
        ##################################
        # handle aliases
        ################################
        if self.data_aliases is not None:
            aliases = dict(self.data_aliases)
            if key in aliases:
                return self.get(aliases[key], removeLowOrder=removeLowOrder, 
                                filterWobble=filterWobble, raw=raw, 
                                fitrange=fitrange, period=period, 
                                arcsec=arcsec, order=order,
                                extras=extras, indexes=indexes)

        ###############################################
        # handle case when key is, e.g. carriage.theta
        ###############################################            
                            
        container = self
        left = key
        while left:
            attr, _, left = left.partition(".")
            if left:
                container = getattr(container, attr)
        if container is not self:
            key = attr    
            return container.get(key, removeLowOrder=removeLowOrder, 
                                    filterWobble=filterWobble, raw=raw, 
                                    fitrange=fitrange, period=period, 
                                    arcsec=arcsec, order=order,
                                    extras=extras, indexes=indexes)                                    


        if getattr(self, "data", None) is None:
            raise ValueError("missing data")

        if opls is not None:
            indexes = self.iopl(opls)            

        extras_dict = {"param":key, 
                       "id":getattr(self, "get_id", lambda : (0, 'unknown',False))(), 
                       "wobble_filtered":False,
                       "loworder_removed": False
                    }

                
        if key == "opl" and hasattr(self, "get_opl"):
            data =  self.get_opl()

        else:
            try:
                param = self.params.get(key)    
            except:
                KeyError("unknown data of name '%s'"%key)                
            data = self.raw_data[key] if raw else self.data[key]

        param = self.params.get(key)
        
        extras_dict.update(unit = param.unit)                
        if arcsec:            
            if param.unit == "rad":
                data = rad2arcsec(data)
                extras_dict.update(unit="arcsec")

        if filterWobble:
            opl = self.get_opl()
            period = self.parse_period(param.period if period is None else period)                          
            worder = param.order if wobbleOrder is None else wobbleOrder

            if period is not None:            
                newdata, amplitude, fit, rge = self.filter_wobble(opl, data, period, worder, fitrange=wobbleFitrange)
                extras_dict.update( wobble_amplitude=amplitude, wobble_fit=fit, wobble_filtered=True, wobble_range=rge)
                data = newdata



        if removeLowOrder:
            opl = self.get_opl()

            data, pol, true_fitrange = self.remove_low_order_of_param(param,opl,data, order=order, fitrange=fitrange)
            extras_dict.update( loworder_polynome = pol, loworder_range = true_fitrange, loworder_removed=True)
            # if string convert the range to numerical
            # remove the distance between            

        if indexes is not None:
            data = data[indexes]    
            
        if extras:
            return data, extras_dict        
        return data

    def _get_extra(self, extraname,  key=None, period=None, arcsec=False):
        if key is None:
            key = [p.name for p in self.params]
        if hasattr(key, "__iter__"):
            return dict( (k,self._get_extra(k, extraname, period=period)) for k in key)

        _, extras = self.get(key, filterWobble=True, period=period, extras=True, arcsec=arcsec)    
        return extras[extraname]
    
    def get_wobble_amplitude(self, key=None, period=None, arcsec=False):
        """ return the wobble amplitude of object value 

        Parameters
        ----------
        key : string
            name of the measurement/parameter

        period : float or string
            opl period of the wobble fitting. If string must be a valid period
            definition like 'wheel' or 'support'
        
        arcsec : bool
            If the parameter is an angle and arcsec is True, return amplitude in 
            arcsec instead of rad
        Outputs
        -------
        amplitude : float
            the wobble amplitude in wathever unit the measurement is

        """
        return self._get_extra("wobble_amplitude", key, period=period, arcsec=arcsec)        

    def get_wobble_fit(self, key=None, period=None):
        """ Return wobble fit object 

        Parameters
        ----------
        key : string
            name of the measurement/parameter

        period : float or string
            opl period of the wobble fitting. If string must be a valid period
            definition like 'wheel' or 'support'
        
        Outputs
        -------
        fit : WobbleFit object
            the object used to perform the fit.
        """
        return self._get_extra("wobble_fit", key, period=period)        

    @staticmethod    
    def parse_period(period):
        """ function to parse period, can be altered from the parent DL class """
        return float(period) if period not in [False,None] else period

    @staticmethod    
    def parse_range(range):
        """ function to parse ranges, can be altered from the parent DL class """
        if range is None:
            range = (None, None)
        try:
            _,_ = range
        except ValueError:
            raise ValueError("range must be a 2 tuple, got '%s'"%s)        
        return range
                    
    @property
    def raw_data(self):
        """ raw data is raw """
        raw_data = getattr(self, "_raw_data", None)
        if raw_data is None:
            return self.data
        return raw_data


class Delirium(DataUtils):
    ## parameters of the D (data) recarray 
    params = parameters.restrict( "time", "opl", "doms", "incl",
                                  "yctr", "zctr", "yend", "zend"                                
                                )

    # DELIRIUM spatial sampling (in m)
    opl_sampling = 0.375 # [m]
    # opl normaly start at 
    opl_offset   = 0.75  # 1.125-0.375 [m]
    # min and max for data completness check 
    opl_min = 1.125     # [m]
    opl_max = 118.875   # [m]
    # spatial sampling tolerance for check
    opl_tol  = 0.03  # [m]

    # define the check to do when loading data 
    data_check = {
        "complete":True, # check if scan complete add data if needed
        "sample": True, # check scan step of 0.375m 
        "nan": True, # check nan values 
        "glitches": True
    }

    # data is filled up by load_data
    data = None
    # flag to check if table is reverse or not 
    reverse = False

    ## list of var name / string to match
    ## in order to read the header
    header_lookup = {
        "tunnel_temp" : "Tunnel Temperature="
    }
    header_txt = []
    header = {}
    def __init__(self, sensors, filepath=None):

        self.dlnum = sensors.num

        self.parse_period = sensors.parse_period
        self.parse_range = sensors.parse_range

        self.f = None
        self.fpath = fpath(filepath) if filepath else None
        self.filepath = filepath                
        ## log sent to stdout by default 
        self.log = Log(context=("DELIRIUM", "DL%d"%sensors.num))
        #####
        # One can add a logfile:
        ## self.log.add_output( "some/path/to/logfile" )

    def reload(self):
        self.f = None
        self.fpath = fpath(self.filepath) if self.filepath else None
        self.data = None    
        self.load_data()
   
    def get_id(self):
        """ return a unique unique tuple identification for the dlstate """    
        return (self.dlnum, self.get_date(), self.reverse)
        #return "DL%d-%s"%(self.dlnum, self.get_date())
        
    def get_date(self):
        """ get the date 'yyyymmdd' of the used delirium file"""        
        if self.f:        
            return get_date(self.filename)
        else:
            return "SIMU%d"%id(self)    

    def get_date2(self):
        """ get the date 'yyyy-mm-dd' of the used delirium file"""        
        if self.f:        
            date = get_date(self.filename)
            return date[0:4] + "-" + date[4:6] + "-" + date[6:8]
        else:
            return "SIMU%d"%id(self)    
            

    def index2opl(self, index):
        """ convert measurement index to the theoritical opl

        Parameters
        ----------
        index : array like
            measurement index (starting from 0)

        Ouputs
        ------      
        opl : array like
            the opl distance [m]
        """
        return self.opl_offset+(np.asarray(index)+1)*self.opl_sampling


    def opl2index(self, opl):
        """ convert opl distance to the closest measurement index

        Parameters
        ----------
        opl : array like
            the opl distance [m]

        Ouputs
        ------      
        index : array like
            measurement index (starting from 0), array of int

        Notes
        -----
        There is no check of physical values boundaries 
        """
        # :TODO: check opl min and max 
        opl = np.asarray(opl)
        index = (np.round((opl-self.opl_offset) / self.opl_sampling))-1
        return index.astype(int)

    # def support2index(self, support):
    #     """ convert support numbers to its closest data measurement index

    #     Parameters
    #     ----------
    #     support : array like
    #         upport number (starting from 1)

    #     Outputs
    #     -------
    #     index : array like
    #         measurement index (starting from 0), array of int
                    
    #     Notes
    #     -----
    #     Support number start from 1
    #     There is no check of physical values boundaries   
    #     """

    #     #if isinstance(offset, basestring):
    #     #    if offset == "ctr":
    #     #        offset = -self.Fctr[0]/1000.*2
    #     #    elif offset == "end":
    #     #        offset = -self.Fend[0]/1000.*2    
    #     #    else:
    #     #        raise ValueError("if string, offset must be 'ctr' or 'end' got '%s'"%offset)    
    #     offset = 0.0        
    #     opl = self.rail.support2opl(support)+offset

    #     return self.opl2index(opl)

    # def index2support(self, index):
    #     """ convert data measurement index to the closest support number

    #     Parameters
    #     ----------
    #     index : array like
    #         measurement index
        

    #     Outputs
    #     -------
    #     support : array like
    #         support number (starting from 1)
                    
    #     Notes
    #     -----
    #     Support number start from 1
    #     There is no check of physical values boundaries   
    #     """
    #     #if isinstance(offset, basestring):
    #     #    if offset == "ctr":
    #     #        offset = -self.carriage.Fctr[0]/1000.*2
    #     #    elif offset == "end":
    #     #        offset = -self.carriage.Fend[0]/1000.*2    
    #     #    else:
    #     #        raise ValueError("if string, offset must be 'ctr' or 'end' got '%s'"%offset)    
    #     offset = 0.0        
    #     opl = self.index2opl(index)+offset
    #     return self.rail.opl2support(opl)

    
    @property
    def filename(self):
        return self.f.name if self.f else "SIMU"
        

    def load_header(self):
        self.header_txt = []
        self.header = dict((k,-999.99) for k in self.header_lookup.keys())
        if not self.fpath:      
            return 

        if not self.f:
            self.f = self.fpath.open()

        while True:
            offset = self.f.tell()
            line = self.f.readline().strip()
            if not line:
                break
            if not line[0:1]=="%":
                # put the file backward
                self.f.seek(offset)
                break

            self.header_txt.append(line)
            for var, search in  self.header_lookup.items():
                where = line.find(search)
                if where>-1:
                    line =  line[where+len(search):]
                    sval = line.split(" ")[0]
                    try:
                        val = float( sval  )
                    except (TypeError, ValueError):
                        self.log("Warning cannot read '%s' header parameter"%var)    
                    self.header[var] = val                            

    def load_data(self):
        """ Load the data from delirium file and make some check on it 

        The check are set thanks to the data_check dictionary attribute:
            'complete' : check if all data is there
            'sample' : check if the data is well sampled
            'nan' :  check for Nan Values 
            'glitches' :  check glitches 

        """        
        ## load the raw data 
        self.load_header()
        if self.f:
            self.log.notice("Loadding raw data of '%s'"%self.filename, 3)            
            print("self.f", self.f)
            data = np.loadtxt(self.f, comments="%", dtype=self.params.get_dtypes())
        else:
            N = int( (self.opl_max-self.opl_min)/self.opl_sampling+1)
            data = np.zeros((N,), dtype=self.params.get_dtypes())
            data['opl'] = np.linspace(self.opl_min, self.opl_max, N)
            self.data = data
            self.log.notice("Fake Delirium created with empy value",3)
            return 

        if not len(data):
            self.log.error("data is empty")
            raise DataError("data is empty")    
        ## set the reverse flag for future isterezis reference
        self.reverse = data["opl"][0]>data["opl"][-1]
        ## now sort the table by opl 
        data.sort(order="opl")

        opl = data["opl"]
        opl_sampling = self.opl_sampling

        ## complete the data to nominal range 
        ## so that the table as always the good size
        if self.data_check["complete"]:            
            self.log.notice("Check data completness", 3)    
            first = self.opl2index(opl[0])

            addrowfirst = 0
            addrowend   = 0
            if first > 0:            
                self.log.warning("Missing data at begining of scan, data start at OPL =%.2f"%opl[0])
                addrowfirst += first
            last = self.opl2index(opl[-1])
            ## the maximum index from the opl_max distance 
            imax = self.opl2index(self.opl_max)
            if last < imax:
                self.log.warning("Missing data at end  of scan, data end at OPL =%.2f"%opl[-1])             
                addrowend += int(imax-last)

            if addrowfirst+addrowend:                                
                newdata = np.recarray( (data.shape[0]+addrowfirst+addrowend, ), dtype=data.dtype)

                if addrowend:
                    newdata[addrowfirst:-addrowend] = data    
                else:
                    newdata[addrowfirst:] = data                                            
                newdata[0:addrowfirst] = data[0]
                newdata["opl"][0:addrowfirst] = self.index2opl(np.arange(first))

                if addrowend:
                    newdata[-addrowend:] = data[-1]
                    newdata["opl"][-addrowend:] = self.index2opl(np.arange(last,imax))

                data = newdata
              
        #  Check that data are sampled by 0.375m ?3cm OPL intervalle.       
        if self.data_check["sample"]:
            self.log.notice("Check data sampling", 3)

            diff = data["opl"][1:] - data["opl"][:-1] 
            
            if (abs(data["opl"][0] - self.opl_min)>self.opl_tol):
                self.log.error('Incorrect range: min = %f m'%(data["opl"][0]))
                raise DataError('Incorrect range: min = %f m'%(data["opl"][0]))
            if (abs(data["opl"][-1] - self.opl_max)>self.opl_tol):
                self.log.error('Incorrect range: max = %f m'%(data["opl"][-1]))
                raise DataError('Incorrect range: max = %f m'%(data["opl"][-1]))
            if ((abs(diff)-self.opl_sampling)>self.opl_tol).any():
                self.log.error('Incorrect sampling, some does not follow %.3f+-%.3f '%(self.opl_sampling, self.opl_tol))
                raise DataError('Incorrect sampling, some does not follow %.3f+-%.3f '%(self.opl_sampling, self.opl_tol))

        mask = np.zeros(data.shape, dtype=bool)             
        if self.data_check["nan"]:
            opl = data['opl']
            self.log.notice("Check NaN values", 3)            
            for key in ["doms", "incl", "yctr", "zctr", "yend", "zend"]:
                v = data[key]
                Nv = len(v)
                test = np.isnan(v)
                ou = np.where(test)[0]
                Nou = len(ou)
                if not Nou: continue

                self.log.notice("Found %d NaN values for parameter %s"%(Nou, key), 1)

                for iou,i in enumerate(ou):
                    if (iou+1)<Nou and ou[iou+1]==i+1:
                        self.log.error("At least two consecutives NaN values. Cannot fixe delirium data")
                        raise DataError("At least two consecutives NaN values. Cannot fixe delirium data")                        
                    if i==0:                        
                        v[i] = v[i+1]
                    elif i==(Nv-1):
                        v[i] = v[i-1]        
                    else:
                        v[i] = np.interp(opl[i], opl[[i-1,i+1]], v[[i-1,i+1]])                                     

                #mask += np.isnan(data[key]) 


        if self.data_check["glitches"]:
            self.log.notice("Checking glitches", 3)

            for key in ["incl", "yctr", "zctr", "yend", "zend"]:
                mask += self.filter_gliches(data[key])

        # :TODO: check what to do with invalid data 
        # in mathlab Replace invalid FOGALE data by extrapolation using second order fit of                        
        self.data = data        
        ## close the file
        if self.f:
            self.f.close()

    def filter_gliches(self, x):
        # :TODO: detect glitches
        # self.log.warning("Glitches check not yet implemented")
        return np.zeros(x.shape, dtype=bool)            

    def get_tunnel_temp(self):
        """ return the tunnel temperature has red in the delirium header 
                
        Outputs
        -------
        tunnel_temp : float
            Tunnel temperature or -999.99 if unknown
        """
        return self.header.get("tunnel_temp", -999.99)    

    def get_data(self):
        """ Get the raw data as written in file 

        """
        if self.data is None:
            self.load_data()
        return self.data

    def data2matrix(self, data, keys=None ):
        """ Return the input data in a matrix """
        if keys is None:
            keys = data.dtype.names
        return np.matrix( [data[k] for k in keys] ).T    
                
    def get_opl(self, period=False):
        """ return the data opl 

        Parametes
        ---------
        period : float, optional
            If given wrap the opl for that period
            
        Ouputs
        -----
        opl : array (vector)
            The optical Path Length
            
        """
        opl = self.get_data()["opl"]    
        if period not in [False, None]:
            period = self.parse_period(period)            
            opl = np.mod(opl, period)*2*np.pi/period
        return opl  

  
