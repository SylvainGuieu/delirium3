from .log import log
from .rail import Rail 
from .delirium import Delirium, DataError, DataUtils
from .carriage import Carriage
from .sensors import Sensors
from .supports import Supports
from .parameters import parameters
from .io import _path2dl
from .computing import (remove_low_order, recarray_matrix_convertion, range2mask)
from .log import Log
import numpy as np 


try:
    unicode
except NameError:
    basestring = (str, bytes)


class DelayLineState(DataUtils):
    """ a delay line at a given moment (built from a delirium file)"""                            
    ## define a dictionary of ranges for easy calling 
    ## The lookup dictionary can be adjusted at __init__ to 
    ## match the individual delay line configuration 
    ## the range values correspond to the ctr position
    ranges_lookup = {
        "full": (None, None),        
        "large": (9.75, 119.25),
        "good": (10, 110), #used for histeresis, wobble fit
        "conservative":(8.9, 88.51) # used in matplotlib version 
                                    # to remove low order 
    }
    
    ## set the parameters
    params = parameters
    data_aliases = []
    for p in Carriage.params:
        data_aliases.append( (p.name,"carriage."+p.name))
    for p in Sensors.params:
        data_aliases.append( (p.name,"carriage.sensors."+p.name))
    for p in Rail.params:
        data_aliases.append( (p.name,"rails."+p.name) )
    for p in Supports.params:
        data_aliases.append( (p.name,"support."+p.name))
    data_aliases.append(("opl", "delirium.opl"))


    def __init__(self, dlnum=None, deliriumfile=None):
        if dlnum is None:
            dlnum = _path2dl(deliriumfile)
            if dlnum is None:
                raise ValueError("You must provide a delay line number or a valide deliriumfile")

        if dlnum<1 or dlnum>6:
            raise ValueError("Delay Line number must be between 1 and 6 got %s"%dlnum)
        self.num = dlnum

        ##
        # build the carriage attitude position from delirium file  
        self.carriage = Carriage(self)
        self.carriage.set_delirium(deliriumfile) 

        ## 
        # build the rail from carriage wheel contact points            
        self.rail = Rail(self)
        self.rail.set_carriage(self.carriage)
                
        ## make a copy to avoid to change the class dictionary but the object dictionary
        self.ranges_lookup = self.ranges_lookup.copy()

        params_lookup = { }
    

    def reload(self):
        self.carriage.reload()
        self.rail.set_carriage(self.carriage)
            
    def set_delirium(self, delirium):
        self.carriage.set_delirium(delirium) 
        self.rail.set_carriage(self.carriage)    


    def get_id(self):
        """ return a unique tuple id """
        return self.carriage.get_id()

    def period(self, period):
        """ return a numerical value for period 

        Parameters
        ----------
        period : None, float or string
            if string must be 'wheel' or 'support', the corresponding period is returned
            if None or value, period is returned as it is

        Output:
        period : float or None     
        """
        if period is None:
            return None
        if isinstance(period, basestring):    
            if period == "wheel":
                r = self.carriage.wheel_radius
                return 2*2*np.pi*r# period in opl
        
            if period == "support":
                ## Not sure where this 1.07 is comming from 
                ## keep for now to have the same behavior than matplotlib
                return Rail.support_sep*2 #/1.07 # Distance between rail supports (x 2)            

            raise KeyError("period name must be 'wheel' or 'support' got '%s'"%period)    

        return float(period)   
    parse_period = period     

    def range(self, range, end=False):
        """ return a 2 tuple numerical range 

        Parameters
        ----------
        range : None, 2 tuple like or string
            if string the corresponding numerical range is returned.
            the string must match one of the item in the dl.ranges_lookup lookup table

            if None or tuple return as it is 
        end : bool, optional
            if True the distance between sensor (in opl) is removed from ranges

        Outputs
        -------
        range : tuple
            (min, max) in opl [m]    
        """
        if range is None:
            return None
        if isinstance(range, basestring):
            try:
                vrange = self.ranges_lookup[range]
            except KeyError:
                raise ValueError("range '%s' not understood, must be one of %s"%(range, ", ".join(['%r'%k for k in self.ranges_lookup.keys()])))        
            else:
                return vrange    
        try:
            _,_ = range
        except:    
            raise ValueError("range must be a a string or a 2 tuple got %s"%range)        
        vrange = range
        if end:
            vrange = tuple(None if r is None else r-self.carriage.sensor_distance*2 for r in vrange)
        return vrange    
    parse_range = range    

    

    def get_opl(self, period=False):
        """ return the data opl 

        Parametes
        ---------
        period : float or string, optional
            If given wrap the opl for that period
            If period is a string, look for that period in the delay line definition
            example 'wheel' or 'support' define a period related to wheel wobble or rail support wobble

        Ouputs
        -----
        opl : array (vector)
            The optical Path Length
            
        """
        return self.delirium.get_opl(period=period)

    def get_date(self):
        """ get the date 'yyyymmdd' of the used delirium file"""        
        return self.delirium.get_date()

    @property
    def reverse(self):
        """ Flag True if delirium taken reversely """
        return self.delirium.reverse
            
    @property
    def delirium(self):
        if self.carriage.delirium is None:
            raise ValueError("Delirium data not loaded")
        return self.carriage.delirium

    @property
    def sensors(self):
        return self.carriage.sensors

    @property
    def supports(self):
        return self.rail.supports    
    
    def get_summary(self, nflat=5):
        _, extras = self.get("psi", filterWobble=True, removeLowOrder=True, 
                             arcsec=True, extras=True)
        cor = self.supports.get_corrections()
        Nv = len(cor[abs(cor['V'])>0])
        Nh = len(cor[abs(cor['H'])>0])
        maxH = np.max(np.abs(cor['H'])) if Nv else 0.0
        maxV = np.max(np.abs(cor['V'])) if Nh else 0.0

        opl = self.get("opl")

        fpoints, angles = self.carriage.get_flat_points(nflat)
        return {
            "date":self.delirium.get_date2(),
            "num":self.num,
            "wobble":extras.get('wobble_amplitude', 0.0),
            "Nv": Nv,
            "Nh": Nh,
            "maxH":maxH, "maxV":maxV,
            "flatpoints": ", ".join("%4.2f"%p for p in opl[fpoints])
        }


    def write_summary(self, wf=None, nflat=5):
        """ write summary in a text format
        
        returned in a text string or writen with the wf function
        e.g.: dl.write_summary( wf=open('path/to/file.html').write )
                
        """
        if wf is None:
            wf = lambda txt: None
        
        txt = """
Delirium of: {date} DL{num}
Number of corrections : 
    Vertical : {Nv}, amp. max={maxV:4.2f} [micron]
    Horizontal : {Nh}, amp. max={maxH:4.2f} [micron]
5 flattest points:
    {flatpoints} [opl]
Wobble Amplitude (psi): {wobble:4.2f} [arcsec]
""".format(**self.get_summary(nflat=nflat))
        wf(txt)
        return txt

    def write_summary_html(self, wf=None, nflat=5):
        """ write summary in a html table 
        
        returned in a text string or writen with the wf function
        e.g.: dl.write_summary_html( wf=open('path/to/file.html').write )

        """

        def warning(summary, key, test):
            val = summary[key]
            if test(val):
                summary[key+"Warning"] = "style='color:red'"
            else:
                summary[key+"Warning"] = ""

        if wf is None:
            wf = lambda txt: None

        summary = self.get_summary(nflat=nflat)
        warning(summary, "wobble", lambda v:v>=self.carriage.maxwobblewarning)
        warning(summary, "maxV", lambda v: abs(v)>=self.supports.maxcorwarning*1e3)
        warning(summary, "maxH", lambda v: abs(v)>=self.supports.maxcorwarning*1e3)

        txt = """<table class='summary'>
<tr><td colspan=2 >{date} DL{num}</td></tr>
<tr><td>Vertical Corrections</td><td {maxVWarning}>{Nv}, amp. max={maxV:4.2f} [micron]</td></tr>
<tr><td>Horizontal Corrections</td><td {maxVWarning}>{Nh}, amp. max={maxH:4.2f} [micron]</td></tr>
<tr><td>Flattest Points</td><td>
    {flatpoints} [opl]</td></tr>
<tr><td>Wobble Amplitude (psi)</td><td {wobbleWarning}>{wobble:4.2f} [arcsec]</td></tr>
</table>
""".format(**summary)
        wf(txt)
        return txt


class DelayLineStates(object):
    """ A list of DelayLineState 

    Can contain any state (any dl, at any moments)    
    """
    def __init__(self, states=[]):
        self.states = []
        for state in states:
            if isinstance(state, basestring):
                self.add_delirium(state)
            else:
                self.add_state(state)        

        self.log = Log(context=("DELIRIUM",))

    def __getitem__(self, index):
        return self.states[index]

    def add_state(self, dlstate):
        """ Add a DelayLineState object to the delayline 

        Parameters
        ----------
        dlstate : DelayLineState
        """
        if not isinstance(dlstate, DelayLineState):
            raise ValueError("expecting a DelayLineStatet object got a %s"%type(state))        
        self.states.append(dlstate)
    
    def add_delirium(self, filename): 
        """ Add a new DelayLineState from a delirium file 

        If the the is corrupted (DataError raised) the dl state 
        is not added.  

        Parmeters
        ---------
        filename : string
            delirium path


        """
        try: 
            dlstate = DelayLineState(deliriumfile=filename)
        except DataError as e:
            self.log.error("file '%s' : %s"%(filename, str(e)))
            self.log.warning("file '%s' will be ignored "%filename)        
        else:                   
            self.states.append(dlstate)
        
    @property
    def size(self):
        """ delay line size """
        return len(self.states)
        

    def get(self, key, removeLowOrder=None, filterWobble=None, 
            raw=False, arcsec=False, period=None, order=None, 
            fitrange=None,  extras=False, inside="recarray", 
            opls=None, supports=None, wobbleFitrange=None, 
            wobbleOrder=None
            ):
            """ see doc of DelayLineState.get """
            ## ask the first one to build the rest
            first = self.states[0].get(key,
                        removeLowOrder=removeLowOrder, filterWobble=filterWobble, 
                        arcsec=arcsec, period=period, order=order, fitrange=fitrange,  
                        extras=extras, inside=inside, raw=raw, opls=opls, supports=supports, 
                        wobbleFitrange=wobbleFitrange, 
                        wobbleOrder=wobbleOrder
                        )
            if extras:
                first, firstextras = first
                allextras = [firstextras]
                
            if isinstance(first, np.recarray):
                output = np.recarray( (self.size,)+first.shape, dtype=first.dtype)
                def set(item,index):
                    output[index] = item

            elif isinstance(first, (np.ndarray, np.matrix, np.float, np.int, np.bool, np.str)):
                output = np.ndarray( (self.size,)+first.shape, dtype=first.dtype)
                def set(item,index):
                    output[index] = item
            elif isinstance(first, (list,dict)):
                output = []
                def set(item,index):
                    output.append(item)
            else:
                raise ValueError("un-recognized output type %s"%type(first))        
            
            set(first,0)
            for index, state in enumerate(self.states[1:], start=1):
                item = state.get(key,
                        removeLowOrder=removeLowOrder, filterWobble=filterWobble, 
                        arcsec=arcsec, period=period, order=order, fitrange=fitrange,  
                        extras=extras, inside=inside, raw=raw, opls=opls, supports=supports, 
                        wobbleFitrange=wobbleFitrange, 
                        wobbleOrder=wobbleOrder
                        )        
                if extras:
                    item, itemextras = item
                    allextras.append(itemextras)
                set(item, index)

            if extras:
                return output, allextras
            return output        
    get.__doc__ = DelayLineState.__doc__                

class DelayLine(DelayLineStates):
    """ DelayLine object with several states for monitoring purpose 

    All DelayLineState child object must be of the same delay line number
    """
    def __init__(self, num, states=[]):
        self.num = num
        DelayLineStates.__init__(self, states)

    def add_state(self, dlstate):
        """ Add a DelayLineState object to the delayline """
        if not isinstance(dlstate, DelayLineState):
            raise ValueError("expecting a DelayLineStatet object got a %s"%type(state))
        if dlstate.num != self.num:
            raise ValueError("expecting a delay line state of delay line %d, got #%s"%(self.num, dlstate.num))            
        DelayLineStates.add_state(self, dlstate)

    def add_delirium(self, filename):
        dlnum = _path2dl(filename)   
        if dlnum != self.num:
            raise ValueError("expecting a delirium file of delay line %d, got #%s"%(self.num, dlnum))           
        DelayLineStates.add_delirium(self, filename)                
        
class DelayLineHysteresis(DataUtils, DelayLine):
    """ DelayLine object with 2 states reverse and direct 
    
    One state must have been constructed from a 'direct' delirium file, 
        the other from a a 'reverse' delirium.
    """

    opl_histeresis_fit_range = "good"    
    nbr_point = 5 # Number of Cat eye flat points to be selected (by default = 5)
    params = parameters.restrict("opl","theta_direct","psi_direct",
                                "theta_reverse","psi_reverse", 
                                "theta_diff", "psi_diff")
    def __init__(self, direct, reverse):

        if isinstance(direct, basestring):
            direct = DelayLineState(deliriumfile=direct)
        if isinstance(reverse, basestring):
            reverse = DelayLineState(deliriumfile=reverse)    

        if direct.reverse == reverse.reverse:
            raise ValueError("The two states must be one direct, one reverse")
        if direct.num != reverse.num:
            raise ValueError("direct and reverse state are not from the same delay line")                
        self.num = direct.num    
        if direct.reverse is True:
            # swap them 
            direct, reverse = reverse, direct             

        DelayLine.__init__(self, direct.num, [direct,reverse])            
        self.parse_range = self.direct.range
        self.parse_period = self.direct.period

        self.compute_histeresis()

        self.log = Log(context=("DELIRIUM", "DL%d"%self.num))
    @property
    def direct(self):
        """ The direct DelayLineState """
        return self[0]
    @property
    def reverse(self):
        """ The reverse DelayLineState """
        return self[1]
    
    def get_id(self):
        d, r = self.direct.get_id(), self.reverse.get_id()
        num = d[0]
        date = d[1]
        if r[1]!=date:
            date += "+"+r[1]
        return (num, date, False)

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
        opl = self.data["opl"]    
        if period not in [False, None]:
            period = self.direct.period(period)            
            opl = np.mod(opl, period)*2*np.pi/period
        return opl  

    def get_diff(self, key, removeLowOrder=None, filterWobble=None, 
            arcsec=False, period=None, order=None, fitrange=None,              
            indexes=None, opls=None, wobbleFitrange=None, wobbleOrder=None):

            keys = ["direct.opl", "direct."+key, "reverse.opl", "reverse."+key]    
            opl_direct, direct, opl_reverse, reverse = self.get(keys, 
                                    inside='list',
                                    removeLowOrder=removeLowOrder, 
                                    filterWobble=filterWobble, 
                                    arcsec=arcsec, period=period, order=order, 
                                    fitrange=fitrange, indexes=indexes, opls=opls, 
                                    wobbleFitrange=wobbleFitrange, wobbleOrder=wobbleOrder
                                    )
            # re interp reverse on direct
            reverse = np.interp(opl_direct, opl_reverse, reverse)  
            return opl_direct, direct-reverse            

    def compute_histeresis(self):

        opl_direct = self.direct.get_opl()
        opl_reverse = self.reverse.get_opl()

        rge = self.parse_range(self.opl_histeresis_fit_range)
        ## get the index points of the range
        index_direct  = self.direct.carriage.iopl(rge)
        index_reverse = self.reverse.carriage.iopl(rge)

        rmlo, fw = False, False
        theta_direct = self.direct.carriage.get("theta",removeLowOrder=rmlo, filterWobble=fw)
        psi_direct   = self.direct.carriage.get("psi",removeLowOrder=rmlo, filterWobble=fw)

        theta_reverse = self.reverse.carriage.get("theta",removeLowOrder=rmlo, filterWobble=fw)
        psi_reverse   = self.reverse.carriage.get("psi",removeLowOrder=rmlo, filterWobble=fw)

        residual_theta_direct, pol_theta_direct,_ = self.direct.carriage.remove_low_order_of_param("theta", opl_direct, theta_direct, order=1, fitrange=rge)
        residual_psi_direct, pol_psi_direct,_ = self.direct.carriage.remove_low_order_of_param("psi", opl_direct, psi_direct, order=1, fitrange=rge)

        residual_theta_reverse, pol_theta_reverse,_ = self.reverse.carriage.remove_low_order_of_param("theta", opl_reverse, theta_reverse, order=1, fitrange=rge)
        residual_psi_reverse, pol_psi_reverse,_ = self.reverse.carriage.remove_low_order_of_param("psi", opl_reverse, psi_reverse, order=1, fitrange=rge)

        mean_theta_direct = np.mean(theta_direct[index_direct])
        mean_psi_direct   = np.mean(psi_direct[index_direct])

        mean_theta_reverse = np.mean(theta_reverse[index_reverse])
        mean_psi_reverse   = np.mean(psi_reverse[index_reverse])

        ##
        # compute the difference with 2 methods 
        # one is just the mean difference _m
        # the other the fit constant term difference _c
        theta_diff_c = pol_theta_direct[1]-pol_theta_reverse[1]
        theta_diff_m = mean_theta_direct - mean_theta_reverse
        
        psi_diff_c = pol_psi_direct[1]-pol_psi_reverse[1]
        psi_diff_m = mean_psi_direct - mean_psi_reverse

        theta_mean = 0.5*(pol_theta_direct[1]+pol_theta_reverse[1])
        psi_mean   = 0.5*(pol_psi_direct[1]+pol_psi_reverse[1])

        ## remove the low order as it was written on Matlab
        ## remove just the mean on thepsi and linear fit of the direct on theta
        theta_direct  = theta_direct - (theta_mean + pol_theta_direct[0]*opl_direct)
        theta_reverse = theta_reverse - (theta_mean + pol_theta_direct[0]*opl_reverse)

        psi_direct = psi_direct - psi_mean
        psi_reverse = psi_reverse - psi_mean

        angle = np.sqrt( psi_direct*psi_direct + theta_direct*theta_direct)

        flat_point_indexes = angle.argsort()[0:self.nbr_point]
        # resort the indexes so they are in opl order        
        flat_point_indexes.sort()

        ## interpolate the reverse value to the direct 
        theta_reverse = np.interp(opl_direct, opl_reverse, theta_reverse)        
        psi_reverse = np.interp(opl_direct, opl_reverse, psi_reverse)

        theta_diff = theta_direct - theta_reverse
        psi_diff = psi_direct - psi_reverse               
        self.data = np.rec.fromarrays([opl_direct,theta_direct,psi_direct,theta_reverse,psi_reverse, theta_diff, psi_diff], names=[p.name for p in self.params])

        self.theta_diff_c = theta_diff_c
        self.psi_diff_c = psi_diff_c

        self.theta_diff_m = theta_diff_m
        self.psi_diff_m = psi_diff_m

        self.theta_diff_mean = np.mean(theta_diff[index_direct])
        self.theta_diff_std  = np.std(theta_diff[index_direct])        

        self.psi_diff_mean = np.mean(psi_diff[index_direct])
        self.psi_diff_std  = np.std(psi_diff[index_direct])        

        self.flat_point_indexes = flat_point_indexes





