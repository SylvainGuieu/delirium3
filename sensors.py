import numpy as np
from .computing import recarray_matrix_convertion, range2mask
from .io import _path2dl
from .parameters import parameters, Parameters
from .delirium import Delirium, DataUtils
from .log import Log

""" 
WPS roll angle correction of the Fogale detector 
angles are in radian 
"""
fogale_angle_lookup = {

    1: dict(
            ctr = -1.82e-3, 
            end = -9.75e-3
            ),
    2: dict(
            ctr = -13.43e-3,
            end = -9.03e-3
           ),
    3: dict(
            ctr = -10.04e-3,
            end = -6.19e-3
           ),
    4: dict(
            ctr = 1.73e-2,
            end = 1.6e-3
         ),
    5: dict(
            ctr = 7.5e-3,
            end = 8.1e-3
        ),
    6: dict(
            ctr = 4.3e-3,
            end = -0.9e-3
     )
}

class _Sensor_(object):
    def get_data(self):
        return self.data

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


class Inclinometer(_Sensor_,DataUtils):
    params = parameters.restrict("opl", "incl")
    data = None

    def __init__(self, data=None, sensors=None):
        if sensors:
            self.num = sensors.num
            self.parse_range = sensors.parse_range
            self.parse_period = sensors.parse_period
            self.get_id = sensors.get_id

        if data is not None:
            self.set_data(data)
        else:
            ## set empty dummy data
            self.set_data( np.array( [(0.0,0.0)], dtype=[("opl","f"),("incl","f")] ).view(np.recarray))

    def set_data(self, data):        
        self.data = data
    
class Fogale(_Sensor_,DataUtils):
    """ Represent one Fogale detector object """
    #params = parameters.restrict("opl","fy","fz")
    params = Parameters([parameters.get("opl"), 
                         parameters.get("fy").rename("y"),
                         parameters.get("fz").rename("z")
                        ]
                    )
    data = None
    def __init__(self, phi=0.0, data=None, name="", sensors=None):
        if sensors:
            self.num = sensors.num
            self.parse_range = sensors.parse_range
            self.parse_period = sensors.parse_period
            self.get_id = sensors.get_id

        self.phi = phi
        self.D2Dc = np.matrix([
                      [np.cos(self.phi),-np.sin(self.phi)],
                      [np.sin(self.phi), np.cos(self.phi)]
                     ])
        self.D2Dc_keys = ( ["y","z"], None) 
        ## log sent to stdout by default 
        self.log = Log(context=("DELIRIUM", "FOGALE %s"%name))
        self.name = name
        if data is not None:
            self.set_data(data)
        else:
            ## set empty dummy data
            self.set_data( np.array( [(0.0,0.0,0.0)], dtype=[("opl","f"),("y","f"), ("z","f")] ).view(np.recarray) )

    def d2dc(self, D):
        """ measured Data to WPS Roll correction error Data 

        Parameters
        ----------
        D : recarray 
            the measured raw data with keys 'yctr','zctr','yend','zend'

        Outputs:
        Dc : recarray
            the corrected data, with the same column name/order than input data   
        """
        return recarray_matrix_convertion(self.D2Dc, D, *self.D2Dc_keys, copy=True)

    def set_data(self, D):        
        self._raw_data = D                    
        self.log.notice("Correct raw data from WPS roll Error")
        Dcorrected = self.d2dc(D)         
        self.data = Dcorrected        



class Sensors(DataUtils):
    ## parameters used by Sensors object
    params = parameters.restrict("time","opl","doms","incl",
                                 "yctr","zctr","yend", "zend"                                
                                )
    ## repartition of data, they are taken from sensors sub-objects 
    # data_aliases = [
    #     ("opl" ,("ctr","opl")),
    #     ("yctr",("ctr", "y")), 
    #     ("zctr",("ctr", "z")),
    #     ("yend",("end", "y")), 
    #     ("zend",("end", "z")), 
    #     ("incl",("incl", "incl"))
    # ]
    data_aliases = [
        ("opl", "ctr.opl"),
        ("yctr","ctr.y"), 
        ("zctr","ctr.z"),
        ("yend","end.y"), 
        ("zend","end.z"), 
        ("incl","incl.incl")
    ]

    ## data is filled by set_delirium
    data = None

    def __init__(self, carriage, delirium=None):
        """ Sensors object 
        
        Represent the 2 Fogale sensors on carriage.

        The data are taken from a delirium object (or file) and are corrected from
        the WPS roll error

        Available data parameters are "opl","incl","yctr","zctr","yend","zend"


        """
        self.num = carriage.num
        self.parse_range = carriage.parse_range
        self.parse_period = carriage.parse_period
        ## log sent to stdout by default 
        self.log = Log(context=("DELIRIUM", "DL%d"%self.num))

        if delirium is not None:
            self.set_delirium(delirium)
        else:
            self.ctr = Fogale(name="ctr",sensors=self)
            self.end = Fogale(name="end",sensors=self)  
            self.incl = Inclinometer(sensors=self) 

    def reload(self):
        self.delirium.reload()
        self.set_delirium(self.delirium)             

    def set_delirium(self, delirium):
        """ Set the delirium file to this sensors and perform Roll error correction 

        Parameters
        ----------
        delirium : string
            file path pointing to a delirium file

        """
        if isinstance(delirium, Delirium):
            self.delirium = delirium                        
            if delirium.dlnum != self.num:
                raise ValueError("Setting a delirium #%d on sensor #%d "%(delirium.dlnum, self.num))

        else:
            if delirium:        
                dlnum = _path2dl(delirium)
                if dlnum != self.num:
                    raise ValueError("Setting a delirium #%d on sensor #%d "%(dlnum, self.num))
            self.delirium = Delirium(self, delirium)

            
        D = self.delirium.get_data()
        if D is None:
            self.log.warning("Cannot get delirium data")
            return 
        self.get_id = self.delirium.get_id
        
        angles = fogale_angle_lookup[self.num]
        self.ctr = Fogale(angles["ctr"], np.rec.fromarrays([D["opl"], D["yctr"], D["zctr"]], names="opl,y,z"), name="ctr", sensors=self)
        self.end = Fogale(angles["end"], np.rec.fromarrays([D["opl"], D["yend"], D["zend"]], names="opl,y,z"), name="end", sensors=self)    
        self.incl = Inclinometer(np.rec.fromarrays([D["opl"], D["incl"]], names="opl,incl"), sensors=self)
        #self.data = D
        return 0    
    
    def get_opl(self, period=False):
        return self.delirium.get_opl(period=period)    
    get_opl.__doc__ = Delirium.get_opl.__doc__    
    
    



class Sensors_(DataUtils):
    ## parameters used by Sensors object
    params = parameters.restrict("time","opl","doms","incl",
                                 "yctr","zctr","yend", "zend"                                
                                )
    ## data is filled by set_delirium
    data = None

    def __init__(self, carriage, delirium=None):
        """ Sensors object 
        
        Represent the 2 Fogale sensors on carriage.

        The data are taken from a delirium object (or file) and are corrected from
        the WPS roll error

        Available data parameters are "opl","incl","yctr","zctr","yend","zend"


        """
        self.num = carriage.num
        self.parse_range = carriage.parse_range
        self.parse_period = carriage.parse_period


        params_lookup = {
            1: dict(
                    phi_ctr = -1.82e-3, 
                    phi_end = -9.75e-3
                    ),
            2: dict(
                    phi_ctr = -13.43e-3,
                    phi_end = -9.03e-3
                   ),
            3: dict(
                    phi_ctr = -10.04e-3,
                    phi_end = -6.19e-3
                   ),
            4: dict(
                    phi_ctr = 1.73e-2,
                    phi_end = 1.6e-3
                 ),
            5: dict(
                    phi_ctr = 7.5e-3,
                    phi_end = 8.1e-3
                ),
            6: dict(
                    phi_ctr = 4.3e-3,
                    phi_end = -0.9e-3
             )
        }
        for k,v in params_lookup[self.num].items():
            setattr(self, k, v)

        ## build the Correction WPS Roll Error Matrix
        ## D2Dc stand for data to data corected
        self.D2Dc  = np.matrix([
                      [np.cos(self.phi_ctr),-np.sin(self.phi_ctr),0,0],
                      [np.sin(self.phi_ctr), np.cos(self.phi_ctr),0,0],
                      [0,0, np.cos(self.phi_end),-np.sin(self.phi_end)],
                      [0,0, np.sin(self.phi_end), np.cos(self.phi_end)]
                     ]) 
        self.D2Dc_keys = ( ["yctr","zctr","yend","zend"], None)               

        ## log sent to stdout by default 
        self.log = Log(context=("DELIRIUM", "DL%d"%self.num))


        ## Live this at the end. Set the delirium on this Sensor
        if delirium:
            self.set_delirium(delirium)

    def reload(self):
        self.delirium.reload()
        self.set_delirium(self.delirium)        
            
    @property
    def raw_data(self):
        """ raw data is the delirium data """
        return self.delirium.data            

    def set_delirium(self, delirium):
        """ Set the delirium file to this sensors and perform Roll error correction 

        Parameters
        ----------
        delirium : string
            file path pointing to a delirium file

        """
        if isinstance(delirium, Delirium):
            self.delirium = delirium                        
            if delirium.dlnum != self.num:
                raise ValueError("Setting a delirium #%d on sensor #%d "%(delirium.dlnum, self.num))

        else:
            if delirium:        
                dlnum = _path2dl(delirium)
                if dlnum != self.num:
                    raise ValueError("Setting a delirium #%d on sensor #%d "%(dlnum, self.num))
            self.delirium = Delirium(self, delirium)

            
        D = self.delirium.get_data()
        if D is None:
            self.log.warning("Cannot get delirium data")
            return 
        self.get_id = self.delirium.get_id
        
        self.delirium_data = D
                    
        self.log.notice("Correct raw data from WPS roll Error")
        Dcorrected = self.d2dc(D)         
        self.data = Dcorrected
        return 0

    def get_data(self):
        """ Return de sensors recarray data

        values are corrected from WPS roll error         
        """
        return self.data    

    def get_opl(self, period=False):
        return self.delirium.get_opl(period=period)    
    get_opl.__doc__ = Delirium.get_opl.__doc__

    def d2dc(self, D):
        """ measured Data to WPS Roll correction error Data 

        Parameters
        ----------
        D : recarray 
            the measured raw data with keys 'yctr','zctr','yend','zend'

        Outputs:
        Dc : recarray
            the corrected data, with the same column name/order than input data   
        """
        return recarray_matrix_convertion(self.D2Dc, D, *self.D2Dc_keys, copy=True)



