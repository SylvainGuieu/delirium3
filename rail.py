
from .supports import Supports
from delirium3 import DataUtils
from .parameters import parameters
import numpy as np

class Rail(DataUtils):   
    ## 
    #
    # support number = support_sep*2*opl + support_offset*2
    # support number start from 1 !!!
    support_sep    =  0.75    # in physical [m]
    support_offset =  5.69/2. # in physical [m]

    ## filter carriage wobble before
    ## reconstructing the rail, normaly True
    filterCarriageWobble = True

    ## remove carriage low order before 
    ## reconstructing the rail 
    ## Should be True to remove the carriage tilts 
    removeCarriageLowOrder = True         

    data = None

    ## parameters of the rail recarray     
    params = parameters.restrict("x", "y", "z")

    def __init__(self, dl, carriage=None):
        self.num = dl.num     

        self.parse_period = dl.period
        self.parse_range  = dl.range

        self.reconstructed = False
        self.carriage = None
        ##
        ## reconstruct the rail with carriage attitude 
        ## informations (wheel contact points) 
        self.supports = Supports(self)


        if carriage:
            self.set_carriage(carriage)            

    def reload(self):
        self.carriage.reload()
        self.set_carriage(self.carriage)

    def set_carriage(self, carriage):
        if carriage.num != self.num:
            raise ValueError("Setting a carriage %d on rail %d "%(carriage.num, self.num))

        self.reconstruct(carriage)
        self.carriage = carriage        
        self.supports.compute_deformations(self)
        self.get_id = carriage.get_id
        self.supports.get_id = self.get_id  

    @property
    def raw_data(self):
        """ raw_data is data """
        return self.data
        

    def support2opl(self, support):
        """ convert support numbers to opl distance

        Parameters
        ----------
        support : array like
            upport number (starting from 1)

        Outputs
        -------
        opl : array like
            the opl distance [m] 

        Notes
        -----
        There is no check of physical values boundaries           
        """
        support = np.asarray(support)
        opl = support*(self.support_sep*2)-(self.support_offset*2)
        return  opl

    def opl2support(self, opl):
        """ convert opl to the closest support number


        Parameters
        ----------
        opl : array like
            the opl distance [m]

        Output
        ------
        support : array like
            support number (starting from 1) as a array of int            

        Notes
        -----
        Support number start from 1
        There is no check of physical values boundaries
        """
        opl = np.asarray(opl)
        support = np.round ((opl+(self.support_offset*2))/(self.support_sep*2))
        return support.astype(int)
    

    def reconstruct(self, carriage, filterCarriageWobble=None, removeCarriageLowOrder=None):
        """ reconstruct the rail from attitude carriage 

        Parameters
        ----------
        carriage : Carriage object
            the carriage object must have a delirium file attached and 
            a attitude computed wich is normaly always the case

        filterCarriageWobble :  boolean or None, optional
            filter attitude wobble before recontructing the rail
            if None the self.filterCarriageWobble attribute is taken

        removeCarriageLowOrder :  boolean or None, optional
            filter low orders before recontructing the rail
            if None the self.removeCarriageLowOrder attribute is taken

        """
        filterCarriageWobble = self.filterCarriageWobble if filterCarriageWobble is None else filterCarriageWobble
        removeCarriageLowOrder = self.removeCarriageLowOrder if removeCarriageLowOrder is None else removeCarriageLowOrder
        Wctr, Wend = carriage.get_wheels_contact_points(filterWobble=filterCarriageWobble, removeLowOrder=removeCarriageLowOrder)
        
        ## index separation between measurement with ctr and  end sensor
        ## should be 12
        index_shift = np.int(np.round(carriage.wheel_distance*2/carriage.sensors.delirium.opl_sampling))
        
        rail = (Wctr[:,index_shift:]+Wend[:,0:-index_shift])/2;

        rail = np.rec.fromarrays(rail, names='x,y,z')
        #rail['y'],_,_ = self.remove_low_order_of_param("y", rail['x'], rail['y'], fitrange="conservative")
        #rail['z'],_,_ = self.remove_low_order_of_param("z", rail['x'], rail['z'], fitrange="conservative")
        self.data = rail     
        self.reconstructed = True                                                   

    def get_data(self):  
        """ return the reconstructed rail recarray """      
        if not self.reconstructed:
            raise ValueError("rail not reconstructed")
            return self.reconstruct()
        return self.data        

    def get_opl(self, period=False):        
        """ return the rail opl 

        Parametes
        ---------
        period : float, optional
            If given wrap the opl for that period
            
        Ouputs
        -----
        opl : array (vector)
            The optical Path Length
            
        """
        opl = self.x    
        if period not in [False, None]:
            period = self.parse_period(period)            
            opl = np.mod(opl, period)*2*np.pi/period

        return opl  

    @property
    def x(self):
        if not self.reconstructed:
            raise ValueError("rail not reconstructed")
        return self.data['x']

    @property
    def y(self):
        if not self.reconstructed:
            raise ValueError("rail not reconstructed")
        return self.data['y']
        
    @property
    def z(self):
        if not self.reconstructed:
            raise ValueError("rail not reconstructed")
        return self.data['z']
