from .computing import (remove_low_order, recarray_matrix_convertion, 
                        range2mask, attitude2coord, get_model)
from .sensors  import Sensors
from .parameters import parameters
from .log import Log
from .delirium import DataUtils
import numpy as np

class Carriage(DataUtils):
    """ Carriage contain the physical carriage parameters
    and coordinate matrix conversions.
    Its data is the attitude computed from the ctr and end sensor measurement (corrected for WPS roll error)

    available data parameters are "opl", "horizontal", "vertical", "phi", "theta", "psi"

         
Side View
---------

          196                 Vaxis = Z / \   |--------                             \
         >.....<                         .    |        |                             \
          |   |                          |    |        |                       theta  \
             +++                         .    |       /        +++             --------
           +     +-----------------------|-------------------+     +
          +       +                      .                  +       +
    \ /    +     +-----------------------|-------------------+     +
 48  |       +++      Uc                 .             Ue      +++                          Uc = ctr Fogale detector
    / \ ==========================================================================--> X     Ue = end Fogale detector 
                                         | 
              .       .                  .             .       .  
              |       |                  |             |       | 
              .       .     900          .       800   .       .
              |       |                  |             |       |
              .                          .                     .  
              |     1125                 |           1125      |  
                                         .

   <.... tunnel center                                               tunnel  end ....>  


Top View 
--------                                 |
                                         .
      \ /                                |
       |  ==########=====================.===================##########=================== Support Rail =====
       |                                 |
       |                                 .
 520   |                                 |                                        \/                                          /
       |   ...........Uc.................................Ue...................................Wire............               /                             
       |                                 .                                        92                                        /   psi 
       |  ==########=====================|===================#########===================== Guiding Rail =====--> X         --------
      / \                                .                                        /\             
                                         |
                                         .
                                         |
                                        \ / Haxis = -Y 

Front View
----------
                                /|\ Vaxis = Z            \
                                 .                        \
                  |--------|     |                    phi  \  
                  | o    o |                          ------
            |_____| o    o |_____|                  
       \/   |--------------------|      
       48   |               U    |
    ........Q....................Q.....>  Haxis = -Y                                       
       /\    
                            >  92 <
           >        520           <
    

    opl[m] = 1.5 #support - 5.69                                                 
    Roll angle :  phi   when lift end support wheel around guiding wheel
    Pitch angle : theta when lift center wheel
    Yaw angle   : psi   when end guiding wheel going near support rail 


    


    """
        
    wheel_distance  = 2.25 # physical separation between wheels [m]
    sensor_distance = 1.7  # physical separation between sensors [m]
    wheel_radius    = 0.196 # [m] exact design value = 195.81 mm
    rail_distance   = 0.520  # separation between rails [m]
    

    # FOGALE (sensors) positions [mm] compare to carriage reference
    #  !! in the x,y,z coordinates  (not the rail,Haxis,Vaxis) !!
    #  [sensor position along rail axis, 
    #   sensor position along Horizontal axis,
    #   sensor position along Vertical axis
    #  ]
    Fend = [ 800,92,48]
    Fctr = [-900,92,48]
    
    # Positions of wheel contact points
    Wend = [ 1000*wheel_distance/2.,0,0]
    Wctr = [-1000*wheel_distance/2.,0,0]
      
    # 
    A2Fend = attitude2coord(Fend);
    A2Fctr = attitude2coord(Fctr);

    ##################################################
    # Conversion matrix from Attitude to Measured data
    # (the inverse of this matrix is actually used)
    #   D = A2D x A 
    #
    # | opl  |    | 1  0  0  0   0  0  |    | opl |
    # | incl |    | 0  0  0  1   0  0  |    | H   |
    # | ystr | =  | 0  1  0 -48  0 -900| X  | V   |
    # | zctr |    | 0  0  1  92  900 0 |    | phi | 
    # | yend |    | 0  1  0 -48  0  800|    |theta|  
    # | zend |    | 0  0  1  92 -800  0|    | psi |
    #
    ###############################################
    A2D = np.array(
          [[1,0,0,0,0,0],
           [0,0,0,1,0,0],
           [0,1,0]+list(A2Fctr[1]),
           [0,0,1]+list(A2Fctr[2]),
           [0,1,0]+list(A2Fend[1]),
           [0,0,1]+list(A2Fend[2]),
           ]
        )
    # row / columns data keys for A2D interaction matrix
    A2D_keys = (
                ["opl","horizontal","vertical","phi","theta","psi"],
                ["opl","incl","yctr","zctr","yend","zend"]
               )
    #################################################
    #
    # D2A is the inverse matrix used to convert the measured delirium
    # to attitude coordinates 
    # 
    D2A = np.linalg.inv(A2D)
    D2A_keys = A2D_keys[::-1]

        

    ####################################################
    #
    # A2Wctr transform attitude to ctr wheel contact point
    # coordinate.
    # NB : a minus sign is apply here to y as it is on the
    #      oposite direction than H. Hold matlab program did it after
    #
    #                                         | opl |
    # | x |    | 1  0  0  0     0     0 |     | H   |                   
    # | y |  = | 0 -1  0  0     0  1125 | X   | V   | 
    # | z |    | 0  0  1  0  1125     0 |     | phi |
    #                                         |theta|
    #                                         | psi |      

    ## below, set None to leave the converted array as a normal array (not a recarray)
    A2Wctr_keys = (["opl","horizontal","vertical","phi","theta","psi"],None)         
    A2Wctr = np.matrix(np.array([np.eye(3),attitude2coord(Wctr).T]).reshape( (6,3) )).T    
    A2Wctr[1,1] = -1   # y = -Haxis
    A2Wctr[1,5] = -A2Wctr[1,5] # y = -Haxis
    ####################################################
    #
    # A2Wend transform attitude to end wheel contact point
    # coordinate.    
    # NB : a minus sign is apply here to y as it is on the
    #      oposite direction than H. Hold matlab program did it after
    #
    #                                         | opl |
    # | x |    | 1  0  0  0     0     0 |     | H   |                   
    # | y |  = | 0 -1  0  0     0 -1125 | X   | V   | 
    # | z |    | 0  0  1  0 -1125     0 |     | phi |
    #                                         |theta|
    #                                         | psi |
    #             
    A2Wend = np.matrix(np.array([np.eye(3),attitude2coord(Wend).T]).reshape( (6,3) )).T
    A2Wend[1,1] = -1 # y = -Haxis
    A2Wend[1,5] = -A2Wend[1,5] # y = -Haxis

    ## set None to leave the converted array as a normal array (not a recarray) 
    A2Wend_keys = (["opl","horizontal","vertical","phi","theta","psi"],None)            
    
    ## parameters of the attitude recarray 
    params = parameters.restrict("opl", "horizontal", "vertical", 
                             "phi", "theta", "psi")

    data = None
    raw_data = None

    sensors_data = None 
    
    maxwobblewarning = 3.0 # [arcsec] warning print if wobble above this value
    def __init__(self, dl, delirium=None):
        self.num = dl.num
        self.parse_range = dl.range
        self.parse_period = dl.period 
        #########
        ### alter here any parameters that are different from 
        ### one delay line to an other 
        #
        #        
        ## Fogal sensor object belong to carriage 
        self.sensors = Sensors(self, delirium=delirium)

        ## log sent to stdout by default 
        self.log = Log(context=("DELIRIUM", "DL%d"%dl.num))
        if delirium:
          self.set_delirium(delirium)

    def reload(self):
        self.sensors.reload()
        self.set_delirium(self.sensors.delirium)

    def set_delirium(self, delirium):
        """ set the delirium file to this carriage, the attitude is computed

        Parameters
        ----------
        delirium : Delirium or string path
        """
        self.sensors.set_delirium(delirium)
        #if self.sensors.data is not None:
        self.compute_attitude()

        self.get_id = self.delirium.get_id

    def d2a(self, D):
        """ measured data to attitude conversion 

        note, measure data must be corrected from FOGAL rotation

        The convertion matrix is taken from carriage D2A matrix 
        Parameters
        ----------
        D : recarray 
            the measured data with keys 'opl','inc',yctr','zctr','yend','zend'

        Outputs
        -------
        A : recarray 
            attitude data : 'opl', 'horizontal', 'vertical', 'phi', 'theta', 'psi'             
        """
        return recarray_matrix_convertion(self.D2A, D, *self.D2A_keys)
    def a2d(self, A):
        """ attitude to measure data conversion 

        The convertion matrix is taken from carriage A2D matrix 
        Parameters
        ----------
        A : recarray 
            attitude data : 'opl', 'horizontal', 'vertical', 'phi', 'theta', 'psi'             
            

        Outputs
        -------
        D : recarray 
            the data as measured with keys 'opl','inc',yctr','zctr','yend','zend'
        """
        return recarray_matrix_convertion(self.A2D, A, *self.A2D_keys)


    def  a2w(self, A, pos):
        """ attitude to guiding wheel contact point conversion 

        The convertion matrix is taken from carriage A2Wctr or A2Wend matrix 

        Parameters
        ----------
        A : recarray 
            attitude data of dim N : 'opl', 'horizontal', 'vertical', 'phi', 'theta', 'psi'             
        pos : string
            must be 'ctr' or 'end'. wheel position of centered or end sensor  

        Outputs
        -------
        W : array 
            the guiding wheel contact points of dim 3 x N           

        """
        if pos not in ['ctr','end']:
            raise ValueError("pos must be 'ctr' or 'end' got '%s'"%pos)
        C = self.A2Wctr if pos=='ctr' else self.A2Wend
        keys = self.A2Wctr_keys if pos=='ctr' else self.A2Wend_keys              
        return recarray_matrix_convertion(C, A, *keys)
          
    def compute_attitude(self, removeLowOrder=True, filterWobble=False):
        """ compute the attitude from the delirium data
        
        This will convert the FOGALE detector measurement into : opl, horizontal, vertical and
        the 3 angles of the carriage.
        The result is stored in a recarray inside the .data attribute.

        """
        if self.sensors.delirium is None:
          raise ValueError("No delirium file attached")  
        self.log.notice("Building the carriage attitude coordinate and orientation", 3)    
        #D = self.sensors.get_data()
        #D = self.sensors.get([p.name for p in self.sensors.params], removeLowOrder=True, filterWobble=True)

        # sensors has the data corrected from fogal tilt angle
        # at this point filterWobble should be False, the wobble will be filtered when 
        # recontructing the rail
        
        ## Old way, the low order are remove independantly to ctr and end sensors 
        ## That cause problems for the histeresis, because histeresis information 
        ## are contained inside the loworder differences between ctr and end  for theta and psi        
        # D = self.sensors.get( self.D2A_keys[0], removeLowOrder=removeLowOrder, filterWobble=filterWobble)                
        
        ## New way calculate a the string model on the ctr sensor and apply that model to 
        ## the end (sensors) with the separation offset
        opl, incl = self.sensors.get(["opl", "incl"], removeLowOrder=removeLowOrder, filterWobble=filterWobble, inside='list')

        data = {}
        for k1,k2 in [('yctr','yend'), ('zctr','zend')]:
            _, e = self.sensors.get(k1, removeLowOrder=True, extras=True, filterWobble=filterWobble)
            p1 = e['loworder_polynome']
            order = len(p1)-1
            ## distance between sensors in opl [m]
            delta = (self.Fctr[0]-self.Fend[0])*2*1e-3

            if order==2:
                a,b,c = p1
                p2 = [a, b-2*a*delta, c+a*delta**2-b*delta]
            elif order==1:
                a,b = p1
                p2 = [a, b-delta]
            elif order==0:
                p2 = p1
            else:
                raise RuntimeError("Low Order must be <=2 got %d for '%s'"%(order,k1))           

            m1 = get_model(opl, p1)
            m2 = get_model(opl, p2)
            data[k1] =  self.sensors.get(k1, filterWobble=filterWobble)-m1
            data[k2] = self.sensors.get(k2, filterWobble=filterWobble)-m2
        D = np.rec.fromarrays( [opl, incl, data['yctr'], data['zctr'], data['yend'], data['zend']], names="opl,incl,yctr,zctr,yend,zend")
                        
        ##
        # save the sensors data 
        # for plot purpose        
        self.sensors_data = D

        A = self.d2a(D)        
        self.data = A

    def get_data(self):
        """ return the attitude recarray """
        if not self.data:
          self.compute_attitude()
        return self.data

    def get_opl(self, period=False):
        """ return the opl 

        eventualy wrap the opl table with the optional given period.
        
        Parametes
        ---------
        period : float, string, optional
            If given wrap the opl for that period
            If string must be a defined period like 'wheel' or 'support'     

        Ouputs
        -----
        opl : array (vector)
            The optical Path Length    
        """
        return self.sensors.get_opl(period=period)
    get_opl.__doc__ = Sensors.get_opl.__doc__

    def get_wheels_contact_points(self, filterWobble=True, removeLowOrder=False):
        """ return the wheel contact points positions computed from attitude

        Parameters
        ----------
        filterWobble : boolean, optional 
            remove the wobble if True (default) before computing
            the rail contact point
        removeLowOrder : boolean , optional
            remove the low order if True before computing
            the rail contact point

        Ouputs
        ------
        Wctr : 3xN array
            center wheel x,y,z positions 
        Wend : 3xN array
            end wheel x,y,z positions 
        """        
        A = self.get([p.name for p in self.params], filterWobble=filterWobble,removeLowOrder=removeLowOrder)                                                      
        return self.a2w(A, "ctr"), self.a2w(A, "end")

    def get_flat_points(self, N=-1, key=None, arcsec=True):
        rmlo, fw = True, False
        if not key:
            theta = self.get("theta", order=0, removeLowOrder=rmlo, filterWobble=fw, arcsec=arcsec)
            psi   = self.get("psi",   order=0, removeLowOrder=rmlo, filterWobble=fw, arcsec=arcsec)
            angle = np.sqrt( psi*psi + theta*theta)        
        else:
            angle = np.abs(self.get(key, order=0, removeLowOrder=rmlo, filterWobble=fw, arcsec=arcsec))
        order = angle.argsort()[:N]            
        return order, angle[order]




                


    @property
    def delirium(self):
        """associated delirium file """
        return self.sensors.delirium
    





