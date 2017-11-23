""" Define classes for parameter 

	At the end of the file, parameters used for delirium 
	are defined.
"""
import numpy as np


class Parameter(object):
    """ Define data parameters properties """
    def __init__(self, name,dtype="f8", unit="", 
                 period=None, order=0, description="", label=''):
        self.name  = name
        self.dtype = dtype
        self.unit  = unit
        self.description = description
        self.period = period 
        self.order  = order
        self._label = label

    def wrap(self, x, period):
        if period is None:
            return x
        return np.mod(x, period)*2*np.pi/period

    def get_dtype(self):
        return (self.name, self.dtype)    

    def rename(self, newname):
        return self.__class__(newname, dtype=self.dtype, unit=self.unit, 
                              period=self.period, order=self.order, description=self.description, label=self._label
                            )    

    @property
    def label(self):
        return self._label if self._label else self.name
    @label.setter
    def label(self, label):
        self._label = label

            

class Parameters(list):
    """ list of parameters """
    def __init__(self,*args):
        list.__init__(self, *args)
        for i,param in enumerate(self):
            if not isinstance(param, Parameter):
                self[i] = Parameter(*param)

    def get_dtypes(self):
        return [item.get_dtype() for item in self]

    def get_keys(self):
        return [item.name for item in self]

    def get(self, name):
        for item in self:
            if item.name == name:
                return item
        raise KeyError("cannot find parameter '%s'"%name)

    def restrict(self, *args):
        return Parameters([self.get(name) for name in args])


parameters = Parameters([
            #name, format, unit, typical period, fitting degree, description 

            ("time",  "f8", "hhmmss.s", None,    None, "time stamp"), #[hhmmss.s] 
            ("opl" ,  "f8", "m",        None,    None, "opl  distance"), #[m]   #1
            ("doms",  "f8", "mm",       "wheel", 2), #[mm]     
            ("incl",  "f8", "rad",      "wheel", None, "inclination"), #[rad] #2
            ("yctr",  "f8", "mm",       "wheel",  1, "horizontal position of center sensor"), #[mm]  #3
            ("zctr",  "f8", "mm",       "support",2, "vertical position of center sensor"), #[mm]  #4 
            ("yend",  "f8", "mm",       "wheel",  1, "horizontal position of end sensor"), #[mm]  #5 
            ("zend",  "f8", "mm",       "wheel",  2, "vertical position of end sensor"),   #[mm]  #6
            ("fy","f8","mm", "wheel", 1, "horizontal position of fogale sensor"), #[mm]  #5 
            ("fz","f8","mm", "wheel", 2, "vertical position of fogale sensor"),   #[mm]  #6
            ("horizontal", "f8", "mm", "wheel", 1), 
            ("vertical"  ,   "f8", "mm", "support", 2), 
            ("phi",   "f8", "rad", "wheel", 0),
            ("theta", "f8", "rad", "wheel", 1),
            ("psi",   "f8", "rad", "wheel", 0), 
            ("x","f8","m", "wheel", None, "rail opl"),
            ("y","f8","mm","wheel", 1, "rail deformation in horizontal"),
            ("z","f8","mm","wheel", 2, "rail deformation in vertical"),

            #("y","f8","mm","wheel", 1, "rail deformation in horizontal"),
            #("z","f8","mm","wheel", 2, "rail deformation in vertical"),
            ("supports","i8", "num", None, None, "Support number"),
            #("H","f8","mm", "support", None, "support horizontal deformation"), # not used !
            #("V","f8","mm", "support", None, "support vertical deformation"),   # not used
            ("Hcorrection", "f8","mm", "support", None, "support horizontal correction"),
            ("Vcorrection", "f8","mm", "support", None, "support vertical correction"),
            ## histerezis parameter 
            ("psi_direct",    "f8", "rad", None, None),
            ("psi_reverse",   "f8", "rad", None, None),
            ("theta_direct",  "f8", "rad", None, None),
            ("theta_reverse", "f8", "rad", None, None),

            ("psi_diff", "f8", "rad", None, None),            
            ("theta_diff", "f8", "rad", None, None),

        ] 
    )


