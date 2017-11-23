""" A bunch of misc functions """
import numpy as np
import re

def get_date(filename):
        """ try to read the date from filename 

        Parameters
        ----------
        filename : string 
        
        Outputs
        -------
        date : string
            date of the form 'yyyymmdd' or '' if failed                         
        """    
        m = re.match("^.*_([0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9])[^0-9]", filename)
        if m:
            return m.groups()[0]
        return ''


def attitude2coord(V):
    """ transform a 3 vector of FOGALE position to a 3x3 coordinate matrix 

    Note: 
        y_axis = -H_axis
        z_axis =  V_axis 

    """
    return  np.array(
            [ 
              [0, V[2], -V[1]],
              [-V[2], 0, V[0]],
              [ V[1], -V[0], 0]
            ]
            )

def rad2arcsec(angle):
    """ convert an angle in rad to arcsec """
    f_rad2arcsec = 180./np.pi*3600
    return angle*f_rad2arcsec



def range2mask(x, range):
    """ from an array of value and a range make a mask of good data

    Parameters
    ----------
    x : array like
    range : None or 2 tuple 
        (xmin, xmax)
        xmin and xmax can be None -> no limits 

    Outputs
    -------
    mask : array of boolean
     (x>=xmin) & (y<=xmax)        
    """
    if range is None:
        return mask
    mask = np.ones(x.shape, dtype=bool)
    xmin, xmax = range
    if xmin is not None:
        mask *= x>=xmin
    if xmax is not None:
        mask *= x<=xmax
    return mask
            


def remove_low_order(x, y, order, mask=None, xrange=None):
    """ remove low order variation on data

    Parameters
    ----------
    x : array like
        data along the x axis
    y : array like
        data along the y axis
    order: int
        fit order

    mask : array of boolean, optional
        the mask of good pixel from where to fit. 
    xrange : 2 tuple, optional                       
        (xmin, xmax) range of the fit  

    Returns
    -------
    y_without_low_order : array like y 
        the y data without the low orders
    p : array of size order+1
        The polynome coefficient used to fit the low order      
    """
    if order is None:
        return y, []
    if mask is None:
        mask = np.ones(x.shape, dtype=bool)

    if xrange is not None:
        mask *= range2mask(x, xrange)
        
    p = np.polyfit(x[mask],
                    y[mask], order)   
    return y-get_model(x, p), p   


def get_model(x, p):
    """ from x array and a polynome return model = x*p[0]**n + x*p[1]**(n-1) + ... + p[-1]

    Parameters
    ----------
    x : array like 
    p : array like 
        the polynome coeeficient, the highest is first 

    Returns
    -------
    model : array like x
        model = x*p[0]**n + x*p[1]**(n-1) + ... + p[-1]     

    """
    x = np.asarray(x)
    if not len(p):
        return np.zeros(x.shape, dtype='f8')

    model = np.ones(x.shape, dtype='f8')*p[-1]        
    for exp,c in enumerate(p[0:-1][::-1], start=1):
        model = model +  c*x**exp
    return model    


def recarray_matrix_convertion(C, reca, keys, newkeys=None, copy=False):
    """ take a conversion matrix and a recarray and return a new recarray

    Parameters
    ----------
    C : array like
        the conversion matrix of dim MxN 
    reca : recarray like
        the input recarray dim N x Ndata. 
    keys : list
        keys list of dimension N. They are the column taken from the recarray 
        to perform the multiplication
    newkeys : list or None
        the output keys list of dimension M
        If None, the output is left has an regular array 
    copy : bool, optional
        if copy = True, the output recarray is a copy of the reca input.
        altered columns are defined by keys.
        Note that if copy is True, the matrix is expected to be rectangular

    Outputs
    -------
    new_array :  recarray, or array 
        the converted rearray of dimension Ndata       
        or array of dimension M x Ndata (if newkeys is None and copy=False)
    """
    M,N = C.shape
    if N!=len(keys):
        raise ValueError("keys list not aligned with conversion matrix %d!=%d"%(N, len(keys)))

    if copy:
        if N!=M:
            raise ValueError("with copy=True the matrix must be rectangular")    
        if newkeys and newkeys!=keys:
            raise ValueError("with copy=True newkeys must be None or identical to keys")            

    if newkeys and M != len(newkeys):
        raise ValueError("newkeys list not aligned with conversion matrix %d!=%d"%(M, len(newkeys)))
    newM = np.asmatrix(C) * np.asmatrix([reca[k] for k in keys])    
    if newkeys is None and not copy:
        return np.asarray(newM)
   
    if copy:
        newRec = reca.copy()
        for i,key in enumerate(keys):
            newRec[key] = newM[i,:]
    else:    
        dtype = C.dtype
        dtypes = [(k,dtype) for k in newkeys]
        newRec = np.recarray( reca.shape, dtype=dtypes)
        for i,key in enumerate(newkeys):
            newRec[key] = newM[i,:] 

    return newRec
