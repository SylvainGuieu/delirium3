""" define functions to fit the wobble """
import numpy as np
from .log import log
from .computing import get_model
from scipy.optimize import leastsq



class WobbleFit(object):
    results = None


    

    def __init__(self, **params):
        for k,v in params.items():
            if not k in self.__class__.__dict__:
                raise KeyError("unknown paramter '%s'"%k)
            setattr(self, k, v)             
    
    def new(self, x, y, **params):
        dct = self.__dict__
        new = self.__class__(**dict(dct,**params))
        new.setxy(x, y)
        return new            

    def setxy(self, x, y):
        self.x, self.y = np.asarray(x), np.asarray(y)                

    def fit(self):
        raise NotImplemented("")

    def model(self, x):
        raise NotImplemented("")

    def residuals(self):
        """ Return data residuals of the fit """
        return self.y-self.model(self.x)    

    def plot(self,x, mk="b-", axis=None, figure=None, **kwargs):
        """ plot the fitted model according to a x array 

        Parameters:
        -----------
        x : array like 
        mk : string, optional 
            marker short string description, default is 'b-'
        axis : matplotlib axis, optional
            if not given take the subplot 1,1 of figure
        figure : matplotlib figure, optional
            if not given take plt.figure()
        **kwargs : optional
            other plot parameter         
        """
        from matplotlib.pylab  import plt
        if axis is None:
            if figure is None:
                figure = plt.figure()
            
            axis = figure.add_subplot(111)
        return axis.plot(x, self.model(x), mk, **kwargs)


class WobbleFitPol4(WobbleFit):
    """ Fit the wobble with a polynome of degree 4 

    the x array must be wrapped from 0 to 2*pi. 
    xstart and xend parameterd gives the start end end value of x, default is 0 and 2*pi

    """
    xstart = 0
    xend = 2*np.pi 

    def fit(self):
        x, y0, xstart, xend = (getattr(self,k) for k in ["x","y","xstart","xend"])

        weight = 2*len(x)

        matx = np.matrix([np.ones(x.shape, x.dtype), x, x**2, x**3, x**4])

        # xtra constraint: value at 0 = value x(end)
        extras = np.array([[0, xstart-xend, xstart**2-xend**2, xstart**3-xend**3, xstart**4-xend**4]]).T
        
        matx = np.hstack((matx,extras)).T
        y = np.matrix(np.hstack((y0,0)))

        matx[-1,:] = matx[-1,:]*weight

        coeff = np.linalg.inv(matx.T*matx)*matx.T*y.T
        coeff = coeff.A.squeeze()[::-1]

        self.results = coeff
        return coeff

    def model(self, x):
        return get_model(x, self.results)


    
class WobbleFitSin(WobbleFit):
    """ Fit the wobble with a sinus 

    The period is known, 
    the x array must be wrapped from 0 to 2*pi over a period
    """                



    def fit(self):
        x, y = self.x, self.y

        guess_mean = np.mean(y)
        guess_std = 3*np.std(y)/(2**0.5)
        guess_phase = 0

        optimize_func = lambda c: c[0]*np.sin(x+c[1]) + c[2] - y    
        est_std, est_phase, est_mean = leastsq(optimize_func, [guess_std, guess_phase, guess_mean])[0]

        self.results = est_std, est_phase, est_mean
        return self.results

    def model(self, x):
        est_std, est_phase, est_mean = self.results
        return est_std*np.sin(x+est_phase) + est_mean            


class WobbleFitSin2(WobbleFit):
    """ Fit the wobble with a sinus 

    The period is known, 
    the x array must be wrapped from 0 to 2*pi over a period
    """                
    
    def fit(self):
        x, y = self.x, self.y

        guess_mean = np.mean(y)
        guess_std = 3*np.std(y)/(2**0.5)
        guess_phase = 0

        optimize_func = lambda c: c[0]*np.sin(x+c[1]) + c[2] - y    
        est_std, est_phase, est_mean = leastsq(optimize_func, [guess_std, guess_phase, guess_mean])[0]

        self.results = est_std, est_phase, est_mean
        return self.results

    def model(self, x):
        est_std, est_phase, est_mean = self.results
        return est_std*np.sin(x+est_phase) + est_mean            



