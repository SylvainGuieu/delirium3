""" 
    Define all the plot for delirium, dl, supports object.

    The plots are defined in a separate file so the script can still
    run if matplotlib library is not available. 
"""

from matplotlib.pylab import plt
import matplotlib.dates as mdates
import datetime
from .supports import Supports
from .delirium import Delirium
from .rail import Rail
from .sensors import Sensors, Inclinometer, Fogale
from .carriage import Carriage
from .dl import DelayLineStates,DelayLineState, DelayLineHysteresis, DelayLine     
import os
from .parameters import parameters
from . import computing
from .log import Log
import numpy as np

try:
    unicode
except NameError:
    basestring = (str, bytes)

log = Log(context=("DELIRIUM", "PLOT"))

# figure saved format
figure_format = "png"

def id2title(id):
    return "DL{0} - {1} - {2}".format(id[0], id[1], ["R","D"][bool(id[2])])

def id2file(id):
    return "dl{0}_{1}_{2}".format(id[0],id[1],  ["R","D"][bool(id[2])])


def getfigure(figure=None, clear=False, show=False, save=False):
    save = save
    if not isinstance(figure, plt.Figure):
        figure = plt.figure(figure)
    if clear:
        figure.clear()
    if show or save:            
        def finalyse(name, save=save, show=show):
            if show:
                figure.show()
                figure.canvas.draw()
            fout = ""
            if save: 
                if isinstance(save, basestring):
                    fout = save
                    f = save
                else:
                    fout = getattr(save, "name", "")
                    #fout = name + "." + figure_format   
                    f = save
                figure.savefig(f)
                log.notice("Figure %s saved to %s"%(name, fout))
            return fout    

    else:
        def finalyse(name):
            return name                
    return figure, finalyse        

def getaxes(axes=None, figure=None, clear=False, fclear=False, show=False, save=False):
    if axes is None:
        figure, finalyse = getfigure(figure, fclear, show, save)
        axes = figure.add_subplot(111)

    elif not isinstance(axes,(plt.Subplot, plt.Axes)):
        figure, finalyse = getfigure(figure,fclear, show, save)
        if isinstance(axes, tuple):
            axes = figure.add_subplot(*axes)  
        else:    
            axes = figure.add_subplot(axes)                    
    else:
        def finalyse(name):
            axes.figure.show()
            axes.figure.canvas.draw()

    if clear:
        axes.clear()   
    return axes, finalyse

def getaxesgrid(grid, axes=None, figure=None, clear=False, fclear=False, show=False,save=False,  N=None):
    if N is None:
        N = grid[0]*grid[1]
    figure, finalyse = getfigure(figure, fclear, show, save)

    if axes is None:        
        axess = [getaxes(grid+(i,), figure, clear)[0] for i in range(1,N+1)]
    else:
        if not hasattr(axes, "__len__") or len(axes)!=N:
            raise ValueError("need %d axis for this graph"%N)               
        axess = [getaxes(a, figure, clear)[0] for a in axes]
    return axess, finalyse   


class PlotCollection(object):
    """ A collection of plot """

    def __get__(self, data, cl=None):
        """ This is a trick to use e.g. :  supports.plot.correction() """
        sid = id(self)
        if id(self) in getattr(data, "__plots__", {}):
            return data.__plots__[sid]
        if not hasattr(data, "__plots__"):
            data.__plots__ = {}
        new = self.__class__()
        new.__dict__.update(self.__dict__)  
        data.__plots__[sid] = new
        new.data = data
        return new

    def figure(self, figure=None, clear=False, show=False, save=False):
        return getfigure(figure=figure, clear=clear, show=show, save=save)
     
    def axes(self, axes=None, figure=None, clear=False, fclear=False, show=False, save=False):
        return getaxes(axes=axes, figure=figure, clear=clear, fclear=fclear, show=show, save=save)
    
    def axesgrid(self, grid, axes=None, figure=None, clear=False, fclear=False, show=False,save=False,  N=None):
        return getaxesgrid(grid, axes=axes, figure=figure, clear=clear, fclear=fclear, show=show, save=save)        

class PlotSupports(PlotCollection):
    
    def corrections(self, direction, axes=None, figure=None, aclear=False, fclear=False, show=False,save=False, unit="micron", treshold=None, **kwargs):
        """ plot support correction with bars

        Parameters
        ----------
        direction :  string
            must be 'H' or 'V' for vertical or horizontal 
        unit :  string
            correction unit 'micron', 'mm' or 'm'    
        treshold : value or None(for default)
            plot the treshold and -treshold limits 

        **kwargs : other option for the bar plot

        Outputs:
        axes : thes axes                        
        """
        axes, finalyse = self.axes(axes,figure,aclear,fclear,show, save)
        direction = direction.upper()

        supports = self.data
        c = supports.get_corrections( supports.supports, unit=unit, treshold=0.0)
        treshold  = supports.correction_treshold if treshold is None else treshold 
        treshold *= supports.get_scale(unit)                

        axes.bar(c['support'], c[direction], **kwargs)
        axes.axhline(  treshold , color="red")
        axes.axhline( -treshold, color="red")
        axes.set_ylabel("Correction in %s direction [%s]"%(direction,unit))
        axes.set_xlabel("Support number")
        axes.set_title(id2title(supports.get_id()))

        return finalyse("corrections_"+direction+"_"+id2file(supports.get_id())), (axes,)


    
    def deformations(self, direction, axes=None, figure=None, aclear=False, fclear=False,show=False,save=False, unit="micron",  treshold=None, **kwargs):
        axes, finalyse = self.axes(axes,figure, aclear, fclear,show,save)
        direction = direction.upper()
        if direction in ["Z","V", "VERTICAL", "Vcorrection"]:
            drail, dsupport = "z", "Vcorrection"
        elif direction in ["Y","H", "HORIZONTAL", "Hcorrection"]:
            drail, dsupport = "y", "Hcorrection"    

        supports = self.data
        rail_data = supports.rail_data
        scale = supports.get_scale(unit)
        treshold = supports.correction_treshold if treshold is None else treshold  

        supports_opl = supports.get_opl()

        supval = supports.get(dsupport)
        axes.bar(supports_opl, supval*scale, label="Corrections")
        test = abs(supval)>treshold
        if any(test):
            axes.bar(supports_opl[test], supval[test]*scale, color="red")
                            
        axes.plot(rail_data["x"], rail_data[drail]*scale, label="Deformation", **kwargs)

        axes.axhline(  treshold*scale, color="red")
        axes.axhline( -treshold*scale, color="red")

        lim = max(treshold, 0.007)*scale*4  

        axes.set_ylabel("rail y [%s]"%unit)    
        axes.set_xlabel("opl [m]")
        axes.set_ylim( -lim, lim)
        axes.set_title(id2title(supports.get_id()))
        return finalyse("deformations_"+direction+"_"+id2file(supports.get_id())), (axes,)
        #return finalyse(DeliriumFile("deformations_"+direction+"_"+id2file(supports.get_id()))), (axes,)



class PlotParam:
    """ some generic function to plot simple parameter """
    def param(self, key, order=None, fitrange=None, removeLowOrder=False, filterWobble=False, wobbleFitrange=None, wobbleOrder=None, period=None, arcsec=False,  axes=None, figure=None, aclear=False, 
                   fclear=False, show=False,save=False, opls=None, indexes=None, **kwargs):
        """ Plot the parameter value vs opl

        Parameters
        ----------
        key : string
            the parameter name : 'ycrt', 'zctr', 'inc', 'horizontal', 'vertical', 'phi', 'theta', 'psi', 'y', 'z'

        order : float or None
            if not None alter the default low order removal

        fitrange : string, 2xtuple or None
            if not None change the default range for the low order fit removal.
            Can be a (min, max) tuple in opl [m]
            Or a string that define a range for each delay line : 'good', 'full', 'conservative'
        
        removeLowOrder: bool, optional
            Remove the lowest order polynome. The order is pre-defined for each parameters (see .params attribute)
            Of course, has no effect if key is 'opl' or 'x'

        
        filterWobble: bool, optional
            If True, the wobble is filtered. The wobble period is pre-defined for each measurement (see .params attribute)
            This has no effect on 'opl'
        
        wobbleFitrange :  string, 2xtuple or None
            if not change the default range for the wobble function when 
            removing the low order   

        wobbleOrder : int or None
            if not None, change the default order to fit the low order 
            before fitting the wobble  
        
        period : string or float or None, optional
            if not None alter the default period of requested params for wobble fit.
            Period can be a string representing a physical thing : 'wheel', 'support'
            or a float in opl [m] unit.
            This has no effect if filterWobble is False     

        indexes : array like of int
            Indexes of returned values 
            is indexes is not None, opls must be None

        opls : 2 tuple or float
            if float return the value for the closest opl
            if 2 tuple (min,max) return the values in between or 
            equal min and max.
            if opls is not None, indexes must be None

        arcsec : bool
            transform radiant to arcsec, only works on angle parameters    
        
        **kwargs : other option for the bar plot                
        """
        axes, finalyse = self.axes(axes,figure, aclear, fclear,show, save)
        
        container = self.data

        y, extras = container.get(key, arcsec=arcsec, removeLowOrder=removeLowOrder, filterWobble=filterWobble, 
                                       wobbleFitrange=wobbleFitrange, wobbleOrder=wobbleOrder, order=order, 
                                       fitrange=fitrange, extras=True, opls=opls, indexes=indexes)
        x = container.get("opl",  opls=opls, indexes=indexes)
        oplparam = container.params.get("opl")
        param = container.params.get(key)

        axes.plot( x, y, 'b+', **kwargs)
        axes.set( xlabel="{0.label} [{0.unit}]".format(oplparam), ylabel="{0.label} [{unit}]".format(param, unit=extras.get("unit",param.unit))  )
        
        id = ""
        if hasattr(container, "get_id"):       
            id = id2title(container.get_id())
            axes.set_title(id)

        return finalyse("parameter_"+key+"_"+id), (axes,)


class PlotFit:
    """ some generic function to plot fit results """

    def low_order_fit(self, key, order=None, fitrange=None, arcsec=False,  axes=None, figure=None, aclear=False, fclear=False, show=False,save=False, **kwargs):
        """ Plot the low order fit for a parameter and its residuals

        Parameters
        ----------
        key : string
            the parameter name : 'ycrt', 'zctr', 'inc', 'horizontal', 'vertical', 'phi', 'theta', 'psi', 'y', 'z'

        order : float or None
            if not None alter the default low order for the requested param.

        fitrange : string, 2xtuple or None
            if not None change the default range for the low order fit.
            Can be a (min, max) tuple in opl [m]
            Or a string that define a range for each delay line : 'good', 'full', 'conservative'

        arcsec : bool
            transform radiant to arcsec, only works on angle parameters    

        axes : list like  or None
            must contain 2 axes or 2 axes position e.g : [211,212]   
        **kwargs : other option for the bar plot                
        """
        (axes1,axes2), finalyse = self.axesgrid((2,1), axes, figure, aclear, fclear, show, save)

        kopl = 'x' if key in ['y','z'] else 'opl'
        obj = self.data
        param = obj.params.get(key)
        oplparam = obj.params.get(kopl)
        x = obj.get(kopl)

        yfiltered, extras = obj.get(key, arcsec=arcsec, removeLowOrder=True, order=order, fitrange=fitrange, extras=True)
        yraw = obj.get(key, arcsec=arcsec)

        axes1.plot( x, yraw, 'b+', **kwargs)
        axes1.set( xlabel="{0.label} [{0.unit}]".format(oplparam), ylabel="{0.label} [{unit}]".format(param, unit=extras.get("unit",param.unit))  )
        try:
            pol = extras['loworder_polynome']
        except KeyError:
            log.warning("Cannot find low order fit polynome for parameter '%s'"%key)    
        else:
            xm = np.linspace(x.min(), x.max(), 200)
            axes1.plot( xm, computing.get_model(xm, pol), "r-", label=", ".join("%2.2E"%p for p in pol))
            axes1.legend(fontsize='small')
        axes1.set_title(id2title(obj.get_id()))    
            
        axes2.plot( x, yfiltered, 'b+', **kwargs)
        axes2.set( xlabel="{0.label} [{0.unit}]".format(oplparam), ylabel="{0.label} [{unit}] residuals".format(param, unit=extras.get("unit",param.unit)))


        rge = extras.get('loworder_range', None)
        if rge:
            for x,a in zip(rge*2, [axes1]*2+[axes2]*2):
                a.axvline(x, linestyle="dashed", color="gray")
        
        ## share the x axis   
        axes1.get_shared_x_axes().join( axes1, axes2)   
        return finalyse("fit_"+key+"_"+id2file(obj.get_id())), (axes1,axes2)
                    

    def wobble_fit(self, key,  period=None, removeLowOrder=True, axes=None, fitrange=None, order=None, arcsec=True, figure=None, aclear=False, fclear=False, show=False,save=False, **kwargs):
        """ plot the wobble fit and its effect on data


        Parameters
        ----------
        key : string
            the parameter name : 'ycrt', 'zctr', 'inc', 'horizontal', 'vertical', 'phi', 'theta', 'psi', 'y', 'z'

        period : string, float or None, optional
            If None the default period for that parameter is used. 
            period can be a string representing a physical thing : 'wheel', 'support'
            or a float in opl [m] unit.

        removeLowOrder : bool, optional
            If True (default) remove the low order in the second plot

        axes : list like  or None
            must contain 2 axes or 2 axes position e.g : [211,212]  
                
        **kwargs :
            other regular option
        
        Outputs
        -------
            axes1 : plt.Axes
                axes where the fit is represented
            axes2 : plt.Axes
                filtered and non filtered measurement plot.                            
        """
        (axes1,axes2), finalyse = self.axesgrid((2,1), axes, figure, aclear, fclear, show, save)

        obj = self.data

        yfiltered, extras = obj.get(key, filterWobble=True, period=period, removeLowOrder=True, 
            wobbleFitrange=fitrange,  arcsec=arcsec, extras=True, wobbleOrder=order)
        
        fit = extras['wobble_fit']


        if not removeLowOrder:
            yfiltered2 = obj.get(key, filterWobble=True, removeLowOrder=False, arcsec=arcsec)
        else:
            yfiltered2 = yfiltered

        yraw = obj.get(key, removeLowOrder=removeLowOrder, arcsec=arcsec)



        param = obj.params.get(key)
        period = param.period if period is None else period

        x = obj.get_opl( period )
        axes1.figure.subplots_adjust(hspace=0.3)
        axes1.plot(x/(2*np.pi), yraw, "k+", **kwargs)
        axes1.plot(fit.x/(2*np.pi), fit.y, "b+", label="used in fit")
        axes1.set_ylim( np.min(fit.y), np.max(fit.y) )

        speriod = "'%s'=%.2f m"%(period, obj.parse_period(period)) if isinstance(period, basestring) else "%.2fm"%period
        axes1.set_xlabel("wrapped opl with period of %s"%(speriod))
        axes1.set_ylabel("%s [%s]"%(key, extras['unit']))
        axes1.legend(loc="lower right", fontsize="small")        
        
        axes1.set_title(id2title(obj.get_id())+"- Wobble Amplitude = {amplitude:2.2e} {unit}".format(*obj.get_id(), amplitude=extras.get('wobble_amplitude', 0.0), unit=extras['unit']))
        try:
            fit = extras['wobble_fit']
        except:
            log.warning("Cannot found fit wobble information for parameter %s"%key)
        else:            
            x = np.linspace(0, 2*np.pi,500)
            axes1.plot(x/(2*np.pi), fit.model(x), color="red")

        
        x = obj.get_opl()
        axes2.plot(x, yraw      , color="red", label="raw")
        axes2.plot(x, yfiltered2, color="blue", label="filtered")
        axes2.set_ylabel("%s [%s]"%(key, extras['unit']))
        axes2.set_xlabel("opl [m]")        
        axes2.legend(fontsize="small")

        rge = extras.get("wobble_range", None)
        if rge:            
            for x in rge:
                axes2.axvline(x, linestyle="dashed", color="gray")

        return finalyse("wobble_"+key+"_"+id2file(obj.get_id())), (axes1, axes2)




class PlotRail(PlotCollection, PlotFit,PlotParam):



    def deformations(self, direction, treshold=None, axes=None, figure=None, aclear=False, fclear=False, show=False,save=False, unit="micron", **kwargs):
        axes, finalyse = self.axes(axes,figure, aclear, fclear,show, save)
        direction = direction.upper()
        if direction in ["Z","V", "VERTICAL", "Vcorrection"]:
            drail, dsupport = "z", "Vcorrection"
        elif direction in ["Y","H", "HORIZONTAL", "Hcorrection"]:
            drail, dsupport = "y", "Hcorrection"    

        rail = self.data
        supports = rail.supports
        supports_opl = supports.get_opl()

        treshold = supports.correction_treshold if treshold is None else treshold  
        scale = supports.get_scale(unit)

        supval = supports.get(dsupport)
        axes.bar(supports_opl, supval*scale, label="%s corrections"%dsupport)
        test = abs(supval)>treshold
        if any(test):
            axes.bar(supports_opl[test], supval[test]*scale, color="red")
        

        axes.plot(rail.get_opl(), rail.get(drail)*scale, label="%s deformations"%drail)

        axes.axhline(  treshold*scale, color="red")
        axes.axhline( -treshold*scale, color="red")
        lim = max(treshold, 0.007)*scale*4  


        axes.legend(fontsize='small')

        axes.set_ylabel("rail %s [%s]"%(drail, unit))    
        axes.set_xlabel("opl [m]")
        axes.set_ylim( -lim, lim)
        axes.set_title(id2title(rail.get_id()))
        return finalyse("deformation_"+direction+"_"+id2file(rail.get_id())), (axes,)

class PlotDelirium(PlotCollection, PlotFit,PlotParam):
    pass

class PlotSensors(PlotCollection, PlotFit,PlotParam):
    def roll_correction(self, key, axes=None, figure=None, aclear=False, fclear=False, show=False,save=False,**kwargs):        
        (axes1,axes2), finalyse = self.axesgrid((2,1), axes, figure, aclear, fclear,show, save)

        sensors = self.data
        opl = sensors.get_opl()
        axes1.plot(opl, sensors.delirium.get(key),color="b", label="raw")
        axes1.plot(opl, sensors.get(key),color="r", label="corrected")
        
        param = sensors.params.get(key)

        axes1.set_ylabel("%s [%s]"%(key, param.unit))                
        axes1.set_xlabel("opl [m]")
        axes1.legend()
        axes1.set_title(id2title(sensors.get_id()))

        axes2.plot(opl, sensors.delirium.get(key, removeLowOrder=True), color="b", label="raw")
        axes2.plot(opl, sensors.get(key, removeLowOrder=True), color="r", label="corrected")
            
        axes2.set_ylabel("%s (order %s removed) [%s]"%(key, param.order, param.unit)) 
        axes2.set_xlabel("opl [m]")
        return finalyse("correction_"+key+"_"+id2file(sensors.get_id())), (axes1,axes2)

class PlotCarriage(PlotCollection, PlotFit, PlotParam):
    def flatness(self, N=15, axes=None, figure=None, aclear=False, fclear=False, show=False,save=False,**kwargs):
        (axtheta,axpsi,axZtheta,axZpsi), finalyse = self.axesgrid((1,4), axes, figure, aclear, fclear,show, save)
                                
        carriage = self.data
        order, angle = carriage.get_flat_points(N)
        
        #order.sort() # sort in opl
        theta = (carriage.get("theta", arcsec=True, order=0, removeLowOrder=True))
        psi   = (carriage.get("psi", arcsec=True, order=0, removeLowOrder=True))
        opl = carriage.get("opl")


        ztheta = theta[order]
        zpsi = psi[order]
        zopl  = opl[order]        

        i = np.arange(len(zopl))

        alpha = 0.3

        axtheta.barh(opl,theta, align='center', color='gray', zorder=10)
        axtheta.barh(zopl, [max(theta)]*len(zopl), align='center', color='red', zorder=9, alpha=alpha)
        axtheta.barh(zopl, [min(theta)]*len(zopl), align='center', color='red', zorder=9, alpha=alpha)

        axtheta.set(xlabel='theta [arcsec]', ylabel="opl [m]",
                    title=id2title(carriage.get_id())+" Profile"
                    )
        axpsi.barh(opl, psi, align='center', color='gray', zorder=10)

        axpsi.barh(zopl, [max(psi)]*len(zopl), align='center', color='red', zorder=9, alpha=alpha)
        axpsi.barh(zopl, [min(psi)]*len(zopl), align='center', color='red', zorder=9, alpha=alpha)

        axpsi.set(xlabel='psi [arcsec]')




        #axtheta.invert_xaxis()
        #axtheta.set(yticks=i, yticklabels=["%4.2f"%o for o in opl])
        #axtheta.yaxis.set_major_locator(plt.NullLocator())
                
        #axtheta.yaxis.tick_right()


        axZtheta.barh(i, ztheta, align='center', color='gray', zorder=10)
        axZtheta.set(xlabel='theta [arcsec]', 
                    title="%d Best"%N                    
            )
        axZpsi.barh(i, zpsi, align='center', color='gray', zorder=10)
        axZpsi.set(xlabel='psi [arcsec]')

        #axZtheta.invert_xaxis()
        axZtheta.set(yticks=list(i)+[N], yticklabels=["%4.2f"%o for o in zopl]+["OPL"])
        axZtheta.set_ylim(0,N)
        axZpsi.set_ylim(0,N)

        axZpsi.yaxis.set_major_locator(plt.NullLocator())
        axZpsi.set(yticks=[0,N], yticklabels=["best", "worst"])
        axZpsi.yaxis.tick_right()
        axZtheta.yaxis.tick_right()

        for ax in [axtheta,axpsi, axZtheta, axZpsi]:
            ax.margins(0.03)
            ax.grid(True)

        fig = axtheta.figure    
        #fig.tight_layout()
        fig.subplots_adjust(wspace=0.4) 
        fig.set_size_inches(14,5.0)   

        return finalyse("flatness_"+id2file(carriage.get_id())), (axtheta,axpsi)


class PlotHysteresis(PlotCollection, PlotFit, PlotParam):
    def histeresis(self, key, axes=None, figure=None, aclear=False, fclear=False, show=False, save=False, arcsec=True, 
                    removeLowOrder=False, filterWobble=False,
                    **kwargs):        
        (axes1,axes2), finalyse = self.axesgrid((2,1), axes, figure, aclear, fclear, show, save)

        if key not in ["theta", "psi"]:
            raise ValueError("Key must be theta or psi")
        

        histeresis = self.data
        rge = histeresis.parse_range(histeresis.opl_histeresis_fit_range)

        unit = "arcsec" if arcsec else "rad" 
        scale = 1./(np.pi/180/3600.) if arcsec else 1

        opl = histeresis.get("opl")
        direct = histeresis.get(key+"_direct", arcsec=arcsec, removeLowOrder=removeLowOrder, filterWobble=filterWobble)
        reverse = histeresis.get(key+"_reverse", arcsec=arcsec, removeLowOrder=removeLowOrder, filterWobble=filterWobble)
        diff = histeresis.get(key+"_diff", arcsec=arcsec)


        axes1.plot( opl, direct ,  label="direct", **dict(dict(color="blue"), **kwargs))
        axes1.plot( opl, reverse, label="reverse", **dict(dict(color="red"), **kwargs))

        ind = histeresis.flat_point_indexes
        axes1.plot( opl[ind], direct[ind], **dict(dict(color="blue", marker="o", linestyle="none"), **kwargs))
        axes1.plot( opl[ind], reverse[ind], **dict(dict(color="red", marker="o", linestyle="none"), **kwargs))


        axes1.set_ylabel("%s [%s]"%(key, unit))
        

        axes1.legend(loc='lower right', fontsize='small')   
        axes2.plot( opl, diff, **dict(dict(color="blue"),**kwargs))
        axes2.plot( opl[ind], diff[ind], **dict(dict(color="blue", marker="o", linestyle="none"), **kwargs))



        axes2.set_ylabel("%s differences [%s]"%(key, unit))
        axes2.set_xlabel("opl [m]")

        delta = getattr(histeresis, key+"_diff_mean",0.0)*scale
        std = getattr(histeresis, key+"_diff_std", 0.0)*scale

        for i in ind:
            axes2.text(opl[i], diff[i]+([std,-std][i%2]*3), "%4.2f"%opl[i], horizontalalignment=['left','right'][i%2], rotation=45)
            #axes2.annotate("%4.2f"%opl[i], (opl[i],diff[i]), (opl[i],diff[i]+([std,-std][i%2]*3)) , rotation=45)    

        axes1.set_title(id2title(histeresis.get_id())+" Hysteresis %.2g+-%.2g  %s"%(delta, std, unit))
        axes2.axhline(delta, linestyle="dashed", color='k')
        for sign in [1,-1]:        
            axes2.axhline( delta+std*sign, linestyle="dashdot" , color='k')

        for v,a in zip(rge*2, [axes1]*2+[axes2]*2):
            a.axvline(v,color="k", linestyle="dashed")


        return finalyse("hysteresis_"+key+"_"+id2file(histeresis.direct.get_id())), (axes1,axes2)


class PlotFogale(PlotCollection, PlotFit, PlotParam):
    pass

class PlotInclinometer(PlotCollection, PlotFit, PlotParam):
    pass
        

class PlotDelayLineState(PlotCollection, PlotFit, PlotParam):
    pass


########################################
#
# Add the plot to all classes 
#
########################################

Supports.plot = PlotSupports()
Delirium.plot = PlotDelirium()
Rail.plot     = PlotRail()
Sensors.plot  = PlotSensors()
Fogale.plot   = PlotFogale()
Inclinometer.plot = PlotInclinometer()
Carriage.plot = PlotCarriage()
DelayLineHysteresis.plot = PlotHysteresis()
DelayLineState.plot = PlotDelayLineState()



########################################
#
# misc plots
#
########################################

def plot_history(data, key, dlnum=0, kind="", unit="arcsec", axes=None, figure=None, aclear=False, 
                 fclear=False, show=False,save=False, **kwargs):
    
    axes, finalyse = getaxes(axes,figure, aclear, fclear,show, save)

    axes.set_ylabel("%s %s [%s]"%(kind,key,unit))

    if not len(data):
        axes.text(0.5,0.5, "No Data")
        axes.set_xlim(0,1)
        axes.set_ylim(0.25,0.75)
        axes.set_title("DL %d"%(dlnum))

        return finalyse("monitoring_%s_%d_%s.png"%(kind,dlnum,key)), (axes,)

    years = mdates.YearLocator()   # every year
    months = mdates.MonthLocator()  # every month
    yearsFmt = mdates.DateFormatter('%Y')
    monthFmt = mdates.DateFormatter('%m')


    date = data['date']
    y = data[key]
    axes.plot(date, y, "k+", **kwargs)

    axes.xaxis.set_major_locator(years)
    axes.xaxis.set_major_formatter(yearsFmt)
    axes.xaxis.set_minor_locator(months)
    #axes.xaxis.set_minor_formatter(monthFmt)
    datemin = datetime.date(date.min().year, 1, 1)
    datemax = datetime.date(date.max().year + 1, 1, 1)
    axes.set_xlim(datemin, datemax)
    
    axes.set_title("DL %d %s -> %s"%(dlnum, date.min(), date.max()))

    # rotates and right aligns the x labels, and moves the bottom of the
    # axes up to make room for them
    axes.figure.autofmt_xdate()
    axes.grid(True)
    return finalyse("monitoring_%s_%d_%s.png"%(kind,dlnum,key)), (axes,)

def plot_temperature(data, key, dlnum=0, kind="", unit="arcsec", axes=None, figure=None, aclear=False, 
                 fclear=False, show=False,save=False, **kwargs):
    
    axes, finalyse = getaxes(axes,figure, aclear, fclear,show, save)
    axes.set_ylabel("%s %s [%s]"%(kind,key, unit))
    axes.set_xlabel("temperature [deg]")

    if not len(data):
        axes.text(0.5,0.5, "No Data")
        axes.set_xlim(0,1)
        axes.set_ylim(0.25,0.75)
        axes.set_title("DL %d"%(dlnum))
        return finalyse("temperature_%s_%d_%s.png"%(kind,dlnum,key)), (axes,)

    axes.plot(data['temp'], data[key], "k+", **kwargs)    
    
    axes.set_title("DL %d  %s -> %s"%(dlnum, data['date'].min(), data['date'].max()))

    return finalyse("temperature_%s_%d_%s.png"%(kind,dlnum,key)), (axes,)


