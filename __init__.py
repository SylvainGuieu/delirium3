""" delirium is a package to explore the delirium data, to compute
the rail deformation and corrections to apply automaticaly or in a more 
interactive python shell.

See the Readme file for more info.

Quick start:
>>> import delirium 
The level of verbose can be set with:
>>> delirium.set_verbose_type('EW')  #Errors and Warnings only 

>>> delirium.processing.run() # to process today delirium
>>> delirium.processing.run("2016-08-08")

## If you want to create the results (and webpage) on your computer:
>>> delirium.processing.setup(websiteroot="/path/to/my/local/directory")
>>> delirium.processing.run("2015-12-24", keep=True) 
# keep=True allows the result directory to stay on disk and never been cleaned up

>>> dl1 = delirium.open_dl(2, "2016-08-08")
>>> dl1.carriage.get("theta") # to get the theta value 
>>> dl1.supports.get("Hcorrection) # to get the horizontal corrections

"""
from .supports import Supports
from .delirium import Delirium, DataError, DataUtils
from .rail import Rail
from .sensors import Sensors, Fogale, Inclinometer
from .carriage import Carriage    
from .parameters import parameters
from .dl import DelayLineState, DelayLineStates, DelayLineHysteresis
from .log import Log, ERROR, WARNING, NOTICE, INFO, DATA
from . import run as processing
from . import plots

run = processing.run
reprocess_all  = processing.reprocess_all
list_all_dates = processing.list_all_dates
list_dates = processing.list_dates
run_interval = processing.run_interval

def open_dl(num, date=None, directory=None, reverse=False, file_index=-1):
    """ Return the delay line object of the last delirium found for the given date and delay line number

    If date is not given,  take the date of today.
    Parameters
    ----------
    num : int
        delay line number 
    date : string, optional
        date in the format 'yyyy-mm-dd'. If None take the today date.
    directory : string, optional
        the root directory containing the deliriums. default is 
            'ftp://utcomm@odyssey3.pl.eso.org/DELIRIUM'
    reverse : bool, optional
        If True, look for the reverse file
    file_index : int, optional
        in case of several file found gives the list index, default is -1 
        Files are ordoned in increasing date
    
    Outputs
    -------
    dl : DelayLineState object
                        
    """
    from . import io
    if directory is None:
        directory = "ftp://utcomm:Bnice2me@odyssey3.pl.eso.org/DELIRIUM"
    dldir = io.DeliriumDirectory(directory)
    files = dldir.deliriumfiles(date, "reverse" if reverse else "direct")
    if not len(files):
        raise IOError("Cannot find %s file for DL %d and  date '%s'  "%("reverse" if reverse else "", num, date))
    files = files[num]
    return DelayLineState(deliriumfile=files[file_index])


def open_histeresis(num, date=None, directory=None, file_index=-1):
    dld = open_dl(num, date=date, directory=directory, file_index=file_index)

    try:
        dlr = open_dl(num, date, directory=directory, reverse=True, file_index=file_index)
    except IOError:
        raise IOError("Cannot find the reverse file for that date %s"%date)    

    return DelayLineHysteresis(dld, dlr)


def update_monitoring(dlnum, start=None, end=None, step=1):
    from . import iomonitoring
    import numpy as np
    from datetime import datetime
    #directory = "ftp://utcomm:Bnice2me@odyssey3.pl.eso.org/DELIRIUM"
    directory = "/Users/sylvain/tmp"
    file = io.fpath(directory+"/dl%d_monitoring.fits"%dlnum)

    try:
        f = file.open("rb")
    except:
        new = iomonitoring.create_monitoring(dlnum)
        new.writeto(file.open("wb"), overwrite=True)
        f = file.open("rb")

    fh = [h.copy() for h in iomonitoring.fits.open(f)]
    if start is None:
        dates = np.array([datetime(d['year'],d['month'],d['day']) for d in fh[1].data])        
        start = dates.max()
        start = "%d-%02d-%02d"%(start.year,start.month,start.day)

    f.close()

    try:
        new = iomonitoring.update_monitoring(fh,  start, end, open_dl )
    except:
        pass
    else:
        new.writeto(file.open("wb"), overwrite=True)
    


def set_verbose_type(msgtypes):
    """ Set the verbose filter type 
    
    Parameters
    ----------
    msgtypes : 4-bites int
          1 1 1 1
          | | | |
       DATA | | |
       NOTICE | |
        WARNING |
            ERROR
    
    Examples
    --------
        set_verbose_type(1+2) will verbose only ERROR and WARNING
        set_verbose_type(4) will verbose only the notices
    """
    global Log

    if isinstance(msgtypes, basestring):
        msgtypes = sum( (s in msgtypes)*bit for s,bit in [('E',ERROR),('W',WARNING),('N',NOTICE),('D',DATA)] )
    if not isinstance(msgtypes, int):
        raise ValueError("expecting int got %s"%msgtypes)    
    Log.msgtypes = msgtypes

