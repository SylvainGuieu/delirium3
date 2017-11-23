"""
Provide function to do daily delirium task or recompute all the 
monitoring value. 

Bottom line:
> run()  # to run today delirium for all DL and log all plot in wiki
> run("2016-08-02") # run delirium of "2016-08-02"
> run("2016-08-02", [1,3]) # run delirium of "2016-08-02" for dl 1 and 3
> run("2015-01-22", plot=False, wiki=False) # execute delirium but do not plot or 
                                            # log anything

> reprocess_all(list_all_dates()) # reprocess everything that is recorded in disk 
                                  # without plotting of course but it will update the
                                  # monitoring files
"""

from .dl import DelayLineState, DelayLineHysteresis
from .log import Log, get_buffer, clear_buffer
from . import io
from .delirium import DataError
from . import plots
import time, datetime
import sys
import os
log = Log(context=("DELIRIUM",))
DEBUG = False

## number of subdiretories that will be saved
maxDays  = 8 # keep always the delirium that are younger or maxDays old 
cashSize = 8 # keep also a max number of cashSize for file that does not match 
             # the criteria above. (useful if one execute the scripte on some old delirium) 



## directory where to found the yyyy-mm/delirium_files
deliriumdir = "ftp://utcomm:Bnice2me@odyssey3.pl.eso.org/DELIRIUM"
#dldir = io.DeliriumDirectory("/Users/sylvain/Dropbox/python/delirium/data/examples/DELIRIUM/")

####
## directory for the wiki monitoring 
# with the structure 
# delirium | 
#          | data 
#                  | daily       -> daily results
#                  | monitoring  -> monitoring files
websiteroot = "ftp://vlti:boulder56@vlti.pl.eso.org//diska/web/vlti/"
#webroot = io.dpath("/Users/sylvain/tmp")


def setup(deliriumroot=deliriumdir, websiteroot=websiteroot):
    """ setup all the key path 
    
    setup is also restablishing ftp connection if time-out
    """
    global dldir, webroot, moduledir, localhtmldir, deliriumdir
    global datadir, dailydir, datesfile, indexfile, monitoringdir

    ##
    # some file needs to be copied from the module to the webpage
    moduledir, _ = os.path.split(__file__)
    localhtmldir = io.dpath(moduledir+"/html")
        
    dldir = io.DeliriumDirectory(deliriumroot)        
    webroot = io.dpath(websiteroot)
   


    deliriumdir = webroot.dpath("delirium")
    datadir = deliriumdir.dpath("data")
    dailydir, monitoringdir = datadir.dbreak("daily", "monitoring")
    
    datesfile = dailydir.fpath("dates.txt")
        
    indexfile = deliriumdir.fpath("index.html")
    indexfile.header = localhtmldir.fpath("index.html").open('r').read
    #indexfile.header = lambda : localhtmldir.fpath("index.html").read().decode()

#setup()

#print("aaa", type(indexfile.header()))
#hist_header = "%% Date  TunnelTemp(deg) Hysteresis(arcsec)\n"
#wobble_header = "%% Date  TunnelTemp(deg) WoobleTheta(arcsec) WooblePsi(arcsec)\n"



def list_all_dates():
    """ list all delirium dates available on the delirium machine """
    fls = dldir.ls("[0-9][0-9][0-9][0-9]-[0-9][0-9]/DL1*.dat")
    dates = [] 
    for fl in fls:
        d, f = os.path.split(fl)
        date = f[11:19]
        date = date[0:4]+"-"+date[4:6]+"-"+date[6:8]
        if not date in dates:
            dates.append(date)
    return dates

def list_dates(start, end=None, step=1):
    """ return a list of dates from a start and a end """
    if step<1:
        raise ValueError("step must be >=1 got %s"%step)

    if end is None:
        end = datetime.date.today().isoformat()
    start, end = (datetime.date( *(int(v) for v in s.split("-"))) for s in [start,end])

    step = datetime.timedelta(step)

    d = start
    dates = []
    while (d)<=end:
        dates.append(d.isoformat())
        d = d+step

    return dates    

        
        
def reprocess_all(dates, plot=False, wiki=False):
    """ reprocess all the delirium for the given list of dates 
    
    use list_all_dates() to get all the available dates from the
    delirium machine.

    e.i. : reprocess_all(list_all_dates())
    """
    for date in dates:
        run(date, plot=plot, wiki=wiki)


def prepare_monitoring_web_structure():
    """ 
    just used to prepare the monitoring directory instead of 
    preparing all the web structure. (used when wiki=False in run)
    """
    global webroot, mdir, ddir, todaydir, monitoringdir
    monitoringdir.build()

    todaydir = io.dpath(dailydir, "dummydir")

def web_release(date):
    """ release the delirium web  directory for beeing kept in the history"""
    keeper = io.fpath(dailydir,date+"/keep")
    if keeper.exists():
        dailydir.rmtree(date+"/keep")

def web_remove(date):
    """ remove the delirium web- directory of that date """
    dailydir.rmtree(date)

def prepare_web_structure(date, nocleanup=False, keep=False):
    """ set all the directory tree and files necessary for the webpage """
    global webroot, mdir, ddir, todaydir, monitoringdir
    nocleanup = True
    ## make the dictionary if they does not exists
    dailydir.build()
    monitoringdir.build()

    ## make the target today directory 
    todaydir = io.dpath(dailydir, date).build()
    ## if does not exists make a dummy data.json file 
    ## so it has the creation date of now 
    todaydir.fpath("data.json").build()
    if keep:
        todaydir.fpath("keep").create()

    ## list all the subdirectory 
    ## check when they have been created (mtime of data.json)     
    mtimes = []
    CANMTIME = True
    if CANMTIME:
        for subdir in (dailydir.dpath(s) for s in dailydir.ls("[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]")):
            try:
                mtime = subdir.fpath("data.json").getmtime()
            except:
                ## this appen when the data.json does not exists
                ## the diretory will not be added
                mtime = 0.0
            mtimes.append( [mtime, datetime.date( *(int(v) for v in subdir.split("-"))) ])
    # order by creation date
        mtimes.sort(key=lambda i:i[0], reverse=True)
        tokeep = []
    else:
        tokeep = [datetime.date(*(int(v) for v in subdir.split("-"))) for subdir in dailydir.ls("[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]")]
        
    toremove = []
    today = datetime.date.today()
    deltatime = datetime.timedelta(maxDays)

    ## select the directories to keep and the one to remove
    # Keep directory if :
    #    - the delirium is newer or as old as maxDays
    #    - the delirium is older it was one of the last *cashSize* directory created
    #    - the directory has a file named keep 
    for i, (mtime, date) in enumerate(mtimes):
        if mtime==0.0: # -> no data.jason file  
            if not nocleanup:
                toremove.append(date)    
            continue
        if nocleanup:          
            tokeep.append(date)
        elif i<cashSize:
            tokeep.append(date)
        elif (today-date)<=deltatime:
            tokeep.append(date)
        elif io.fpath(dailydir,str(date)+"/keep").exists():
            tokeep.append(date)
        else:
            toremove.append(date)    
    tokeep.sort()

    ##
    # write the dates file with the keepers
    with datesfile.open("w") as g:
        g.write("dates = [")
        g.write(", ".join( '"%s"'%d for d in tokeep ))
        g.write("];")

    ##
    # remove the old stuff    
    for date in toremove:
        log.notice("removing Old Directory: %s"%date)
        dailydir.rmtree(str(date))
    
    ## copy the index.html file.
    indexfile.create(clobber=True)

def run_interval(start, end=None, step=1, dls=range(1,7), plot=True, wiki=True, nocleanup=False, keep=False):
    """ run delirium from start date to end date 
    
    Parameters
    ----------
    start : string
        start date yso format 'yyyy-mm-dd'
    end : string, optional
        end date yso format 'yyyy-mm-dd'
        default is today
    step : int, optional
        day step, default is one    
    dls: list of int, optional
        A list of DL number for wich the script will be executed
    plot: bool or string, optional
        can be True or "dm": everything ploted
                "d" : only daily products are plotted 
                "m"  : only monitoring products are plotted        
        default is True
    wiki: bool, optional
        if True the plot and corection are sent to the webpage
        otherwhise correction are just printed on stdout and plots on 
        screen
        default is True
    nocleanup : bool, optional
        if True old directories of the web page will not be removed
        default if False.     
    keep : bool, optional
        if True force this delirium result to always be present on the webpage 
        it will never been cleaned-up unless web_release(date) is ran.        
    see run function for more info.
    """
    for date in list_dates(start, end):        
        run(date, dls=dls, plot=plot, wiki=wiki, nocleanup=nocleanup, keep=keep)

def run(date=None, dls=range(1,7), plot=True, wiki=True, nocleanup=False, keep=False):
    """ run the delirium in a secure way 

    Parameters
    ----------
    date: string, optional
        the iso formated date ('yyyy-mm-dd'). 
        take the today date is none
    dls: list of int, optional
        A list of DL number for wich the script will be executed
    plot: bool or string, optional
        can be True or "dm": everything ploted
                "d" : only daily products are plotted 
                "m"  : only monitoring products are plotted        
        default is True
    wiki: bool, optional
        if True the plot and corection are sent to the wiki-webpage
        otherwhise correction are just printed on stdout and plots on 
        screen
        default is True
    nocleanup : bool, optional
        if True old directories of the web page will not be removed
        default if False.     
    keep : bool, optional
        if True force this delirium result to always be present on the webpage 
        it will never been cleaned-up unless web_release(date) is ran.
    """
    global todaydir, monitoringdir

    if date is None:
        date = datetime.date.today().isoformat()

    if webroot:
        if wiki:
            prepare_web_structure(date,nocleanup=nocleanup)
        else:
            prepare_monitoring_web_structure()    
    
    if wiki:
        products = io.Product(todaydir.fpath("data.json"), date, 
                               dirs={"d":"data/daily/"+date, #daily sub-directory
                                     "m":"data/monitoring"   # monitoring sub-directory
                                     })            
    else:
        products = None

    if plot is True:
        plot = "dm"
    elif not plot:
        plot = ""          
    
    ## get a list of direct files
    fdirect  = dldir.deliriumfiles(date, "direct")
    ## get a list of reverse file if any 
    freverse = dldir.deliriumfiles(date, "reverse")

    for num in dls:
        context = ("DELIRIUM", "DL%d"%num)
        clear_buffer(context[1])
        log = Log(context=context)

        ## look at direct file
        fls = fdirect.get(num, [])        
        if len(fls)>1:
            ## take the last file (they are laready sorted from older to most recent)
            file = fls[-1]
            log.warning("%d files found taking the most recent one"%(len(fls)))
        elif not fls:
            log.error("No delirium file found for dl %d"%num)  
            if "m" in plot:
                ## plot anyway the monitoring
                run_monitoring_plot_dl(num, monitoringdir, products, log)
            continue
        else:
            file = fls[-1]

        dld = run_init_dl(file, num, log)           
        run_daily_correction_dl(dld, todaydir, products, log)
        run_monitoring_dl(dld, monitoringdir, log)
        if "d" in plot:
            run_daily_plot_dl(dld, todaydir, products, log)
        if "m" in plot:    
            run_monitoring_plot_dl(num, monitoringdir, products, log)

        ############
        #
        #  Check is there is a reverse file.
        #    

        flr = freverse.get(num, [])
        if len(flr)>1:
            ## take the last file (already sorted from last to newer)
            file_r = flr[-1]
            log.notice("%d files found taking the last one"%(len(flr)))
        elif not flr:
            log.warning("No reverse delirium file found for dl %d"%num)
            if "m" in plot:
                ## plot anyway the monitoring
                run_monitoring_plot_hysteresis(num, monitoringdir, products, log)                
            continue
        else:
            file_r = flr[-1]    

        dlr = run_init_dl(file_r, num, log)         
        hyst = run_init_hysteresis( dld, dlr, log)

        run_monitoring_hysteresis(hyst, monitoringdir, log)
        if "d" in plot:
            run_daily_plot_hysteresis(hyst, todaydir, products, log)
        if "m" in plot:            
            run_monitoring_plot_hysteresis(num, monitoringdir, products, log)
                
    ##
    # add the data.json in the wiki
    if products:  
        #products.flushjs() 
        products.flushjson()      
        

class dummy(object):
    pass


def run_init_dl(file, num, log):
    ###
    # The correction file to write 
    cor_file = file.replace_ext('_CORR2.txt') # correction text
    context = log.context
    
    if DEBUG:
        dld = DelayLineState(deliriumfile=file)
    else:    
        try:                       
            dld = DelayLineState(deliriumfile=file)
        except:
            log.error("cannot load data of '%s'. Is data corrupted ?"%(file))          
            write_failure(context, file, [cor_file])
            dld = dummy()
            dld.cor_file = cor_file 
            dld.num = num
            dld.status = False
            return dld
        else:
            dld.cor_file = cor_file    
           
    dld.status = True
    return dld

def run_init_hysteresis(dld, dlr, log):
    if DEBUG:
        hyst = DelayLineHysteresis(dld, dlr)
    else:    
        try:
            hyst = DelayLineHysteresis(dld, dlr)
        except Exception as e:
            hyst = dummy()
            hyst.status = False
            log.error("Cannot compute hysteresis got error: '%s'"%(e))
            return hyst
    hyst.status = True
    return hyst          

def run_daily_correction_dl(dld, todaydir, products, log):
    """ put the daily correction file in delirium machine and website """

    num = dld.num
    cor_file = dld.cor_file
    context = log.context
    try:
        file = dld.delirium.fpath
    except:
        file = "unknown"
    if not dld.status:
        log.error("cannot compute corrections. Is data corrupted ?")
        write_failure(context, file, [cor_file]) 
    else:    
        if DEBUG:
            dld.rail.supports.log_corrections(filename=cor_file) 
        else:    
            try:   
                dld.rail.supports.log_corrections(filename=cor_file) 
            except DataError:                
                write_failure(context, file,  [cor_file])
                dld.status = False
            except:
                log.error("cannot compute corrections. Is data corrupted ?")
                write_failure(context, file, [cor_file]) 
                dld.status = False            

    if not products:
        return 

    num = dld.num    
    web_cor_file = todaydir.fpath("DL%d_last_CORR.txt"%num)
    web_cor_file.open("w").write(cor_file.open("r").read())
    #cor_file.copy(web_cor_file)
    log.notice("Correction file copied to %s"%web_cor_file) 
    products.add( num, "Corrections", "d", "txtfile",  web_cor_file.filename)

    web_summary_file = todaydir.fpath("DL%d_summary.html"%num)
    if dld.status: 
        with web_summary_file.open("w") as g:
            dld.write_summary_html(g.write)    
    else:
        with web_summary_file.open("w") as g:
            g.write("Delirium failed")
    products.add(num, "Summary", "d", "htmlfile", web_summary_file.filename)        
            

def run_monitoring_dl(dld, monitoringdir, log):
    """ update the monitoring files """ 
    if not dld.status:
        return ;
    num = dld.num    
    ###
    # update monitoring files 
    ## get the wobble monitoring file 
    wblfile = io.WobbleLog(monitoringdir.fpath("DL%d_wobble.txt"%num)).build()
    ## add the wobble amplitude results at the end of the file
    wblfile.add_from_dl(dld)
    ## get the wobble monitoring file 
    correctionfile = io.CorrectionLog(monitoringdir.fpath("DL%d_correction.txt"%num)).build()
    ## add the wobble amplitude results at the end
    correctionfile.add_from_dl(dld)
        
    return dld

def run_monitoring_hysteresis(hyst, monitoringdir, log):
    """ update the monitoring hysteresis files """ 
    if not hyst.status:
        return ;
    num = hyst.direct.num  
     ## get the Hysteresis  monitoring file 
    hystfile = io.HysteresisLog(monitoringdir.fpath("DL%d_hysteresis.txt"%num)).build()
    ## add the wobble amplitude results 
    hystfile.add_from_hysteresis(hyst)
    

def build_monitoring_fits_file(dlnum):
    pass





def run_daily_plot_dl(dld, todaydir, products, log):
    """ 
    Plot all the daily plot in the website or on screen if products is None 
    """
    if not dld.status:
        return 
    fig = [dld.num*10]
    def inc(fig):
        fig[0] +=1
        return fig[0]

        
    num = dld.num    
    ### Horizontal and Vertical Corrections 
    for k  in ["H", "V"]:
        if products:
            pfile = todaydir.fpath("DL%d_%s_CORR.png"%(num,k))
            with pfile.open("wb") as g:  
                dld.rail.plot.deformations(k, figure=inc(fig), fclear=True, save=g)
            log.notice("Correction figure copied to %s"%pfile)
            products.add(num, "%s corection plot"%k, "d", "img", pfile.filename)
        else:
            dld.rail.plot.deformations(k, figure=inc(fig), fclear=True, show=True)    
    

    for k in ["theta", "psi", "phi"]:
        if products:                            
            pfile = todaydir.fpath("DL%d_%s_theta.png"%(num,k))
            with pfile.open("wb") as g:  
                dld.carriage.plot.wobble_fit(k, figure=inc(fig), fclear=True, save=g)
            #archive.copy(pfile)
            log.notice("Wobble figure copied to %s"%pfile)        
            products.add( num, "Wobble %s plot"%k, "d", "img", pfile.filename)
        else:
            dld.carriage.plot.wobble_fit(k, figure=inc(fig), fclear=True, show=True)
    
    if products:                            
        pfile = todaydir.fpath("DL%d_flatness.png"%(num))
        with pfile.open("wb") as g:  
            dld.carriage.plot.flatness(figure=inc(fig), fclear=True, save=g)
        #archive.copy(pfile)
        log.notice("Flatness figure copied to %s"%pfile)        
        products.add( num, "Flatness plot", "d", "img", pfile.filename)
    else:
        dld.carriage.plot.flatness(figure=inc(fig), fclear=True, show=True)




def run_daily_plot_hysteresis(hyst, todaydir, products, log):
    """ 
    Plot all the daily plot related to hysteresis in the website or on screen if products is None 
    """  
    if not hyst.status:
        return 
            
    num = hyst.direct.num
    fig = [num*10+100]
    def inc(fig):
        fig[0] +=1
        return fig[0]
        
    for k in ["theta", "psi"]:
        ## histeresis plot
        if products:
            pfile = todaydir.fpath("DL%d_%s_histeresis.png"%(num,k))
            with pfile.open("wb") as g:
                hyst.plot.histeresis(k, fclear=True, figure=inc(fig), save=g)
                
                log.notice("Histeresis %s plot updated to  %s"%(k,pfile))
                products.add(num, "Hysteresis %s plot"%k, "d",  "img", pfile.filename)
        else:
            hyst.plot.histeresis(k, fclear=True, figure=inc(fig), show=True)


def run_monitoring_plot_dl(num, monitoringdir, products, log):

    fig = [num*10+200]
    def inc(fig):
        fig[0] +=1
        return fig[0]
            
    
    ## get the wobble monitoring file 
    wblfile = io.WobbleLog(monitoringdir.fpath("DL%d_wobble.txt"%num)).build()
    # get the data
    wdata = wblfile.read_data()
    for k in ["theta", "psi"]:
        if products:
            web_file = monitoringdir.fpath("DL%d_%s_monitoring.png"%(num,k)) 
            with web_file.open("wb") as g:
                plots.plot_history(wdata, k, num, "Wobble", fclear=True, figure=inc(fig), save=g)
                log.notice("Wobble monitoring updated to  %s"%web_file)       
                products.add( num, "Wobble %s monitoring plot"%k, "m", "img", web_file.filename)
        else:
            plots.plot_history(wdata, k, num, "Wobble", fclear=True, figure=inc(fig), show=True)        

        if products:
            web_file = monitoringdir.fpath("DL%d_%s_temp.png"%(num,k)) 
            with web_file.open("wb") as g:
                plots.plot_temperature(wdata, k, num, "Wobble", fclear=True,figure=inc(fig), save=g)
                log.notice("Wobble monitoring updated to  %s"%web_file)       
                products.add( num, "Wobble %s temp plot"%k, "m", "img", web_file.filename)
        else:
            plots.plot_temperature(wdata, k, num, "Wobble", fclear=True,figure=inc(fig), show=True)
    
    ## get the correctionfile monitoring (number of correction per day)
    correctionfile = io.CorrectionLog(monitoringdir.fpath("DL%d_correction.txt"%num)).build()

    ## read the correction data
    wdata = correctionfile.read_data()        
    for k in ["Nv", "Nh"]:  
        if products:  
                web_file = monitoringdir.fpath("DL%d_correction_%s_monitoring.png"%(num,k)) 
                with web_file.open("wb") as g:
                    plots.plot_history(wdata, k, num, "Corrections", unit="#", fclear=True,figure=inc(fig), save=g)
                    log.notice("Coorection monitoring updated to  %s"%web_file)       
                    products.add( num, "Correction %s monitoring plot"%k, "m", "img", web_file.filename)    
        else:
            plots.plot_history(wdata, k, num, "Corrections", unit="#", fclear=True,figure=inc(fig), show=True)

        if products:                
                ##  number of correction vs temperature 
                web_file = monitoringdir.fpath("DL%d_correction_%s_temp.png"%(num,k)) 
                with web_file.open("wb") as g:
                    plots.plot_temperature(wdata, k, num, "Corrections", unit="#", fclear=True,figure=inc(fig), save=g)
                    log.notice("Correction vs temp updated to  %s"%web_file)       
                    products.add( num, "Correction %s temp plot"%k, "m", "img", web_file.filename)                    
        else:
            plots.plot_temperature(wdata, k, num, "Corrections", unit="#", fclear=True,figure=inc(fig), show=True)        

def run_monitoring_plot_hysteresis(num, monitoringdir, products, log):

    fig = [num*10+300]
    def inc(fig):
        fig[0] +=1
        return fig[0]
    

    
    ## get the Hysteresis  monitoring file 
    hystfile = io.HysteresisLog(monitoringdir.fpath("DL%d_hysteresis.txt"%num)).build()
    # get the data
    wdata = hystfile.read_data()
    for k in ["hysteresis"]: 
        if products:
            web_file = monitoringdir.fpath("DL%d_%s_monitoring.png"%(num,k)) 
            with web_file.open("wb") as g:
                plots.plot_history(wdata, k, num, "Hysteresis", unit="arcsec", fclear=True,figure=inc(fig), save=g)
                log.notice("Hysteresis monitoring updated to  %s"%web_file)       
                products.add(num, " %s monitoring plot"%k, "m", "img", web_file.filename)
        else:
            plots.plot_history(wdata, k, num, "Hysteresis", unit="arcsec", fclear=True,figure=inc(fig), show=True)            
            
        if products:    
            web_file = monitoringdir.fpath("DL%d_%s_temp.png"%(num,k)) 
            with web_file.open("wb") as g:
                plots.plot_temperature(wdata, k, num, "Hysteresis", unit="arcsec", fclear=True,figure=inc(fig), save=g)
                log.notice("Hysteresis vs temp updated to  %s"%web_file)       
                products.add(num, "%s temp plot"%k, "m", "img", web_file.filename)
        else:
            plots.plot_temperature(wdata, k, num, "Hysteresis", unit="arcsec", fclear=True, figure=inc(fig), show=True)            



def run_computing(date=None, dls=range(1,7)):

    ## get a list of direct files
    fdirect  = dldir.deliriumfiles(date, "direct")
    ## get a list of reverse file if any 
    freverse = dldir.deliriumfiles(date, "reverse")

    ## status of computation for each dls True for succsess
    status = {} 
    for num in dls:
        context = ("DELIRIUM", "DL%d"%num)
        clear_buffer(context[1])
        log = Log(context=context)
        fls = fdirect.get(num, [])
        if len(fls)>1:
            ## take the last file (they are laready sorted from older to most recent)
            file = fls[-1]
            log.warning("%d files found taking the most recent one"%(len(fls)))
        elif not fls:
            log.error("No delirium file found for dl %d"%num)
            status[num] = False            
            continue
        else:
            file = fls[-1]  

       
        

                
                




def write_failure(context, origfile, files):
    for file in files:
        with file.open("w") as g:
            g.write("!!!! Corection could not be computed. Got the following errors : \n ")
            g.write("script executed at %s\n UT"%(time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime())))
            g.write("on file %s\n"%origfile) 
            g.write("".join(get_buffer(context, "ERROR")))

