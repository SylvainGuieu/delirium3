#import path as _p
from path import *
import time
import os
import numpy as np
from .log import Log
import datetime
import json 
log = Log(context=("DELIRIUM, I/O"))

DIRECT = "direct"
REVERSE = "reverse"


try:
    unicode
except NameError:
    basestring = (str, bytes)


def _path2dl(path):
    """ return the dl number of a given path or None in case of failure """
    directory, filename = os.path.split(path)
    pref, rest = filename.split("_",1)
    if pref[0:2] != "DL":
        return None
        #raise ValueError("cannot read dl number on file %s"%path)
    try:
        dl = int(pref[2:])
    except:
        return None
        #raise ValueError("cannot read dl number on file %s"%path)
    return dl

class DeliriumDirectory(dpath):

    def _parse_time(self,ftime=None):
        if ftime is None:
            ftime = time.localtime()
        elif isinstance(ftime, basestring):
            ftime = time.strptime(ftime, "%Y%m%d")
        elif not isinstance(ftime, time.struct_time):
            raise ValueError("time should be string or struct_time got '%s'"%ftime) 

        return ftime
    def _parse_time(self,ftime=None):
        if ftime is None:
            ftime = time.localtime()
        elif isinstance(ftime, basestring):
            ftime = time.strptime(ftime, "%Y-%m-%d")
        elif not isinstance(ftime, time.struct_time):
            raise ValueError("time should be string or struct_time got '%s'"%ftime) 

        return ftime

    def _time2direct_file_glob(self, ftime):
        """ from a time structure get the glob file """
        return "DL*_FOGALE_%s???.dat"%(time.strftime("%Y%m%d", ftime))

    def _time2reverse_file_glob(self, ftime):
        """ from a time structure get the glob file """
        return "DL*_FOGALE_%s*REV.dat"%(time.strftime("%Y%m%d", ftime))



    def deliriumfiles(self, ftime, ftype=DIRECT, inside=None):
        """ get DELIRIUM files according to a given date

        Parameters
        ----------
        ftime : struc_time, None or string
            if string must be of format "yyyymmdd"
            if None time.localtime() is taken

        ftype : string
            "direct" or "reverse", the type of file      
        """
        ftime = self._parse_time(ftime)         
        
        if ftype==DIRECT:
            glob = self._time2direct_file_glob(ftime)
        elif ftype==REVERSE:
            glob = self._time2reverse_file_glob(ftime)  
        else:
            raise TypeError("unrecognized file type '%s'"%ftype)
        dtime =  time.strftime( "%Y-%m", ftime)   
        glob =  os.path.join( dtime, glob)
        if inside:
            files= self.get(glob, os.path.join(inside, dtime), child=DeliriumPath )
            #files = self.get(glob, os.path.join(inside, dtime), child=DeliriumPath)            
        else:
            files = self.ls(glob, child=lambda x:DeliriumPath(self, x))
                                
        output = {}
        for file in files:
            dl = _path2dl(file)
            if dl in output:
                output[dl].append(file)
            else:
                output[dl] = [file]        
        # sort by modification date
        for fls in output.values():
            if len(fls)<2: continue
            try:
               fls.sort(key=getmtime)
            except:
               log.warning("Cannot sort files by modification date") 
               break

        return output   

class DeliriumPath(fpath):
    header_lookup = {
        "tunnel_temp" : ("Tunnel Temperature=", float),
        "date" : ("Date=", lambda s: time.strptime(s, "%Y-%m-%dT%H:%M:%S"))
    }
    header_txt = []
    header = None 

    def __repr__(self):
        return "delirium'%s'"%self

    def get_header(self):
        if not self.header:
            self.load_header()
        return self.header

    def load_header(self):
        self.header_txt = []
        self.header = dict((k,-999.99) for k in self.header_lookup.keys())

        f = open(self)
        while True:
            offset = f.tell()
            line = f.readline().strip()
            if not line:
                break
            if not line[0:1]=="%":
                # put the file backward
                f.seek(offset)
                break

            self.header_txt.append(line)
            for var, (search, tpe) in  self.header_lookup.items():
                where = line.find(search)
                if where>-1:
                    line =  line[where+len(search):]
                    sval = line.split(" ")[0]
                    try:
                        val = tpe( sval  )
                    except (TypeError, ValueError):
                        log.notice("Warning cannot read '%s' header parameter"%var)    
                    self.header[var] = val                            


class MonitorDirectory(dpath):
    _v_ =  "2"
    @property    
    def hdir(self):
        return HysteresisDiretory(self)
        #return HysteresisDiretory(self.dpath("DL_hysteresis%s"%self._v_))
    @property    
    def wdir(self):
        return WobbleDiretory(self)
        #return WobbleDiretory(self.dpath("DL_wobble%s"%self._v_))
    @property    
    def cdir(self):
        return CorrectionDiretory(self)
        #return WobbleDiretory(self.dpath("DL_wobble%s"%self._v_))
    
        
class MonitoringLog(fpath):
    ## make the create verbose
    def create(self,header=None):        
        fpath.create(self, header)
        log.notice("File %s created"%self)
            
def read_date(date):
    return datetime.date(*[int(v) for v in date.split("-")])

        
class HysteresisDiretory(dpath):    
    def hysteresisfiles(self):
        return self.als("DL*hyst.txt", HysteresisLog)
    def hysteresisfile(self, dlnum):
        return self.als("DL%d*hyst.txt"%dlnum, HysteresisLog)[0]
           
    def hysteresisfigure(self, dlnum):
        return self.als("DL%d*hyst.jpg"%dlnum)[0]   


class HysteresisLog(MonitoringLog): 
    header =  "%% Date  TunnelTemp(deg) Hysteresis(arcsec)\n"

    def __repr__(self):
        return "histeresys'%s'"%self

    def add_from_hysteresis(self, hyst):
        """ Take a delayline object and add an entry in the file """

        date = hyst.direct.delirium.get_date2()
        date = datetime.date(*[int(v) for v in date.split("-")])        
        temp = hyst.direct.delirium.get_tunnel_temp()
        data = self.read_data()

        if len(np.where( (data['date']==date) * (data['temp']==temp) )[0]):
            log.notice("Entry '%s' for hysteresis already existing"%date)
        else: 
            hysteresis =  hyst.psi_diff_m / (np.pi/180/3600.)
            self.add_entry(date, temp, hysteresis)           

    def add_entry(self,  date, temp, hysteresis):
        with open(self, "a") as f:
            f.write( "%s %4.2f %9.2f\n"%(date, temp, hysteresis))
    
    def read_data(self):
        lst = []
        with self.open('r') as f:
            for line in f:
                line = line.strip()
                if line[0]!="%":
                    try:                    
                        tmp =  tuple( t(l.strip()) for t,l in zip([read_date,float,float], [v for v in line.split(" ") if v]) )                    
                        lst.append(tmp)
                    except ValueError:                    
                        log.warning("'%s' : line '%s' corrupted"%(self,line))                    
        return np.array(lst, dtype=[("date", "O"), ("temp", "f8"), ("hysteresis", "f8")])

    def data_of(self, date):
        d = self.read_data()
        return d[d['data']==date]    

class WobbleDiretory(dpath):    
    def wobblefiles(self):
        return self.als("DL*wobble.txt", WobbleLog)
    def wobblefile(self, dlnum):
        return self.als("DL%d_wobble.txt"%dlnum, WobbleLog)[0]
    def wobblefigure(self, dlnum):
        return self.als("DL%d*_wobble.jpg"%dlnum)[0]   



class WobbleLog(MonitoringLog):
    
    header = "%% Date  TunnelTemp(deg) WoobleTheta(arcsec) WooblePsi(arcsec)\n"   
    def __repr__(self):
        return "wobble'%s'"%self


    def add_from_dl(self, dl):
        """ Take a delayline object and add an entry in the file """

        date = dl.delirium.get_date()
        #date = date[0:4] + "-" + date[4:6] + "-" + date[6:8]
        date = datetime.date(int(date[0:4]), int(date[4:6]), int(date[6:8]))
        temp = dl.delirium.get_tunnel_temp()
        data = self.read_data()

        _, e = dl.carriage.get("theta", extras=True, filterWobble=True,arcsec=True)
        theta = e["wobble_amplitude"]
        dl.carriage.log.notice("Wobble Amplitude theta %4.2f arcsec"%theta)
        _, e = dl.carriage.get("psi", extras=True, filterWobble=True,arcsec=True)
        psi = e["wobble_amplitude"]
        dl.carriage.log.notice("Wobble Amplitude psi %4.2f arcsec"%psi)

        if len(np.where( (data['date']==date) * (data['temp']==temp) )[0]):
            log.notice("Entry '%s' for wobble already existing"%date)
        else:            
             
            ## add the entry of today 
            self.add_entry(date, temp, theta, psi)

    def add_entry(self,  date, temp, theta, psi):
        with open(self, "a") as f:
            f.write( "%s %4.2f %8.2f %8.2f\n"%(date.isoformat(), temp, theta, psi))
    
    def data_of(self, date):
        d = self.read_data()
        return d[d['data']==date]

    def read_data(self):
        lst = []
        with self.open('r') as f:
            for line in f:
                line = line.strip()
                if line[0]!="%":
                    try:                    
                        tmp =  tuple( t(l.strip()) for t,l in zip([read_date,float,float,float], [v for v in line.split(" ") if v]) )                    
                        lst.append(tmp)
                    except ValueError:                    
                        log.warning("'%s' : line '%s' corrupted"%(self,line))                    
        return np.array(lst, dtype=[("date", "O"), ("temp", "f8"), ("theta", "f8"),("psi", "f8") ])

class CorrectionDiretory(dpath):    
    def correctionfiles(self):
        return self.als("DL*correction.txt", CorrectionLog)

    def correctionfile(self, dlnum):
        return self.als("DL%d*correction.txt"%dlnum, CorrectionLog)[0]
           
    def correctionfigure(self, dlnum):
        return self.als("DL%d*correction.jpg"%dlnum)[0]   


class CorrectionLog(MonitoringLog):
    header = "%% Date  TunnelTemp(deg) Vcorrection(number) Hcorrection"

    def __repr__(self):
        return "corrections'%s'"%self

    def add_from_dl(self, dl):
        """ Take a delayline object and add an entry in the file """

        date = dl.delirium.get_date()
        #date = date[0:4] + "-" + date[4:6] + "-" + date[6:8]
        date = datetime.date(int(date[0:4]), int(date[4:6]), int(date[6:8]))
        temp = dl.delirium.get_tunnel_temp()
        data = self.read_data()

        cor = dl.supports.get_corrections()
        
        Nv = len(cor[abs(cor['V'])>0])
        Nh = len(cor[abs(cor['H'])>0])        

        if len(np.where( (data['date']==date) * (data['temp']==temp) )[0]):
            log.notice("Entry '%s' for number of coorection already existing"%date)
        else:            
             
            ## add the entry of today 
            self.add_entry(date, temp, Nv, Nh)        

    def add_entry(self,  date, temp, Nv, Nh):
        with open(self, "a") as f:
            f.write( "%s %4.2f %d %d\n"%(date.isoformat(), temp, Nv, Nh))
    
    def data_of(self, date):
        d = self.read_data()
        return d[d['data']==date]

    def read_data(self):
        lst = []
        with self.open('r') as f:
            for line in f:
                line = line.strip()
                if line[0]!="%":
                    try:                    
                        tmp =  tuple( t(l.strip()) for t,l in zip([read_date,float,int,int], [v for v in line.split(" ") if v]) )                    
                        lst.append(tmp)
                    except ValueError:                    
                        log.warning("'%s' : line '%s' corrupted"%(self,line))                    
        return np.array(lst, dtype=[("date", "O"), ("temp", "f8"), ("Nv", "f8"),("Nh", "f8") ])
                        




class Product(fpath):
    """ A class to handle the writing of delirium products in a javascript
    web page.
    """

    def __new__(cl, file, date=None, script_date=None, dirs={}):
        self = fpath.__new__(cl, file)
        self.date = date
        if script_date is None:
            script_date = datetime.datetime.utcnow().isoformat()
        self.script_date = script_date    
        self.products = {}
        self.dirs = dict(dirs)

        self.load()
        return self

    def load(self):        
        if self.exists():
            with fpath(self).replace_ext(".json").open('r') as f:
                try:
                    self.products = json.load(f)
                except ValueError:
                    pass
                else:
                    # Quick patch json keys are always string not integer
                    # change the dlnum string to integers
                    for sdl in list(self.products['data'].keys()):
                        self.products['data'][int(sdl)] = self.products['data'][sdl]
                        del self.products['data'][sdl]

    def add(self, dlnum, name, kind, ptype, *args):
        """ add a product to the list 
    
        Parameters
        ----------
        dlnum : int
            delay line number
        name : string
            product name
        kind : string
            'd' for daily products, 'm' for monitoring
        ptype : string
            'img', 'text', 'file'            
        *args : depend of ptype
        """        
        products = self.products
        
        kkind = "daily_product_names" if kind is 'd' else "monitoring_product_names"

        args = list(args)        
        if ptype in ("img","txtfile", "htmlfile"):
            args[0] = os.path.join(self.dirs.get(kind, ""), args[0])
            
        products.setdefault("data" , {})
        names = products.setdefault(kkind, [])


        products['data'].setdefault(dlnum,{}).update( {name:[ptype]+list(args)})
        if not name in names:
            names.append(name)
    

    def product2json(self):
        """ transform the product to a string json file """
        products = self.products
        text = ""
        #text += """   "daily_product_names": [\n    %s\n],\n"""%(",\n    ".join('"%s"'%n for n in products['daily_product_names']))
        text += """   "daily_product_names": [%s],\n"""%(", ".join('"%s"'%n for n in products.get('daily_product_names', [])))
        #text += """   "monitoring_product_names": [\n    %s\n],\n"""%(",\n    ".join('"%s"'%n for n in products['monitoring_product_names']))
        text += """   "monitoring_product_names": [%s],\n"""%(", ".join('"%s"'%n for n in products.get('monitoring_product_names',[])))

        text += """   "dlsnum": [%s],\n"""%(",".join("%d"%n for n in products.get('data', {})))

        ind = " "*4
        blocks = []
        for dlnum,subproducts in products.get('data', {}).items():
            block = ind+""""%d" : {\n"""%dlnum

            subblocks = []
            for name, args in subproducts.items():
                ptype = args[0]
                args = args[1:]


                subblock = """"%s" : ["%s","""%(name, ptype)

                subblock += ",".join('"%s"'%a for a in args)
                subblock += "]"
                subblocks.append(subblock)
            block += (ind*2)+((",\n"+ind*2).join(subblocks))
            

            block += "\n"+ind+"}" 
            blocks.append(block)   
        
        text += """   "data" : {\n%s\n}"""%(",\n".join(blocks))
        return text              

    def product2js(self):
        """ transform the product to a string javascript executable """
        products = self.products
        text = "daily_product_names = [\n    %s\n];\n"%(",\n    ".join("'%s'"%n for n in products['daily_product_names']))
        text += "monitoring_product_names = [\n    %s\n];\n"%(",\n    ".join("'%s'"%n for n in products['monitoring_product_names']))

        text += "dlsnum = [%s];\n"%(",".join("%d"%n for n in products['data']))

        ind = " "*4
        blocks = []
        for dlnum,subproducts in products['data'].items():
            block = ind+"%d : {\n"%dlnum

            subblocks = []
            for name, args in subproducts.items():
                ptype = args[0]
                args = args[1:]
                subblock = "'%s' : ['%s',"%(name, ptype)

                subblock += ",".join("'%s'"%a for a in args)
                subblock += "]"
                subblocks.append(subblock)
            block += (ind*2)+((",\n"+ind*2).join(subblocks))
            

            block += "\n"+ind+"}" 
            blocks.append(block)   
        
        text += "data = {\n%s\n};"%(",\n".join(blocks))
        return text            

    def flushjs(self):
        with self.open('w') as g:
            g.write("date = '%s';\n"%self.date)
            g.write("script_date = '%s';\n"%self.script_date)
            g.write(self.product2js())

    def flushjson(self):
        with fpath(self).replace_ext(".json").open('w') as g:
            g.write("""{\n   "date": "%s",\n"""%self.date)
            g.write("""   "script_date": "%s",\n"""%self.script_date)
            g.write(self.product2json())
            g.write("\n}")        
                    
