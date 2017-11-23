from __future__ import print_function
import time
import sys

ERROR   = 1
WARNING = 2
NOTICE  = 4 # info and notice are aliases
INFO = 4
DATA = 8 

# Buffer size
BUFFER_SIZE = 100 

# default verbose type 
verbose_type = ERROR+WARNING+NOTICE+DATA
# default level of verbose 
verbose_level = 3

##
# convert the MSG code to string 
msgtype2string_lookup = {ERROR:"ERROR", WARNING:"WARNING", NOTICE:"NOTICE", DATA:"DATA"}

counters = dict( (m,0) for m in msgtype2string_lookup)
buffers  = {}

##
# convert the MSG code to string when is printed to stoud
try:
    import colorama as col
except:
    colorized_msgtype2string_lookup = msgtype2string_lookup
else:
    colorized_msgtype2string_lookup = {
                 ERROR  :col.Fore.RED+"ERROR"+col.Fore.RESET,
                 WARNING:col.Fore.MAGENTA+"WARNING"+col.Fore.RESET,
                 NOTICE :col.Fore.BLUE+"NOTICE"+col.Fore.RESET,
                 DATA   :col.Fore.BLUE+"DATA"+col.Fore.RESET
                }

try:
    unicode
except NameError: #python 3
    basestring = (str, bytes)

def toggle_color(flag):
    global stdout_msgtype2string_lookup, msgtype2string_lookup , colorized_msgtype2string_lookup
    if flag:
        stdout_msgtype2string_lookup = colorized_msgtype2string_lookup
    else:
        stdout_msgtype2string_lookup = msgtype2string_lookup        

toggle_color(True)

def clear_counter(tpe):
    global counters
    counters[tpe] = 0

def clear_counters():
    global counters
    for k in counters:
        counters[k] = 0

def counter(tpe):
    global counters
    counters[tpe] += 1




def get_buffer_keys(context):
    global buffers
    if context is None:
        return []
    if isinstance(context, basestring):
        context = (context,)
    keys = list(buffers.keys())
        
    for c in context:
        for k in list(keys):
            if not c in k:
                keys.remove(k)
    return keys    

def get_buffer(context, tpe):

    keys = get_buffer_keys(context)
    out = []
    for k in keys:
        if tpe in buffers[k]:
            out.extend(buffers[k][tpe])
    return out

def clear_buffer(context, tpe=None):
    keys = get_buffer_keys(context)
    if tpe is None:
        for k in keys:
            buffers.pop(k,None)
    else:
        for k in keys:
            buffers[k].pop(tpe, None)        

def add_to_buffer(context, tpe, msg):
    global buffers
    if context is None:
        return []

    if isinstance(context, basestring):
        context = (context,)
    if not context in buffers:    
        buffers[context] = {}
    if not tpe in buffers[context]:
        buffers[context][tpe] = []
    b =  buffers[context][tpe]
    b.append(msg)                           
    for i in range(len(b)-BUFFER_SIZE):
        b.pop(0)

def clear_buffers():
    global buffers
    for k in buffers:
        buffers[k] = []

def buffer(tpe, msg):
    b = buffers[tpe]
    b.append(msg)
    for i in range(len(b)-BUFFER_SIZE):
        b.pop(0)


class Log(object):
    """A Log set of function for Gravity """
    default_format = """{context} {msgtype} {date}: {msg}\n"""         
    formatlookup  = {DATA:"{msg}"}
    date_format = "%Y-%m-%dT%H:%M:%S"

    msgtypes  = verbose_type
    maxlevel  = verbose_level
    context = None
    def __init__(self, outputs=None, msgtypes=None, maxlevel=None, context=None):
        """ create a log writing object.


        A log object is defined by 
         -outputs (can be stdout a file, a list, ...) 
         -msgtypes allowed message type for each outputs: ERROR, WARNING, NOTICE or DATA
         -maxlevel a maximum level for each output 
         -context : string or function that set the log context
         -default_format :  a default format for all msgtype 
         -formatlookup : a dictionary of format for each types.
                the format accept the following keys:
                    context, msgtype, msglevel, date, count, msg
                the date format can be defined with the data_format attribute    
                              
        Parameters
        ----------
        outputs : list, optional
            a list of : string -> path to a file opened in "w" mode 
                        file   -> f.write is used to output the message
                        callable object -> with a signnature func(msg) where msg is a string
                        tuple -> if a tuple must be (f, maxlevel) or (f, maxlevel, msgtypes)
            
                                where f can be a file or a string.
                                This way one can set a different maxlevel and msgtype 
                                for each output
                                maxlevel and msgtype can be omited or None to set to their 
                                default value                                    
            e.g. :
                log = Log(outputs=[(sys.stdout, NOTICE+DATA), (sys.stderr, ERROR)])                               
            If no output is given the default is stdout.
            
        msgtypes : int(binary) or string optional 
            default msgtypes for each outputs if not defined.
            each msgtypes bits turn on/off the allowed msgtype.
            One can use the sum combination of the defined constants ERROR, WARNING, NOTICE and DATA
            
            The bytes are as follow:
                #byte  int  msgtype 
                1      1    ERROR 
                2      2    WARNING
                3      4    NOTICE
                4      8    DATA  

            e.g. :   msgtypes = ERROR+WARNING  -> will print only the error and warning messages   

            Also msgtypes can be astring containing a combination of the caracters  'E','W','N' or 'D'
             msgtypes = ERROR+WARNING+NOTICE equivalent to msgtypes = "EWN" 
    
        maxlevel : int, optional
            The default maxlevel for each output if not defined.    
            If not given maxlevel=1 

        context : string or func
            a string or function that return a string defining the log context.
            e.g.: 
                class Data:
                    def __init__(self, file):
                        # do something with file
                        self.log = Log( context=lambda : "DATA %s"%(self.file) )
        

        """

        if maxlevel is not None:
            self.maxlevel = maxlevel
        if msgtypes is not None:
            if isinstance(msgtypes, basestring):
                msgtypes = sum( (s in msgtypes)*bit for s,bit in [('E',ERROR),('W',WARNING),('N',NOTICE),('D',DATA)] )
            self.msgtypes   = msgtypes

        if isinstance(context, basestring):
            context = (context,)
        if context is not None:    
            self.context = tuple(context)
                                
        #    self.logcontext = lambda :context
        #elif context:
        #    if not hasattr(context, "__call__"):
        #        raise ValueError("Context must be a string or a function got a %s"%type(context))    
        #    self.logcontext = logcontext  
                
        self.outputs = []
        self._files_ = []
                
        if not outputs:
            self.add_output(sys.stdout.write, None, None)
        else:
            for output in outputs:
                output = output if isinstance(output, tuple) else (output,)       
                self.add_output(*output)

        # make a copy of the class default        
        self.formatlookup = dict(self.formatlookup)        
        self.count = 0
        self.noisy()


    def logcontext(self):
        if self.context:
            return "[%s]"%(" ".join(self.context))
        return ""        

    def add_output(self, wf,  msgtypes=None, maxlevel=None):

        msgtype2string = lambda tpe: msgtype2string_lookup.get(tpe,"") 
        if isinstance(wf, basestring):
            wf = open(wf, "w").write

        #if isinstance(wf, file):
        #    if (wf is sys.stdout) or (wf is sys.stderr):
        #        msgtype2string = lambda tpe: stdout_msgtype2string_lookup.get(tpe,"")         
        #
        #    self._files_.append(wf)
        #    wf = wf.write

        #msgtypes = self.msgtypes if msgtypes is None else msgtypes
        #maxlevel = self.maxlevel if maxlevel is None else maxlevel        
        self.outputs.append( (wf,  msgtypes, maxlevel, msgtype2string) )

    def close(self):
        """ Attempt to close every output files 

        The log stay usable for all other output that are not files
        """
        for f in self._files_:
            if (f is sys.stdout) or (f is sys.stderr): continue
            self.outputs = [(wf, maxlevel, msgtypes, t2s) for wf, maxlevel, msgtypes,t2s in self.outputs if wf!=f.write]
            f.close()

    def quiet(self):
        """ put the log quiet 

        use log.noisy() to put it back
        """
        self._disabled = True

    def noisy(self):
        """ send back the log to normal behavior """    
        self._disabled = False

    def log(self, msg, level=1, msgtype=NOTICE):
        if self._disabled:
            return 
        for wf,  msgtypes , maxlevel, msgtype2string in self.outputs:
            msgtypes = self.msgtypes if msgtypes is None else msgtypes
            maxlevel = self.maxlevel if maxlevel is None else maxlevel

            ## do not do anything if level is too hight            
            if level>maxlevel: continue
            ## if msgtype is not in msgtypes, return 
            if not msgtypes & msgtype: continue

            fmt = self.formatlookup.get(msgtype, self.default_format)
            formated_msg = fmt.format(
                                context=self.logcontext(),
                                msgtype=msgtype2string(msgtype),
                                msglevel=level,                                        
                                date=time.strftime(self.date_format, time.gmtime()),
                                msg =msg, 
                                count=self.count
                                )
            wf(formated_msg)            
            add_to_buffer(self.context, msgtype, formated_msg)
            self.count += 1

    def data(self, msg, level=1):
        self.log(msg, level, DATA)

    def info(self, msg, level=1):
        self.log(msg, level, NOTICE)

    def notice(self, msg, level=1):
        self.log(msg, level, NOTICE)
            
    def warning(self, msg, level=1):
        self.log(msg, level, WARNING)

    def error(self, msg, level=1):
        self.log(msg, level, ERROR)	       

## open a new log         
log = Log()    
