try:
    from astropy.io import fits
except:
    import pyfits as fits

import datetime
import numpy as np

def create_coldef():
    coldef = fits.ColDefs([
            fits.Column(name='year', format='I'),
            fits.Column(name='month', format='I'),
            fits.Column(name='day', format='I'),
            fits.Column(name='temp', format='E', unit="celcius"),
            fits.Column(name='wobble_psi', format='E', unit="arcsec"),
            fits.Column(name='wobble_theta', format='E', unit="arcsec"),
            fits.Column(name='wobble_phi', format='E', unit="arcsec"),
            fits.Column(name='vertical_ncor', format='I'),
            fits.Column(name='horizontal_ncor', format='I'),
            fits.Column(name='support_numbers', format='74I'),
            fits.Column(name='vertical_corrections', format='74E', unit="micron"),
            fits.Column(name='horizontal_corrections', format='74E', unit="micron")
        ]
    )
    return coldef



def monitoring_data():
    return {k:[] for k in ['year', 'month','day','temp', 'wobble_psi', 'wobble_theta', 'wobble_phi',
                           'vertical_ncor', 'horizontal_ncor', 
                           'support_numbers',
                           'vertical_corrections', 'horizontal_corrections']}

def build_dls_data(dls_list):
    data = monitoring_data()
    for dl in dls_list:
        update_dl_data(dl, data)

def update_dl_data(dl, data):
                             
    dt = dl.get_date()
    year, month, day = int(dt[0:4]), int(dt[4:6]), int(dt[6:8])

    # dates = list(zip(data['year'],data['month'],data['day']))

    # if (year, month, day) in dates:
    #     print(dates, "already exists")
    #     return 



    data['year'].append(year)
    data['month'].append(month)
    data['day'].append(day)
            
    data['temp'].append(dl.delirium.get_tunnel_temp())

    _, e = dl.carriage.get("theta", extras=True, filterWobble=True,arcsec=True)
    theta = e["wobble_amplitude"]
    _, e = dl.carriage.get("psi", extras=True, filterWobble=True,arcsec=True)
    psi = e["wobble_amplitude"]
    _, e = dl.carriage.get("phi", extras=True, filterWobble=True,arcsec=True)
    phi = e["wobble_amplitude"]
    data['wobble_psi'].append(psi)
    data['wobble_theta'].append(theta)
    data['wobble_phi'].append(phi)

    cor = dl.supports.get_corrections()

    Nv = len(cor[abs(cor['V'])>0])
    Nh = len(cor[abs(cor['H'])>0])

    data['vertical_ncor'].append(Nv)
    data['horizontal_ncor'].append(Nh)

    data['support_numbers'].append(dl.supports.supports[0:74])
    data['vertical_corrections'].append(dl.supports.Vcorrection[0:74]*1000)
    data['horizontal_corrections'].append(dl.supports.Hcorrection[0:74]*1000)


def add_data_to_table(tbl, data):
    nrow_table = tbl.data.shape[0]
    
    

    #dates = list(zip(tbl.data['year'],tbl.data['month'],tbl.data['day']))
    # check = []
    # for year,month,day in zip(data['year'], data['month'], data['day']):
    #     if (year, month, day) in dates:
    #         print("dates %d-%02d-%02d already exists"%(year, month, day))
    #         check.append(False)
    #     else:
    #         check.append(True)

    # check = np.array(check)
    # for k,v in data.items():
    #     data[k] = np.array(v)[check]

    nrow_new = len(data['day'])
    newtable = fits.BinTableHDU.from_columns(tbl.columns, nrows=nrow_table+nrow_new)
    if not nrow_new:
        print("No data")
        return newtable

    for colname in tbl.columns.names:      
        newtable.data[colname][nrow_table:] = data[colname]

    return newtable

def add_to_table(tbl, dls_list):    
    data = build_dls_data(dls_list)
    return add_data_to_table(tbl, data)

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


def update_monitoring(fh, start, end, opener, force=False):
    dlnum = fh[0].header['DL']

    data = monitoring_data()
    dates = ["%d-%02d-%02d"%(y,m,d) for y,m,d in zip(fh[1].data['year'],fh[1].data['month'],fh[1].data['day'])]


    for date in list_dates(start, end):
        if not force and (date in dates):
            print("Entry for date %s already exists"%date)
            continue

        try:
            dl = opener(dlnum, date)
        except (ValueError, OSError, KeyError) as e:
            print(e)
            continue
        update_dl_data(dl, data)

    new = create_monitoring(dlnum)
    new[1] = add_data_to_table(fh[1], data)
        
    return new
    










def create_monitoring(dlnum):
    """ Create an empty monitoring fits file for the given DL number 

    Parameters
    ----------
    dlnum: int
        delay line number
    """
    dh = {
        "DL":(dlnum, "Delay line Number")
    }


    h = fits.Header(
            [fits.Card(k,v,c) for k,(v,c) in dh.items() ]
    )

    coldef = create_coldef()

    hdus = fits.HDUList([fits.PrimaryHDU(header=h), fits.BinTableHDU.from_columns(coldef)])
    return hdus

