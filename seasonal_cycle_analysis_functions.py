import numpy as np
import xarray as xr
from matplotlib.colors import Normalize

def gavg(idata):
	"""calculate global average
	e.g., x1=gavg(d1['t2m'])"""	
	wgt1=np.cos(np.deg2rad(idata.latitude))*(idata*0+1)
	ga=(wgt1*idata).sum(dim=['latitude','longitude'])/wgt1.sum(dim=['latitude','longitude'])
	return ga
	
def seasonal_amp1D(sig):
    return sig.max()-sig.min()

def seasonal_amp(var,ax): 
    return np.apply_along_axis(seasonal_amp1D,axis=ax,arr=var)


#--- a function that computes the saturation mixing ratio
#--- given pressure and temperature
def sat_mix_mk(p,t):
    esi = np.zeros_like(t)
    ilarge = np.where(t>273.5)
    ismall = np.where(t<=273.5)
    esi[ilarge] = np.exp(54.842763-6763.22/t[ilarge]-4.210*np.log(t[ilarge]) +\
        0.000367*t[ilarge]+np.tanh(0.0415*(t[ilarge]-218.8))*(53.878-1331.22/t[ilarge] - \
        9.44523*np.log(t[ilarge])+0.014025*t[ilarge]))/ 100. # hPa
    
    esi[ismall] = np.exp(9.550426-5723.265/t[ismall]+3.53068*np.log(t[ismall])-0.00728332*t[ismall])/100.0 # hPa
    
    h2o=esi/(p-esi) # the mixing ratio
    ilarge2 = np.where(esi*100>p)
    h2o[ilarge2] = esi[ilarge2]/p[ilarge2] # if esi is too large, use specific humidity
    h2o[h2o>=1] = np.nan 
    return h2o*1e6

#--- a function that helps to setup the colorbars for the figrues
def colorbar_setting(fig,im,ax,vmin,vmax,vnum,clabel,cbaxes):
    cbaxes = fig.add_axes(cbaxes)   
    c = fig.colorbar(im,ax=ax,cax=cbaxes,extend='both') 
    cticks = np.linspace(vmin,vmax,vnum)
    clabels = np.array([np.str(round(i,2)) for i in cticks])
    c.set_ticks(cticks)
    c.set_ticklabels(clabels)
    c.set_label(clabel,fontsize=10) 
    c.ax.tick_params(labelsize=8)
    

class PiecewiseNorm(Normalize):
    def __init__(self, levels, clip=False):
        # the input levels
        self._levels = np.sort(levels)
        # corresponding normalized values between 0 and 1
        self._normed = np.linspace(0, 1, len(levels))
        Normalize.__init__(self, None, None, clip)

    def __call__(self, value, clip=None):
        # linearly interpolate to get the normalized value
        return np.ma.masked_array(np.interp(value, self._levels, self._normed))

