
import sys
import pymisca.header as pyheader
pyheader.set__numpy__thread(locals().get('NCORE',None))
# if "matplotlib" not in sys.modules.keys():
pyheader.mpl__setBackend('agg')
# import matplotlib as mpl; mpl.use('Agg');
import synotil.dio as sdio;reload(sdio)
import synotil as synotil; reload(synotil)
import synotil.qcplots as pkg; reload(pkg)
import synotil.util as sutil;reload(sutil)
import synotil.PanelPlot as spanel;reload(spanel)
import synotil.CountMatrix as scount;reload(scount)
import synotil.jobs as sjob
import synotil.norm as snorm

import pymisca.shell as pysh
import pymisca.util as pyutil; reload(pyutil)
import pymisca.vis_util as pyvis; reload(pyvis)
import cPickle as pk
import funcy


np = pyutil.np; pd = pyutil.pd
plt = pyutil.plt; 
if pyutil.hasIPD:
    get_ipython().magic('matplotlib inline')