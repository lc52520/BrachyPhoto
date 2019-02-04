#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pymisca.ext as pyext
NCORE = 10
execfile(pyext.base__file('headers/header__import.py'))
tks = pyutil.readBaseFile('results/0130__makeTracks-Brachy/tracks.npy').tolist()
figs = pyutil.collections.OrderedDict()


ENT_CUTOFF=2.5
STEP = 100


mdl0 = pyutil.readBaseFile('results/0129__cluster__Brachy-RNA-all/mdl.npy').tolist()
mdl = mdl0.callback.mdls[STEP][-1]
mdl.predict = pyutil.functools.partial(mdl.predict, 
                                              entropy_cutoff=ENT_CUTOFF)
clu0 = clu= mdl.predictClu(mdl0.data, index=mdl0.data.index )
cacheFile = pyutil.cache__model4data(mdl,tdf=mdl0.data,
                                     ofname='cache.npy')
# clu0 = clu= mdl.predictClu(mdl0.data, entropy_cutoff=ENT_CUTOFF,index=mdl0.data.index )
# clu.to_csv()
clu = clu.sort_values('clu')
vdf = scount.countMatrix( mdl0.data)
vdf = vdf.reindex(clu.index)
vdf.heatmap(figsize=[12,12])
plt.title('N=%d'%len(vdf))
figs['allClu'] = plt.gcf() 

vdf = scount.countMatrix( mdl0.data)
vdf = vdf.reindex(clu.query('clu!=-1').index)
vdf.heatmap(figsize=[12,12])
plt.title('N=%d'%len(vdf))
figs['sureClu'] = plt.gcf() 


clu0.to_csv('clu.csv')
pyutil.render__images(figs,)