#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pymisca.ext as pyext
execfile(pyext.base__file('headers/header__import.py'))
figs = pyutil.collections.OrderedDict()

# meta = pyutil.readBaseFile('meta/meta_rna.tsv')
meta = pyutil.readBaseFile('results/0201__dumpMeta__Brachy/rnaCurr.csv')

keyDF = pyutil.readBaseFile('headers/key_brachy.csv')
meta.sort_values(["SpecAcc","gtype","Age","light","temp","ZTime_int"],inplace=True)
# mcurr = mcurr0 = meta.query('SpecAcc=="Bd" & Age in ["Wk2","Wk3"]')
mcurr = mcurr0 = meta.query('SpecAcc=="Bd" & Age in ["Wk2","Wk3","Wk4"]')
# rnaDF = pyutil
# rnaseq = pyutil.readData('/home/feng/envs/Fig_Brachy/rna_brachy_log2p1_aligned.pk')
rnaseq = pyutil.readBaseFile('results/0131__dumpDataRNA__Brachy/rna_log2p1_filtered.pk')

# key ='rnaseq'
# assert key in locals().keys(),'Variable not defined:%s'%key

# rnaseq = dfccc

# dfccc.heatmap()



## sutil.qc_Sort(df=dfccc)

tdf  = rnaseq.reindex(columns=mcurr0.index).dropna()
tdf = sutil.meanNorm(tdf)
rnaseq_all  = tdf
# tks = pyutil.util_obj()

# mcurr = meta.query('SpecAcc=="Bd" & gtype == "Bd21" & Age == "Wk3"')
query = 'SpecAcc=="Bd" & light == "SD" & Age == "Wk3"'

mcurr = meta.query(query)
m1 = mcurr.query('gtype=="Bdelf3"')
m2 = mcurr.query('gtype=="Bd21"')
tdf = pyutil.init_DF(
    (rnaseq.reindex(columns = m1.index).values 
       - rnaseq.reindex(columns = m2.index).values ),rowName = rnaseq.index,colName = m1.index)

tdf = scount.countMatrix(tdf,colMeta = m1)
tdf.relabel('ZTime')
tdf.vlim = [-2,2]

tdf.name_ = query.replace('&','_').replace('==','=').replace('"','').replace(' ','')
tdf.name_ += "_gtype=Bdelf3-Bd21"

rnaseq_wk3sd = rnaseq_wk3sd_elf3ko = tdf
assert ~(tdf.isnull().values.any())



# mcurr = meta.query('SpecAcc=="Bd" & gtype == "Bd21" & Age == "Wk3"')
query = 'SpecAcc=="Bd" & light == "SD" & Age == "Wk2"'
mcurr = meta.query(query)
m1 = mcurr.query('gtype=="Bdelf3"')
m2 = mcurr.query('gtype=="Bd21"')
tdf = pyutil.init_DF(
    (rnaseq.reindex(columns = m1.index).values 
       - rnaseq.reindex(columns = m2.index).values ),rowName = rnaseq.index,colName = m1.index)

tdf = scount.countMatrix(tdf,colMeta = m1)
tdf.relabel('ZTime')
tdf.vlim = [-2,2]
tdf.name_ = query.replace('&','_').replace('==','=').replace('"','').replace(' ','')
tdf.name_ += "_gtype=Bdelf3-Bd21"

rnaseq_wk2sd = rnaseq_wk2sd_elf3ko = tdf


query = 'SpecAcc=="Bd" & light == "LD" & Age == "Wk2"'
mcurr = meta.query(query)
m1 = mcurr.query('gtype=="phyC"')
m2 = mcurr.query('gtype=="Bd21"')

mm = m1.reset_index().merge(m2.reset_index(),on='ZTime_int',how='inner',sort=True)

tdf = pyutil.init_DF(
    (rnaseq.reindex(columns = mm.DataAcc_x).values 
       - rnaseq.reindex(columns = mm.DataAcc_y).values ),
    
    rowName = rnaseq.index,colName = mm.DataAcc_x,)

tdf = scount.countMatrix(tdf,colMeta=meta.reindex(mm.DataAcc_x))
tdf.vlim = [-2,2]
tdf.name_ = query.replace('&','_').replace('==','=').replace('"','').replace(' ','')
tdf.name_ += "_gtype=BdphyC-Bd21"
tdf.relabel('ZTime')

rnaseq_wk2ld_phycko = tdf




query = ''
mcurr = meta.query('SpecAcc=="Bd" & Age == "Wk3"')
# mcurr = meta.query('SpecAcc=="Bd" & light == "LD" & Age == "Wk2"')
# m1 = mcurr.query('gtype=="Bdelf3"')
# m2 = mcurr.query('gtype=="Bd21"')
# tdf = pyutil.init_DF(
#     (rnaseq.reindex(columns = m1.index).values 
#        - rnaseq.reindex(columns = m2.index).values ),rowName = rnaseq.index,colName = m1.index)
tdf = rnaseq.reindex(columns = mcurr.index)

tdf = scount.countMatrix(tdf,colMeta = mcurr)
tdf.relabel(['gtype','ZTime'])
tdf.vlim = [-2,2]

tdf.name_ = query.replace('&','_').replace('==','=').replace('"','').replace(' ','')


rnaseq_wk3sd_joint = tdf


query = 'SpecAcc=="Bd" & gtype=="Bd21" & Age == "Wk2"'
mcurr = meta.query(query)
m1 = mcurr.query('light=="LD"')
m2 = mcurr.query('light=="SD"')
tdf = pyutil.init_DF(
    (rnaseq.reindex(columns = m1.index).values 
       - rnaseq.reindex(columns = m2.index).values ),rowName = rnaseq.index,colName = m1.index)

tdf = scount.countMatrix(tdf,colMeta = m1)
tdf.relabel('ZTime')
tdf.vlim = [-2,2]

tdf.name_ = query.replace('&','_').replace('==','=').replace('"','').replace(' ','')
tdf.name_ += '_light=LD-SD'
print tdf.name
rnaseq_wk2wt = tdf


query = 'SpecAcc=="Bd" & gtype=="Bd21" & Age == "Wk3" & runID!="196R"'
mcurr = meta.query(query)
mcurr = mcurr.query( 'ZTime!="ZT22"')
m1 = mcurr.query('light=="LD"')
m2 = mcurr.query('light=="SD"')
tdf = pyutil.init_DF(
    (rnaseq.reindex(columns = m1.index).values 
       - rnaseq.reindex(columns = m2.index).values ),rowName = rnaseq.index,colName = m1.index)

tdf = scount.countMatrix(tdf,colMeta = m1)
tdf.relabel('ZTime')
tdf.vlim = [-2,2]

tdf.name_ = query.replace('&','_').replace('==','=').replace('"','').replace(' ','')
tdf.name_ += '_light=LD-SD'
print tdf.name
rnaseq_wk3wt = tdf

query = 'SpecAcc=="Bd" & gtype=="Bd21" & light == "LD" & runID!="196R"'
mcurr = meta.query(query)
# mcurr = mcurr.query( 'ZTime!="ZT22"')
m1 = mcurr.query('Age=="Wk3"')
m2 = mcurr.query('Age=="Wk2"')
mm = m1.reset_index().merge(m2.reset_index(),on='ZTime_int',how='inner',sort=True)
tdf = pyutil.init_DF(
    (rnaseq.reindex(columns = mm.DataAcc_x).values 
       - rnaseq.reindex(columns = mm.DataAcc_y).values ),
    
    rowName = rnaseq.index,colName = mm.DataAcc_x,)

tdf = scount.countMatrix(tdf,colMeta = m1)
tdf.relabel('ZTime')
tdf.vlim = [-2,2]

tdf.name_ = query.replace('&','_').replace('==','=').replace('"','').replace(' ','')
tdf.name_ += '_Age=Wk3-Wk2'
print tdf.name
rnaseq_wk3ld_wk2ld = tdf


query = 'SpecAcc=="Bd" & Age == "Wk3" & light=="LD" & runID!="196R"'
mcurr = meta.query(query)
m1 = mcurr.query('gtype=="Bdelf3"')
m2 = mcurr.query('gtype=="Bd21"')
tdf = pyutil.init_DF(
    (rnaseq.reindex(columns = m1.index).values 
       - rnaseq.reindex(columns = m2.index).values ),rowName = rnaseq.index,colName = m1.index)

tdf = scount.countMatrix(tdf,colMeta = m1)
tdf.relabel('ZTime')
tdf.vlim = [-2,2]

tdf.name_ = query.replace('&','_').replace('==','=').replace('"','').replace(' ','')
tdf.name_ += '_gtype=Bdelf3-Bd21'
print tdf.name
rnaseq_wk3ld_elf3ko = tdf



query = 'SpecAcc=="Bd" & gtype == "Bd21" & light=="LD" & runID!="196R" '
mcurr = meta.query(query)
m1 = mcurr.query('Age=="Wk3"')
m2 = mcurr.query('Age=="Wk2"')

mm = m1.reset_index().merge(m2.reset_index(),on='ZTime_int',how='inner',sort=True)

tdf = pyutil.init_DF(
    (rnaseq.reindex(columns = mm.DataAcc_x).values 
       - rnaseq.reindex(columns = mm.DataAcc_y).values ),
    
    rowName = rnaseq.index,colName = mm.DataAcc_x,)

# tdf = pyutil.init_DF(
#     (rnaseq.reindex(columns = m1.index).values 
#        - rnaseq.reindex(columns = m2.index).values ),rowName = rnaseq.index,colName = m1.index)

tdf = scount.countMatrix(tdf,colMeta = m1)
tdf.relabel('ZTime')
tdf.vlim = [-2,2]

tdf.name_ = query.replace('&','_').replace('==','=').replace('"','').replace(' ','')
tdf.name_ += '_Age=Wk3-Wk2'
print tdf.name
rnaseq_wk3ld_wk2ld = tdf

query = 'SpecAcc=="Bd"  & light=="LD" & Age== "Wk3" & runID=="196R"'
mcurr = meta.query(query)
m1 = mcurr.query('gtype=="ppd1-1"')
m2 = mcurr.query('gtype=="Bd21"')

mm = m1.reset_index().merge(m2.reset_index(),on='ZTime_int',how='inner',sort=True)

tdf = pyutil.init_DF(
    (rnaseq.reindex(columns = mm.DataAcc_x).values 
       - rnaseq.reindex(columns = mm.DataAcc_y).values ),
    
    rowName = rnaseq.index,colName = mm.DataAcc_x,)

# tdf = pyutil.init_DF(
#     (rnaseq.reindex(columns = m1.index).values 
#        - rnaseq.reindex(columns = m2.index).values ),rowName = rnaseq.index,colName = m1.index)

tdf = scount.countMatrix(tdf,colMeta = m1)
tdf.relabel('ZTime')
tdf.vlim = [-2,2]

tdf.name_ = query.replace('&','_').replace('==','=').replace('"','').replace(' ','')
tdf.name_ += '_gtype=ppd1-Bd21'
print tdf.name
rnaseq_wk3ldppd1_wk3ldwt = tdf

# assert 0



query = 'SpecAcc=="Bd" & Age.str.contains("day")'
mcurr = meta.query(query)
# m1 = mcurr.query('gtype=="Bdelf3"')
# m2 = mcurr.query('gtype=="Bd21"')
# tdf = pyutil.init_DF(
#     (rnaseq.reindex(columns = m1.index).values 
#        - rnaseq.reindex(columns = m2.index).values ),rowName = rnaseq.index,colName = m1.index)
tdf = sutil.meanNorm(scount.countMatrix(rnaseq.reindex(columns=mcurr.index)))
tdf.columns= mcurr.index
tdf = scount.countMatrix(tdf, colMeta = mcurr)

# tdf.relabel('ZTime')
tdf.vlim = [-2,2]

# tdf.name_ = query.replace('&','_').replace('==','=').replace('"','').replace(' ','')
# tdf.name_ += ''
print tdf.name
rnaseq_age_meanNorm = tdf


mcurr=  meta.copy()
query = 'SpecAcc=="Bd" & Age=="Wk4" & runID =="143R" & light=="LD"'
mcurr = mcurr.query(query)

m1 = mcurr.query('gtype=="Bdphyc"')
m2 = mcurr.query('gtype=="Bd21"')

mm = m1.reset_index().merge(m2.reset_index(),on='ZTime_int',how='inner',sort=True)

tdf = pyutil.init_DF(
    (rnaseq.reindex(columns = mm.DataAcc_x).values 
       - rnaseq.reindex(columns = mm.DataAcc_y).values ),
    
    rowName = rnaseq.index,colName = mm.DataAcc_x,)

tdf = scount.countMatrix(tdf,colMeta=meta.reindex(mm.DataAcc_x))
tdf.vlim = [-2,2]
tdf.name_ = query.replace('&','_').replace('==','=').replace('"','').replace(' ','')
tdf.name_ += "_gtype=BdphyC-Bd21"
tdf.relabel('ZTime')

rnaseq_wk4ldphyc_wk4ldwt = tdf


rnaseq_wk2wt.name_ = 'LD/SD_Wk2_wt'
rnaseq_wk3wt.name_ = 'LD/SD_Wk3_wt'
rnaseq_wk2sd_elf3ko.name_ = 'SD_Wk2_elf3-1/wt'
rnaseq_wk3sd_elf3ko.name_ = 'SD_Wk3_elf3-1/wt'
rnaseq_wk2ld_phycko.name_ = 'LD_Wk2_phyC-4/wt'
rnaseq_wk4ldphyc_wk4ldwt.name_ = 'LD_Wk4_phyC-4/wt'
rnaseq_wk3ld_wk2ld.name_ = 'LD_Wk3/Wk2_wt'



inc = '''
GI
PIF2
PIF3
PIF7
FT1
FT2
BdAP1
PRR37
CCA1
ELF4
LUX
RVE7
RVE6
#BdPRR95
BdCO13
BdPRR59
BdGHD7L2
BBX281
BdMAD5
#BdMAD16
#CDF2
#BdLNK4
BdLNK1
'''.strip().splitlines()
inc = [x for x in inc if not x.startswith("#")]
keyDFC = keyDF.query('BioName in @inc')
keyDFTrack = scount.countMatrix(
            keyDFC.get(['BioName']),height=1.0,
        )

# tracks = pyutil.dictFilter(oldd=locals(),keys='''
# rnaseq_wk2wt
# rnaseq_wk3wt
# rnaseq_wk2sd_elf3ko
# rnaseq_wk3sd_elf3ko
# rnaseq_wk2ld_phycko
# rnaseq_wk4ldphyc_wk4ldwt
# rnaseq_wk3ld_wk2ld
# keyDFTrack
# '''.strip().splitlines())
tracks = pyutil.dictFilter(locals(),
                           keys = filter(lambda x:x.startswith('rnaseq'),
                                         locals().keys()
                                        ) + ['keyDFTrack']
                          )

tracks = pyutil.util_obj(**tracks)
np.save('tracks.npy', tracks)

# tracks = 