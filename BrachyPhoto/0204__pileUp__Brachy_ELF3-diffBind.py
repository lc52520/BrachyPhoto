#!/usr/bin/env python
# -*- coding: utf-8 -*-
NCORE=6
import pymisca.ext as pyext
execfile(pyext.base__file('headers/header__import.py'))
figs = pyutil.collections.OrderedDict()
# execfile('/home/feng/meta/header__script2figure.py')
# import synotil.chipShot
keyDF = pyutil.readBaseFile('headers/key_brachy.csv')
# gconf = pyutil.envSource('/home/feng/ref/config/config_Bd21-3.sh',silent=1);
# gconf = pyutil.util_obj(**gconf)
# peak2gene = pyutil.readData('/home/feng/envs/Fig_Brachy/peak2gene.tsv')


plt.rcParams['font.size'] = 20.
plt.rcParams['xtick.labelsize'] = 22.
plt.rcParams['ytick.labelsize'] = 22.
plt.rcParams['axes.titlepad'] = 24.
plt.rcParams['legend.fontsize'] = 22.

# bwCurr = bwMeta.query('runID=="182C"').query('bamFinal.str.contains("TAIR10")')
# bwCurr = bwMeta
# bwCurr.bname=bwCurr.bname.str.split('_').str.get(0)

# bwCurr['bnameShort']=bwCurr.bname.str.split('_').str.get(0)
# bwCurr = bwCurr.query('runID=="188CR"')




# #####
# rnaScore = pyutil.readData('/home/feng/envs/Fig_Brachy/score__rnaseq__AUC.csv')
# xlab = 'transcriptionally perturbed\n over day-night cycle'
# query = 'score__rnaseq__AUC>1.0'
# ind1 = rnaScore.query(query).index


# (pyext.base__file('results/0130__callDiffTarg__Brachy-ELF3/DONE'))

chipClu = pyutil.readBaseFile('results/0130__callDiffTarg__Brachy-ELF3/clu.csv')
chipClu.columns = ['clu']

# peak2gene = pyutil.readBaseFile('results/0130__makePeakWindows__Brachy-ELF3/peak2geneFile.tsv')
# geneDF  = sdio.peak2gene(peak2gene, chipClu)
# indCHIP = rnaIndex = geneDF.query('clu==1').feat_acc.unique() 
# indCHIP = tks.rnaseq.index & indCHIP

chipClu.query('clu==1')
# chipTrack = tks.rnaseq.eval("index.isin(@rnaIndex)")




# chipDF = pyutil.readData('/home/feng/envs/Fig_Brachy/score__1024__chipTarg__ELF3.csv')
# chipDF['isTarg'] = chipDF['score__chipDiff-GT-0dot5andsd__chip-GT-0dot45']
# geneDF = sdio.peak2gene(peak2gene,chipDF.query('isTarg'))
# ylab = 'chipSeq_differentially_bound'
# ind2 = geneDF.feat_acc.unique()



peakAcc = chipClu.query('clu==1').index.str.replace('_\d+$','').unique()
# peakAcc = chipDF.query('isTarg').index.str.replace('_\d+$','').unique()
# bedFile  = '/home/feng/envs/Fig_Brachy/per_score-GT-0dot6188C_RESEQ-combined.bed'
bedFile = pyext.base__file('results/0130__makePeakWindows__Brachy-ELF3/windowFile.bed')
bedFile =  pyutil.grepFileByKeys(bedFile,peakAcc)




# sutil.job__nearAUG
bwCurr= pyutil.readBaseFile('results/0201__dumpMeta__Brachy/bwCurr.csv')
bwCurr = bwCurr.query('bname.str.contains("ELF3",case=0)')
bwFiles = bwCurr.RPKMFile.map(pyext.base__file)
res = sjob.figs__peakBW(bwFiles=bwFiles,
                        center_summit=1,
                       peakFile=bedFile,
                        detailByChip=0,
                      outerRadius=1500,
                      innerRadius=300,
)
figD = res[0]
axs =figD.values()[0].axes
print (axs[1].get_xlim())
# axs[1].set_xlim(-1600,1600)
figs.update(figD)

pyutil.render__images(figs,exts=['png','svg'])