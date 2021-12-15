#!/usr/bin/env python
# coding: utf-8

# In[2]:


from pyopenms import *
dig = ProteaseDigestion()
dig.getEnzymeName()
entries=[]
f=FASTAFile()
f.load("E:/fci/year 4/computational biology techniques/section/tasks/yeast.fasta",entries)
c=0
while c<len(entries)-1:
    f=AASequence.fromString(entries[c].sequence)
    print("sequence: ",c)
    result = []
    dig.digest(f, result)
    for e in result:
        print(e.toString())
    print(len(result)) 
    c=c+1
peptides=[AASequence.fromString(s.toString()) for s in result]
for peptide in peptides:
    tsg=TheoreticalSpectrumGenerator()
    spec1=MSSpectrum()
    p=Param()
    p.setValue("add_b_ions","false")
    p.setValue("add_metainfo","true")
    tsg.setParameters(p)
    tsg.getSpectrum(spec1,peptide,1,1)
    print("Spectrum 1 of",peptide,"has",spec1.size(),"peaks")
    for ion,peak in zip(spec1.getStringDataArrays()[0],spec1):
        print(ion.decode(),"is generated at m/z",peak.getMZ())


# In[ ]:





# In[ ]:




