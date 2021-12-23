#!/usr/bin/env python
# coding: utf-8

# In[ ]:


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
#searchdb="E:/fci/year 4/computational biology techniques/section/tasks/yeast.fasta"        
protein_id = ProteinIdentification()
protein_id.setIdentifier("IdentificationRun1")
protein_hit = ProteinHit()
protein_hit.setAccession("sp|MyAccession")
protein_hit.setSequence("PEPTIDERDLQMTQSPSSLSVSVGDRPEPTIDE")
protein_hit.setScore(1.0)
protein_hit.setMetaValue("target_decoy", b"target")
protein_id.setHits([protein_hit])
targets=[]
Decoys=[]
FASTAFile().load("E:/fci/year 4/computational biology techniques/section/tasks/yeast.fasta",targets)
decoy=DecoyGenerator()
for entry in targets:
    rec=FASTAEntry(entry)
    rec.identifier="Decoy-"+rec.identifier
    seq=AASequence().fromString(rec.sequence)
    rec.sequence=decoy.reverseProtein(seq).toString()
    Decoys.append(rec)
FASTAFile().store("yeastFileEditig.fasta",targets+Decoys)    


# In[ ]:





# In[ ]:




