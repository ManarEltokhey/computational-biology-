#!/usr/bin/env python
# coding: utf-8

# In[3]:


import pyopenms
from pyopenms import *
import pyopenms as SQ
sum=0
seq=SQ.AASequence.fromString("VAKA")
for aa in seq:
    sum+=aa.getMonoWeight()
print(sum)
seq=SQ.AASequence.fromString("VAKA")
total=seq.getMonoWeight()
print(total)


# In[ ]:




