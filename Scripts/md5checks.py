#!/usr/bin/env python
# coding: utf-8

# In[126]:


md5Dict = {}
with open('logs/md5checks.log', "r") as file:
            Lines = file.readlines()
            
            for line in Lines:
                Filename = line.split(":")[0]
                check= line.split(":")[1]
                check = check.strip("\n")
                check = check.strip(" ")
                
                md5Dict[Filename] = check
            
for (file, check) in md5Dict.items():
    if check != 'OK':
        print(f'md5check failed for {file}')
    else:
        print(f'md5check was successful for {file}')
                

            
              

            

