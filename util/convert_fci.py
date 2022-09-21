# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 17:11:31 2019


clean no-physical elements and add symmetry elements

    spin up = 1
    spin down = 0 

Now compatible for all FCIDUMP file

A whole block is missing

@author: izzy_
"""
import sys

tostr = lambda u:u[0]+' '+str(u[1])+' '+str(u[2])+ ' '+str(u[3])+ ' '+str(u[4])+'\n'; 

def nothing(output, isIn,coeff,m,n,p,q):
    
    keyGen = lambda a,b,c,d:str(a)+' '+str(b)+' '+str(c)+' '+str(d)
    key = keyGen(m,n,p,q)
    if isIn.has_key(key) == False:
        output.append(tostr([coeff,m,n,p,q]))
        isIn[key] = True
#main 

filename = sys.argv[1:]

for name in filename:
    
    print name 
    
    input_name = name
    output_name = name.split('_')[0] + "_IN"
    
    fp_out = open(output_name,"w") 
    
    seg = "0   0   0   0"
    raw = [] 
    result  = [] 
    
    spin = ["uu","dd","ud","du","u","d","c"] 
    group = []
    
    mk = False 
    
    with open(input_name) as f:
        for line in f:
            if mk == True:
                temp = line.split()

                #if temp[0]=="0.0000000000000000E+00":
                if temp[1] == "0" and temp[2] == "0" and temp[3] == "0" and temp[4] == "0" and float(temp[0]) == 0.0: 

                    group.append(raw)
                    raw = [] 
                else:
                    temp[1] = int(temp[1])
                    temp[2] = int(temp[2])
                    temp[3] = int(temp[3])
                    temp[4] = int(temp[4])
                    raw.append(temp) 
            if line == ' /\n':
                mk = True 
    
    group.append(raw)
    a = group[:3]
    c = [group[2]]
    b = group[3:]
    group = a + c + b
    
    cand = zip(spin,group) 
    univer = []
    
    isIn = {}
    
    for spin, orb in cand:
        
        output_temp = [] 
        
        if len(spin)==2:
            first_spin = 1 if spin[0] == 'u' else 0
            second_spin = 1 if spin[1] == 'u' else 0
            
            #isIn = {} 
            
            for item in orb:
                coeff = item[0]
                m = item[1]<<1
                n = item[2]<<1
                p = item[3]<<1
                q = item[4]<<1
                
                m = m+first_spin
                n = n+first_spin
                p = p+second_spin
                q = q+second_spin
                
                """
                if m==n and n==p and p==q:
                    nothing(output_temp,isIn, item[0], m,n,p,q)
                """ 
                
                """
                    ud & du section should be in the same. 
                """                
                
                if m!= p and n!=q:
                    
                    nothing(output_temp,isIn, item[0], m,n,p,q)
                    nothing(output_temp,isIn, item[0], n,m,q,p)
                    nothing(output_temp,isIn, item[0], p,q,m,n)
                    nothing(output_temp,isIn, item[0], q,p,n,m)
                             
                if n!=p and m!=q: 
                    
                    nothing(output_temp,isIn, item[0], n,m,p,q)
                    nothing(output_temp,isIn, item[0], p,q,n,m)
                    nothing(output_temp,isIn, item[0], m,n,q,p)
                    nothing(output_temp,isIn, item[0], q,p,m,n)
        
        else:        
            if spin[0] == 'c':
                for item in orb:
                    #ceoff = item[0]
                    output_temp.append(tostr([item[0],0,0,0,0]))  
                    
            else:
                first_spin = 1 if spin[0] == 'u' else 0 
                
                for item in orb: 
                    ceoff = item[0]
                    m = item[1] <<1
                    n = item[2] <<1
                    
                    m = m+ first_spin
                    n = n+ first_spin
                        
                    output_temp.append(tostr([item[0],m,n,0,0]))  
                    if m!=n:
                        output_temp.append(tostr([item[0],n,m,0,0]))     

        univer = univer + output_temp
        
        print (str(len(output_temp)))                
        fp_out.writelines(str(len(output_temp)) +'\n' )
        fp_out.writelines(output_temp)      
    
    fp_out.close()