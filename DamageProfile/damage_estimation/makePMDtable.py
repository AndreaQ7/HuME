#!/usr/bin/python

# This script uses the files generated in the previous steps to report in a table the information over each microrganism on:
# - Edit distance
# - Damage in 3' and 5' ends
# The script has to be run in the form:
#
# python makePMDtable.py > $table_damage.tsv 


import os
import numpy
import pandas as pd
import sys
import collections
import pandas

# Get NM values
def NMfreq(editDcol):
    dict_edits = {}
    freq_array= []
    for elem in editDcol:
        if elem not in dict_edits.keys():
            dict_edits[elem] = 1
        else:
            dict_edits[elem] += 1
    for k,v in sorted(dict_edits.items()):
        #print(k,v)
        freq_array.append(v)
    #print(freq_array)
    return freq_array


list_orgs= os.walk('/path/to/PMDtools/output/folder')

# Get list of organisms iterating into PMDtools subfolders
def getOrgs(listorgs):
    list_organisms=[]
    for d,p,f in listorgs:
        if 'damageCalculation' in d:
            for elem in os.listdir(d):
                if 'NM.bed' in elem:
                    org = elem.split('pmd1_')[1].split('.NM')[0]
                    if org not in list_organisms:
                        list_organisms.append(org)
    return list_organisms

list_org=getOrgs(list_orgs)


#print(list_org)

orgdict = dict()
#print(list_org)

# Create a dictionary in which each organism is associated to all the samples in which is present.
for org in list_org:
    if 'sp.' not in org:
        species = org.split('.')[1]
        sample = org.split('.')[0]
    else:
        species = '.'.join(org.split('.')[1:3])
        sample = org.split('.')[0]
    if species in orgdict:
        orgdict[species].append(sample)
    else:
        orgdict[species] = [sample]

#print(orgdict)
sorted_orgdict={}
for species,samples in orgdict.items():
    sorted_orgdict[species] = sorted(orgdict[species])

print(sorted_orgdict)
listObj = list()

dict_samples = {}
dict_CT5 = {}
dict_GA3 = {}
dict_delta = {}


# for each couple sample-microrganism extract the edit distance and damage to 5' and 3'
for species,samples in sorted_orgdict.items():
    sample_list=[]
    CT5_list = []
    GA3_list = []
    delta_list = []
    #orgObj = orgAlignment(species)
    for sample in samples:
        path_bed='/path/to/PMDtools/folder/' + species + '/damageCalculation/' + 'pmd1_' + sample + '.' + species + '.NM.bed'
        path_temp='/path/to/PMDtools/folder/'+ species + '/damageCalculation/' + 'PMD_' + sample + '.' + species + '_temp.txt'
    #print(path_bed)
    #print(path_temp)
        with open(path_temp, 'rt') as damage:
            for l in damage.readlines():
                if 'GA3' not in l:
                    l = l.rstrip()
                    #print(species)
                    CT5 = l.split('\t')[1]
                    GA3 = l.split('\t')[21]
                    break
    #print(species)
    #print(CT)
    #print(GA)
        editDcol = []
    #print(name)
        with open(path_bed, 'rt') as bed:
            for l in bed.readlines():
            #print(l)
                l=l.rstrip()
                editDval = l.split(' ')[1]
                #print(editDval)
                editDcol.append(int(editDval))
            freq_array=NMfreq(editDcol)
            if len(freq_array) >1:
                #print(freq_dict)
                #print(path)
                #print(freq_array)
                numEdit=numpy.array(freq_array)
                diffEdit = numpy.diff(numEdit)
                #print(diffEdit)
                negative_val=(diffEdit<0)*diffEdit
                sum_vals= numpy.sum(negative_val)
                #print(sum_vals)
                delta= -(float(numpy.sum(diffEdit)) / numpy.sum(numpy.absolute(diffEdit)))
                #print(org, delta, sep='\t')
                #if delta < 1:
                #print(diffEdit)
                #print(freq_array)
                #print(sample, species, delta, sep='\t')
            elif len(freq_array)==1:
                #print(freq_array)
                delta = 1.0
                #print(org, delta, sep='\t')
            else:
                delta='-'
                #print(org, delta, sep='\t')
        #print(species, sample, delta, CT5, GA3)
        sample_list.append(sample)
        CT5_list.append(str(CT5))
        GA3_list.append(str(GA3))
        delta_list.append(str(delta))
        dict_samples[species] = sample_list
    dict_CT5[species] = CT5_list
    dict_GA3[species] = GA3_list
    dict_delta[species] = delta_list
    #print(species, sample, delta, CT, GA, sep='\t')
    #deltaList.append(delta)
    #CTlist.append(CT)
    #GAlist.append(GA)
#print(dict_samples)

# Create table merging the info contained in previous dictionaries
n = 0
for species in orgdict.keys():
    if n == 0:
        line1 = '\t'.join(['Species', 'Stat', '\t'.join(dict_samples[species])])
        line2 = '\t'.join([species, 'CT5 damage',  '\t'.join(dict_CT5[species])])
        line3 = '\t'.join(['', 'GA3 damage',  '\t'.join(dict_GA3[species])])
        line4 = '\t'.join(['', 'delta',  '\t'.join(dict_delta[species])])
        print(line1, line2, line3, line4, sep='\n')
        n+=1
    else:
        line1 = '\t'.join([species, 'CT5 damage', '\t'.join(dict_CT5[species])])
        line2 = '\t'.join(['', 'GA3 damage', '\t'.join(dict_GA3[species])])
        line3 = '\t'.join(['', 'delta', '\t'.join(dict_delta[species])])
        print(line1,line2,line3, sep='\n')
