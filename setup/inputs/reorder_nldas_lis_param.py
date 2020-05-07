from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import os
import netCDF4 as ncd
import math
#from iteration_utilities import duplicates

dataset1 = Dataset('/glade/u/home/ridwan/work/SUMMA/hhdw1/lis_noah_params/param/nldasAttributes_vector.v0.nc','r')
#dataset2 = Dataset('/glade/u/home/ridwan/work/SUMMA/hhdw1/lis_noah_params/param/nldasAttributes_vector.vcheck.nc','r')

hru_old = dataset1.variables['hruId']
hru_old_array = hru_old[:]
hru_old_list = hru_old_array.tolist()

lat_old = dataset1.variables['latitude']
lat_old_array = lat_old[:]
lat_old_list = lat_old_array.tolist()

long_old = dataset1.variables['longitude']
long_old_array = long_old[:]
long_old_list = long_old_array.tolist()


#hru_new = dataset2.variables['hruId']
#hru_new_array = hru_new[:]
#hru_new_list = hru_new_array.tolist()

hru_dim = dataset1.dimensions['hru'].size
hru_new = np.zeros((hru_dim,1))
hru_final = np.zeros((hru_dim,1))
new_id = np.zeros((hru_dim,1))
z=np.zeros((hru_dim,1))
for x in range(0,hru_dim):
  nd=hru_old_list[x]
  z[x]=((nd+1)/464) + (224-nd%224-1)*464
#  hru_new[0]=z
#     print(int(x/2))
  hru_new[x]=z[x]
#kk=np.array(hru_new)
hru_new_list=hru_new.tolist()
#hru_new_list.sort()
#print(hru_new)
#r=np.argsort(kk)
#print(r)  
#for x in range(0,hru_dim):
#  h = hru_new_list.index(z[x])
#  print(h)
#  new_id[x]=math.floor(hru_new[h])
#  if new_id[x]==new_id[x-1]:
#     new_id[x]=new_id[x-1]+1
#  hru_final[x]=h
#hru_final_list=hru_final.tolist()

#hru_final_list.sort( )
#print(hru_final_list)

#unique = []
#for ele in hru_new_list:
#    if ele not in unique:
#        unique.append(ele)
#print(unique)
#print(len(unique))
#print(hru_new_list.sort())
np.savetxt('result.txt', hru_new_list, fmt='%.1f')
hru_new_list.sort()
np.savetxt('result1.txt', hru_new_list, fmt='%.1f')
#lat_new_list = np.zeros((hru_dim,1))
#long_new_list = np.zeros((hru_dim,1))
#for x in range(0,hru_dim):
#    h = hru_new_list[x]
#    h1 = hru_old_list.index(h)
#    lat_new_list[x] = lat_old_list[h1]
#    long_new_list[x] = long_old_list[h1]

#print(lat_new_list)
#print(long_new_list)
a = np.loadtxt("result.txt",skiprows=0) 
print(a)
#r=np.argsort(a)
for x in range(0,hru_dim):
#    h = r[x]
#    thrsh = math.floor(min(a))
    jj = round((a[x]-50)*10)
    
#    if new_id[x] > new_id[x-1]:
#       new_id[x] = new_id[x]
#    else:
#       new_id[x] = new_id[x]+1
    hru_final[x] = jj
hru_final_list=hru_final.tolist()
np.savetxt('result5.txt', hru_final, fmt='%.1f')
unique = []
for ele in hru_final_list:
    if ele not in unique:
        unique.append(ele)
print(unique)
print(len(unique))
#hru_final_list.sort()
#np.savetxt('result3.txt', hru_final_list, fmt='%.1f')
#    h1 = hru_old_list.index(h)
#    lat_new_list[x] = lat_old_list[h1]
#    long_new_list[x] = long_old_list[h1]
#print(r)
#b = np.loadtxt("result1.txt",skiprows=0)
#print(b) 
#file1.close()
fileName = '/glade/u/home/ridwan/work/SUMMA/hhdw1/lis_noah_params/param/nldasAttributes_vector.v6.nc'
root_grp = Dataset(fileName,'w', format='NETCDF4')
#ifile = fileName
#infile = Dataset(ifile,'r+',format='NETCDF4')
hru = root_grp.createDimension('hru',76088)
cond = root_grp.createVariable('hruId',int,('hru')) 
cond.units = '-'
cond.longname = 'hru id'
cond[:] = hru_final

#cond1 = infile.createVariable('longitude',double,('hru'))
#cond1.units = 'decimal degree east'
#cond1.longname = 'Longitude of HRU\'s centriod point'
#cond1 = long_new_list

root_grp.close()
