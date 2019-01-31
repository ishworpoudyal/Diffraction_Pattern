#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 10:50:35 2017

@author: ipoudyal
"""

#this is a changed version according to pdb file explanation of
### www.cgl/ucsf/edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
#import argparse


#parser = argparse.ArgumentParser()
#parser.add_argument("--file","-f",type=str, required=True)
#args = parser.parse_args()

X = []
Y = []
Z = []
Type = []


pdb = '5llw.pdb'
with open(pdb,'r') as pdbfile:
    for line in pdbfile:
        if (line[:4] == "ATOM" ):
           X.append(float(line[31:38].strip()))
           Y.append(float(line[39:46].strip()))
           Z.append(float(line[47:54].strip()))
           Type.append(line[77:78].strip())

NATOMS = len(X)
#print(Type)
#sys.exit()
XA = sum(X)/len(X)
YA = sum(Y)/len(Y)
ZA = sum(Z)/len(Z)

XX = []
YY = []
ZZ = []

for i in range(NATOMS):
    XX.append(X[i]-XA)
    YY.append(Y[i]-YA)
    ZZ.append(Z[i]-ZA)

XMIN = min(XX)
XMAX = max(XX)
YMIN = min(YY)
YMAX = max(YY)
ZMIN = min(ZZ)
ZMAX = max(ZZ)

length = max(XMAX-XMIN,YMAX-YMIN,ZMAX-ZMIN)

print('%d ATOMS READ FROM PDB-FILE'%NATOMS)
print(" ")
print("TRANSLATION VECTOR TO BRING TO ORIGIN")
print ('X,Y,Z : %.4f %.4f %.4f'%(-XA,-YA,-ZA))
print(" ")
print("DIMENSION OF BOX TO HOLD ONE MOLECULE:")
print('X: FROM %.4f, TO %.4f'%(XMIN,XMAX))
print('Y: FROM %.4f, TO %.4f'%(YMIN,YMAX))
print('Z: FROM %.4f, TO %.4f'%(ZMIN,ZMAX))

print("Box size for this pdb is %.1f A" %length)
    

