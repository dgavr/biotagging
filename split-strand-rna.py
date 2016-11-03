#!/usr/bin/env python

####### Strand splitting script
###### © Daria Gavriouchkina
###### Sauka-Spengler lab

import pysam
import sys

sam=pysam.Samfile("-", "r")
name=sys.argv[1]

sense= pysam.Samfile(name+'_+.sam', "wh", header=sam.header)
antisense= pysam.Samfile(name+'_-.sam', "wh", header=sam.header)

prob_ct=0
n,n_um,n_ps,n_pr,n_ss,n_sr=0,0,0,0,0,0
for fwd in sam.fetch():
    n+=1
    if n%100000==0:
        print "{0} reads ... paired: {1}/{2}, unpaired: {3}/{4}".format(n,n_ps,n_pr,n_ss,n_sr,n_um)
    if fwd.is_unmapped:
        n_um+=1
        continue
    if fwd.is_proper_pair:
        rev=sam.next()
        if fwd.is_read1 and rev.is_read2:
            if not fwd.is_reverse and rev.is_reverse:
                n_ps+=1
                #print "sens +"
                sense.write(fwd)
                sense.write(rev)
            elif fwd.is_reverse and not rev.is_reverse:
                n_pr+=1
                #print "sens -"
                antisense.write(fwd)
                antisense.write(rev)
        else:
            #print "We have a problem!!!"
            prob_ct+=1
    else:
        if fwd.is_reverse:
            n_sr+=1
            antisense.write(fwd)
        else:
            n_ss+=1
            sense.write(fwd)


print "unmapped: {0}".format(n_um)
print "paired sens: {0}".format(n_ps)
print "paired antisens: {0}".format(n_pr)
print "single sens: {0}".format(n_ss)
print "single antisens: {0}".format(n_sr)
print "problems=", prob_ct
