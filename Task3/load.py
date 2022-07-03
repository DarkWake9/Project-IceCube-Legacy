##########################################################################

###THIS FILE CONTAINS ALL THE MTHODS USED IN TASK/METHOD 2 FOR EASY ACCESS

##########################################################################
import numpy as np
import pandas as pd
import scipy as sp


##########################################################################################################################################################################

#GENERATES A SYNTHETIC CATALOGUE OF RA, DEC COORDINATES UNIFORMLY DISTIBUTED IN THE SAME RANGE AS IN OBSERVED POPULATION

def gen_cat(minra, maxra, mindec, maxdec):

    msra2 = [np.random.uniform(minra, maxra) for i in range(441)]
    msdec2 = [np.random.uniform(mindec, maxdec) for i in range(441)]
    return(msra2, msdec2)


##########################################################################################################################################################################

#METHOD TO REURN THE SPACE ANGLE(in degrees) BETWEEN 2 PAIRS OF RA AND DEC COORDINATE VECTORS, AND IT'S COSINE VECTORS

def hvovec(lon1, lat1, lon2, lat2):

    #Convert decimal degrees to Radians:
    lon1 = np.radians(lon1)
    lat1 = np.radians(lat1)
    lon2 = np.radians(lon2)
    lat2 = np.radians(lat2)

    #Implementing Haversine Formula: 
    dlon = np.subtract(lon2, lon1)
    #dlat = np.subtract(lat2, lat1)

    a = np.add(np.multiply(np.sin(lat1), np.sin(lat2)), np.multiply(np.multiply(np.cos(lat1), np.cos(lat2)), np.cos(dlon)))

    return (np.rad2deg(np.arccos(a)), a)

def hvcvec(lon1, lat1, lon2, lat2):

    #Convert decimal degrees to Radians:
    lon1 = np.radians(lon1)
    lat1 = np.radians(lat1)
    lon2 = np.radians(lon2)
    lat2 = np.radians(lat2)

    #Implementing Haversine Formula: 
    dlon = np.subtract(lon2, lon1)
    #dlat = np.subtract(lat2, lat1)

    a = np.add(np.multiply(np.sin(lat1), np.sin(lat2)), np.multiply(np.multiply(np.cos(lat1), np.cos(lat2)), np.cos(dlon)))

    return a

##########################################################################################################################################################################

#CALCULATES THE SPACE ANGLE VECTORS AND RETURNS THEM BY REPLACING EXTENDED VALUES WITH "NONE"

def t2a2b(icra, icdec, msra, msdec, extensions, which='b'):
    if which == 'both':
        ft2a = []
        ft2b = []
        for x in range(3):
            st2a = []
            st2b = []
            lg = int(len(icra[x])/len(msra))
            p = len(msra)    
            for k in range(lg):
                ilo = icra[x][k * p  :p * k + p]
                ila = icdec[x][k * p  :p * k + p]
                lo =[]
                la = []
                #t2a = []
                #t2b = []
                for j in range(p):#441
                    lo = [msra[(i + j)%p] for i in range(0,p)]
                    la = [msdec[(i + j)%p] for i in range(0,p)]
                    hvs = hvovec(ilo, ila, lo, la)       #2A, 2C
                    st2a.append(hvs[0])
                    st2b.append(hvs[1])
                    if k == lg - 1:
                        st2a[-1][extensions[x]:] = None
                        st2b[-1][extensions[x]:] = None

            l = np.ravel(st2a)
            m = np.ravel(st2b)
            ft2a.append(l)
            ft2b.append(m)
        return(ft2a, ft2b)
    else:
        ft2a = []
        for x in range(3):
            st2a = []
            lg = int(len(icra[x])/len(msra))
            p = len(msra)    
            for k in range(lg):
                ilo = icra[x][k * p  :p * k + p]
                ila = icdec[x][k * p  :p * k + p]
                lo =[]
                la = []
                for j in range(p):#441
                    lo = [msra[(i + j)%p] for i in range(0,p)]
                    la = [msdec[(i + j)%p] for i in range(0,p)]
                    st2a.append(hvcvec(ilo, ila, lo, la))
                    if k == lg - 1:
                        st2a[-1][extensions[x]:] = None

            l = np.ravel(st2a)
            for i in range(len(l)):
                if l[i] == None:
                    l.pop(i)
            ft2a.append(l)
        return ft2a

#CALCULATES THE COS(SPACE ANGLE) VECTORS AND RETURNS THEM BY REMOVING THE EXTENDED VALUES

def t2b(icra, icdec, msra, msdec, extensions):
    p = len(msra)
    cosdisb = []
    ft2a = []
    for x in range(3):
        st2a = []
        lg = int(len(icra[x])/len(msra))
        p = len(msra)    
        for k in range(lg):
            ilo = icra[x][k * p  :p * k + p]
            ila = icdec[x][k * p  :p * k + p]
            lo =[]
            la = []
            for j in range(p):#441
                lo = [msra[(i + j)%p] for i in range(0,p)]
                la = [msdec[(i + j)%p] for i in range(0,p)]
                st2a.append(hvcvec(ilo, ila, lo, la))
                if k == lg - 1:
                    st2a[-1][extensions[x]:] = None

        l = np.ravel(st2a)
        #REMOVING THE ERRONEOUS DATA AS
        for i in range(len(l)):
            if l[i] == None:
                l.pop(i)
        ft2a.append(l)
    return ft2a

##########################################################################################################################################################################

#GENERATES THE HISTOGRAM AND THE BACKGROUND AND THE SIGNAL VALUES FROM THE GIVEN COSINE DISTRIBUTION

def t2bs(ft2b):
    nbins = np.linspace(np.cos(np.radians(7.35)), np.cos(0),  7)

    freq = [np.histogram(ft2b[i], nbins)[0] for i in range(3)]

    freq1 = [freq[i][:-1]for i in range(3)]
    midpt = [(nbins[i+1] + nbins[i])/2.0 for i in range(len(nbins) - 1)]
    midpt = np.array(midpt[:-1])
    mean = [np.dot(freq1[i], midpt)/float(np.sum(midpt)) for i in range(3)]
    signal = [freq[i][-1] for i in range(3)]
    return(mean, signal, freq)

##########################################################################################################################################################################

#TAKES ICRA, ICDEC ANGERR, MSRA, MSDEC, EXTENSIONS AND RETURNS THE no.of MATCHES r

def t2c(icra, icdec, icangerr, msra, msdec, extensions):
    p = len(msra)
    match_count = []
    for x in range(3):
        r = 0
        lg = int(len(icra[x])/len(msra))
        lo = []
        la = []
        for k in range(lg):
            ilo = icra[x][k * p  :p * k + p]
            ila = icdec[x][k * p  :p * k + p]
            lo =[]
            la = []

            for j in range(p):#441
                        lo = [msra[(i + j)%p] for i in range(0,p)]
                        la = [msdec[(i + j)%p] for i in range(0,p)]
                        hvsang = hvovec(ilo, ila, lo, la)[0]
                        if k == lg - 1:
                            hvsang[extensions[x]:] = None
                        for l in range(p):
                            if hvsang[l] != None:
                            #if hvsang[l] != None and l < extensions[x]:
                                if abs(hvsang[l]) <= icangerr[x][k*p + l]:
                                    #r.append([k,j,l, hvsang[l], icangerr[k*p + l]])
                                    r=r + 1
        match_count.append(r)
    return match_count

#####################################################################################

def t2cfromt2a(sp_ang, icangerr, licra, lmsra):
    for x in range(3):
        lg = int(licra[x]/lmsra)
        p = lmsra
        match_count = 0
        lo = []
        la = []
        for k in range(lg):
            for j in range(p):#441
                for l in range(p):
                        if sp_ang[x][k*p+l] != None:
                            if abs(sp_ang[x][k*p+l]) <= icangerr[x][k*p + l]:
                                #r.append([k,j,l, hvsang[l], icangerr[k*p + l]])
                                match_count+=1
    return match_count

##########################################################################################################################################################################

''' THIS METHOD CALCULATES THE HAVERSINE ANGLES (IN BATCHES) AND IF THAT it satisfies $|cos(\theta)-cos(5)|  < 0.5 x (1-cos X)$

    Then the variable bgdcount (which counts the #background events) is increased by 1'''


def hvocos(lon1, lat1, lon2, lat2):
    cos5 = np.cos(np.deg2rad(5))
    #Convert decimal degrees to Radians:
    lon1 = np.radians(lon1)
    lat1 = np.radians(lat1)
    lon2 = np.radians(lon2)
    lat2 = np.radians(lat2)

    #Implementing Haversine Formula: 
    dlon = np.subtract(lon2, lon1)
    #dlat = np.subtract(lat2, lat1)
    #a = np.add(np.multiply(np.sin(lat1), np.sin(lat2)), np.multiply(np.multiply(np.cos(lat1), np.cos(lat2)), np.cos(dlon)))
    a = np.abs(np.subtract(np.add(np.multiply(np.sin(lat1), np.sin(lat2)), np.multiply(np.multiply(np.cos(lat1), np.cos(lat2)), np.cos(dlon))), cos5 * np.ones(len(lat1))))

    return a


def t2d(icra, icdec, icangerr, msra, msdec, extensions):
    canger = [0.5 * (1-np.cos(np.deg2rad(i))) for i in icangerr]
    background = []    
    for x in range(len(icra)):
        lg = int(len(icra[x])/len(msra))
        p = len(msra)
    #background = []
        bgcount = 0
        lo = []
        la = []
        for k in range(lg):
            ilo = icra[x][k * p  :p * k + p]
            ila = icdec[x][k * p  :p * k + p]
            lo =[]
            la = []

            for j in range(p):#441
                        lo = [msra[(i + j)%p] for i in range(0,p)]
                        la = [msdec[(i + j)%p] for i in range(0,p)]
                        hvscos = hvocos(ilo, ila, lo, la)
                        if k == lg - 1:
                            hvscos[extensions[x]:] = None
                        lo = []
                        la = []
                        for l in range(p):
                            if hvscos[l] != None:
                                if hvscos[l] < canger[x][k*p+l]:
                                    bgcount+=1
                                    #background.append([int(k*p+j),int(l)])
        background.append(bgcount)
    return background

#DONOT USE EITHER OF THE TWO METHODS OF t2d AND hvocos WITHOUT THE OTHER

##########################################################################################################################################################################