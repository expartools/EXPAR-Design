#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division

"""
    TmDeltaG.py
    1/10/2008
    ~~~~~~~~
    
    A module for calculating the Tm and DeltaG
"""


__version__ = '1.0'
__author__ = 'Wubin Qu <quwubin@gmail.com> @ZCGLAB @BMI @CHINA'
__license__ = 'GPL'

import re
import math
import LightSeq as LS
import ThermodynamicsParameters as TP

# References

# SantaLucia, J. (1998) A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics, Proceedings of the National Academy of Sci-ences of the United States of America, 95, 1460-1465.

# von Ahsen, N., Wittwer, C.T. and Schutz, E. (2001) Oligonucleotide melting tempera-tures under PCR conditions: Nearest-neighbor corrections for Mg2+, deoxynu-cleotide triphosphate, and dimethyl sulfoxide concentrations with comparison to alternative empirical formulas, Clinical Chemistry, 47, 1956-1961.

def calDeltaG(qseq, hseq, monovalent_conc, divalent_conc, dNTP_conc): 
    """
    Calculate the free Gibbs energy 
    """
    # Initial
    deltaG = 0

    #Calculate deltaH and deltaS
    (deltaH, deltaS) = calDeltaHS(qseq, hseq)

    # Calculate the free Gibbs energy
    tao = 273.15 + 37 # Constant temperature tao in Kelvin

    length = len(hseq)

    # Many thanks for the anonymous referee who help me fix the bug in last version.
    monovalent_conc = monovalent_conc + divalentToMonovalent(divalent_conc, dNTP_conc)
    monovalent_conc = monovalent_conc / 1000

    deltaS_adjust = deltaS + 0.368 * (length - 1) * math.log(monovalent_conc, math.e)

    deltaG = (deltaH * 1000 - tao * deltaS_adjust) / 1000
    return deltaG

def calDeltaHS(qseq, hseq): 
    """
    Calculate deltaH and deltaS
    """
    # Initial
    deltaH = 0
    deltaS = 0

    qseq = re.sub('[atcgn]+', '', qseq)
    hseq = re.sub('[atcgn]+', '', hseq)

    init_begin = 'init' + hseq[0]
    init_end = 'init' + hseq[-1]

    initH = TP.dH_full[init_begin] + TP.dH_full[init_end]
    initS = TP.dS_full[init_begin] + TP.dS_full[init_end]

    deltaH = deltaH + initH
    deltaS = deltaS + initS

    # Get the complement of a sequence
    hseq_complement = LS.seq(hseq).complement

    # Calculate the free Gibbs energy
    length = len(hseq)
    for i in range(length):
        if i > length - 2:
            break
        dinuc = qseq[i : (i+2)] + hseq_complement[i : (i+2)]
        if TP.dH_full.has_key(dinuc) and TP.dS_full.has_key(dinuc):
            deltaH = deltaH + TP.dH_full[dinuc]
            deltaS = deltaS + TP.dS_full[dinuc]
    return (deltaH, deltaS)

def calTm(qseq, sseq, monovalent_conc, divalent_conc, \
          oligo_conc, dNTP_conc):
    """ Calculate Tm value of amplicon"""

    length = len(qseq)

    oligo_conc = oligo_conc / 1000000000

    (delta_H, delta_S) = calDeltaHS(qseq, sseq)
    delta_H = delta_H * 1000

    # Many thanks for the anonymous referee who help me fix the bug in last version.
    monovalent_conc = monovalent_conc + divalentToMonovalent(divalent_conc, dNTP_conc)
    monovalent_conc = monovalent_conc / 1000

    delta_S = delta_S + 0.368 * (length - 1) * math.log(monovalent_conc, math.e)

    Tm = delta_H / (delta_S + 1.987 * math.log(oligo_conc / 4, math.e)) - 273.15

    return Tm

def divalentToMonovalent(divalent, dntp):
    if divalent==0: 
        dntp = 0

    if divalent < 0 or dntp <0: 
        print 'Error'
        exit()

    if divalent < dntp:
        # According to theory, melting temperature doesn't depend on divalent cations
        divalent = dntp

    return 120 * (math.sqrt(divalent - dntp))



def main():
    qseq = 'ACCACAGGATTCACTTGATGAACGA'
    sseq = qseq
    mono = 50
    diva = 1.5
    oligo = 50
    dntp = 0.25
    Tm = calTm(qseq, sseq, mono, diva, oligo, dntp)
    print '%.1f' % Tm
    qseq = 'CACAGGATTCACTTGATGAACGA'
    Tm = calTm(qseq, qseq, mono, diva, oligo, dntp)
    print '%.1f' % Tm
    qseq = 'CAGGATTCACTTGATGAACGA'
    Tm = calTm(qseq, qseq, mono, diva, oligo, dntp)
    print '%.1f' % Tm
    qseq = 'GGATTCACTTGATGAACGA'
    Tm = calTm(qseq, qseq, mono, diva, oligo, dntp)
    print '%.1f' % Tm
    qseq = 'ATTCACTTGATGAACGA'
    Tm = calTm(qseq, qseq, mono, diva, oligo, dntp)
    print '%.1f' % Tm
    qseq = 'TCACTTGATGAACGA'
    Tm = calTm(qseq, qseq, mono, diva, oligo, dntp)
    print '%.1f' % Tm
    qseq = 'ACTTGATGAACGA'
    Tm = calTm(qseq, qseq, mono, diva, oligo, dntp)
    print '%.1f' % Tm
    qseq ='AGTTCCAG'
    sseq='CTGGAACTGGGGGACTCCTGGAACT'
    #Tm = calTm(qseq,sseq,1000,0,1000,0.25)
    [h,s]=calDeltaHS(qseq, sseq)
    print 'dH and dS',h,s
    print 'tm is',  calTm(qseq,'CTGGAACTTC',1000,0,1000,0.25)
    #print Tm

if __name__ == '__main__':
    main()
