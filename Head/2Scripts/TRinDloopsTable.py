#!/home/aragret/anaconda3/bin/python

import re
import sys

word = 'Sequence:'

with open(str(sys.argv[1]),'r') as infile, open('../../Body/2Derived/TRinDloops.txt', 'w') as outfile:
    outfile.write('Species\tStart\tEnd\tCopyNumber\tPercentMatches\tConsensus\tRepeatsRegion\n')
    for line in infile:

        # вытаскиваем название
        if word in line:
            a = line[10:].strip()
            name = re.sub(' ', '_', a)

        # ищем строки, где в конце стоит одна из букв ATGC
        # шлифуем от всего лишнего, разделяем по пробелу: до него консенсусная посл., после сам регион повтора
        if re.search('[ATGC]{5}$', line) != None:
            b = re.sub('[\d.]', '', line)
            c = b[13:]
            d = re.split(' ', c)
            tr = '\t'.join(d)
            bb = re.sub('[ATGC]', '', line)
            cc = re.split(' ', bb)

            #вытаскиваем координаты повторов, copy number, percent matches
            dd = cc[0:2]
            StEnd = '\t'.join(dd)
            ff = cc[3:4]
            cnumb = ''.join(ff)
            gg = cc[5:6]
            percent = ''.join(gg)
            outfile.write('{0}\t{1}\t{2}\t{3}\t{4}'.format(name, StEnd, cnumb, percent, tr))
