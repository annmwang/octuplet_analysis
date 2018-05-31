"""
nathan.py :: A script to sync up raw GBT and FIT data.
The input data is assumed to follow specific formats,
so the suggestions is to run something like:

./StereoRoadAnalysis.x -i ${DATA} -o hist_${RUN}.root -p ${PDO} -n | tee nathan.txt
grep ALL nathan.txt | tr -d "A-Z" > nathan_all.txt
grep FIT nathan.txt | tr -d "A-Z" > nathan_fit.txt
grep GBT nathan.txt | tr -d "A-Z" > nathan_gbt.txt
wc -l nathan*.txt
python decodeGBT_32bit.py -i Run${RUN}_GBT_raw.dat -o ignore.txt -d nathan_gbt.txt                   > nathan_gbt_raw.txt
python decodeFIT_uvr.py   -i Run${RUN}_FIT_raw.dat -o ignore.txt -d nathan_fit.txt --sl -f -r ${RUN} > nathan_fit_raw.txt
# Be sure to remove Ann's start/end emojis from the txt files!
python nathan.py

Its a lot, I know. Thanks for your patience.
"""

import os
import sys
import time

def main():

    packet_len_fit = 14
    packet_len_gbt = 121

    raw_fit = open("nathan_fit_raw.txt").readlines()
    raw_gbt = open("nathan_gbt_raw.txt").readlines()
    events  = open("nathan_all.txt")    .readlines()

    start_fit  = 0
    start_gbt  = 0
    start_time = time.time()

    tag = "Event %i "
    outdir = "nathan"
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    for i,ev in enumerate(events):

        if i % 100 == 0 and i > 0:
            progress(time.time() - start_time, i, len(events))
        if not ev:
            continue
        if ev.startswith("#"):
            continue
        ev = ev.strip()
        if not ev:
            continue

        splitted = ev.split()
        gbt  = int(splitted[0])
        fits = sorted([int(sp) for sp in splitted[1:]])

        # get the gbt packet and write to file
        outgbt = open(os.path.join(outdir, "gbt_%08i.txt" % (gbt)), "w")
        lines_gbt, end_gbt = find(raw_gbt, tag % gbt, packet_len_gbt, start_gbt)
        for li in lines_gbt:
            outgbt.write(li)
        outgbt.close()
        start_gbt = end_gbt

        # get the gbt packet, format it, and write to file
        outsim = open(os.path.join(outdir, "sim_%08i.txt" % (gbt)), "w")
        for li in format(lines_gbt):
            outsim.write(li)
        outsim.close()

        # get the fit packets and write to file
        for fit in fits:
            lines_fit, end_fit = find(raw_fit, tag % fit, packet_len_fit, start_fit)
            outfit = open(os.path.join(outdir, "gbt_%08i_fit_%08i.txt" % (gbt, fit)), "w")
            for li in lines_fit:
                if "Event" in li:
                    li = "%s (GBT Event %i)\n" % (li.strip(), gbt)
                outfit.write(li)
            outfit.close()
            start_fit = end_fit

    print

def find(lines, tag, length, start=0):
    for i in xrange(start, len(lines)):
        if tag in lines[i]:
            return lines[i : i+length], i

def format(gbtpacket):

    # ignoring the header makes this easier
    outpk = [ gbtpacket[0] ]
    gbtpk = [ li.strip() for li in gbtpacket[1:]]

    # loop
    for i,li in enumerate(gbtpk):

        new_li = li
        
        # get the BC of the other ADDC
        if i%4 == 0 and i%8 < 4:
            new_li = new_li[:4] + gbtpk[i+4][4:]
        
        # add the fifo trailer
        new_li = new_li + (" 21" if i%8 < 4 else " 20")
        
        outpk.append(new_li+"\n")

    return outpk

def progress(time_diff, nprocessed, ntotal):
    nprocessed, ntotal = float(nprocessed), float(ntotal)
    rate = (nprocessed+1)/time_diff
    msg = "\r > %6i / %6i | %2i%% | %8.2fHz | %6.1fm elapsed | %6.1fm remaining"
    msg = msg % (nprocessed, ntotal, 100*nprocessed/ntotal, rate, time_diff/60, (ntotal-nprocessed)/(rate*60))
    sys.stdout.write(msg)
    sys.stdout.flush()

if __name__ == "__main__":
    main()
