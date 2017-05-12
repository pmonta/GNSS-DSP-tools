#!/bin/sh

DATA=$1
DEST_DIR=$2
mkdir -p ${DEST_DIR}

# L1 1584.754875 MHz

<${DATA} packet2wav_3ch 1 | track-gps-l1.py       /dev/stdin 69984000  -9334875  21  2400.0    817.50  >${DEST_DIR}/track-gps-l1-prn21.dat
<${DATA} packet2wav_3ch 1 | track-glonass-l1.py   /dev/stdin 69984000  17245125  -3 -1200.0    362.82  >${DEST_DIR}/track-glonass-l1-m3.dat
<${DATA} packet2wav_3ch 1 | track-galileo-e1b.py  /dev/stdin 69984000  -9334875  24   250.0   2838.00  >${DEST_DIR}/track-galileo-e1b-prn24.dat
<${DATA} packet2wav_3ch 1 | track-beidou-b1i.py   /dev/stdin 69984000 -23656875  34  -600.0    562.20  >${DEST_DIR}/track-beidou-b1i-prn34.dat

# L2 1227.727125 MHz

<${DATA} packet2wav_3ch 2 | track-gps-l2cm.py     /dev/stdin 69984000   -127126  29  1120.0   4208.80  >${DEST_DIR}/track-gps-l2cm-prn29.dat
<${DATA} packet2wav_3ch 2 | track-glonass-l2.py   /dev/stdin 69984000  18272874  -2 -1800.0    470.98  >${DEST_DIR}/track-glonass-l2-m2.dat
<${DATA} packet2wav_3ch 2 | track-glonass-l3i.py  /dev/stdin 69984000 -25702126   9 -1800.0   9429.00  >${DEST_DIR}/track-glonass-l3i-prn9.dat
<${DATA} packet2wav_3ch 2 | track-galileo-e5bi.py /dev/stdin 69984000 -20587126  24   200.0   7919.00  >${DEST_DIR}/track-galileo-e5bi-prn24.dat
<${DATA} packet2wav_3ch 2 | track-beidou-b2i.py   /dev/stdin 69984000 -20587126  14  -600.0   1682.90  >${DEST_DIR}/track-beidou-b2i-prn14.dat

# L5 1191.641625 MHz

<${DATA} packet2wav_3ch 3 | track-gps-l5i.py      /dev/stdin 69984000 -15191625  25 -1600.0   9696.00  >${DEST_DIR}/track-gps-l5i-prn25.dat
<${DATA} packet2wav_3ch 3 | track-galileo-e5ai.py /dev/stdin 69984000 -15191625  24   200.0   7919.00  >${DEST_DIR}/track-galileo-e5ai-prn24.dat
