#!/bin/sh

DATA=$1
DEST_DIR=$2

# L1 1584.754875 MHz

<${DATA} packet2wav_3ch 1 | ./acquire-gps-l1.py /dev/stdin 69984000 -9334875 >${DEST_DIR}/m-gps-l1.dat
<${DATA} packet2wav_3ch 1 | ./acquire-glonass-l1.py /dev/stdin 69984000 17245125 >${DEST_DIR}/m-glonass-l1.dat
<${DATA} packet2wav_3ch 1 | ./acquire-galileo-e1b.py /dev/stdin 69984000 -9334875 >${DEST_DIR}/m-galileo-e1b.dat
<${DATA} packet2wav_3ch 1 | ./acquire-galileo-e1c.py /dev/stdin 69984000 -9334875 >${DEST_DIR}/m-galileo-e1c.dat
<${DATA} packet2wav_3ch 1 | ./acquire-beidou-b1i.py /dev/stdin 69984000 -23656875 >${DEST_DIR}/m-beidou-b1i.dat

# L2 1227.727125 MHz

<${DATA} packet2wav_3ch 2 | ./acquire-gps-l2cm.py /dev/stdin 69984000 -127126 >${DEST_DIR}/m-gps-l2cm.dat
<${DATA} packet2wav_3ch 2 | ./acquire-glonass-l2.py /dev/stdin 69984000 18272874 >${DEST_DIR}/m-glonass-l2.dat
<${DATA} packet2wav_3ch 2 | ./acquire-glonass-l3i.py /dev/stdin 69984000 -25702126 >${DEST_DIR}/m-glonass-l3i.dat
<${DATA} packet2wav_3ch 2 | ./acquire-glonass-l3q.py /dev/stdin 69984000 -25702126 >${DEST_DIR}/m-glonass-l3q.dat
<${DATA} packet2wav_3ch 2 | ./acquire-galileo-e5bi.py /dev/stdin 69984000 -20587126 >${DEST_DIR}/m-galileo-e5bi.dat
<${DATA} packet2wav_3ch 2 | ./acquire-galileo-e5bq.py /dev/stdin 69984000 -20587126 >${DEST_DIR}/m-galileo-e5bq.dat
<${DATA} packet2wav_3ch 2 | ./acquire-beidou-b2i.py /dev/stdin 69984000 -20587126 >${DEST_DIR}/m-beidou-b2i.dat

# L5 1191.641625 MHz

<${DATA} packet2wav_3ch 3 | ./acquire-gps-l5i.py /dev/stdin 69984000 -15191625 >${DEST_DIR}/m-gps-l5i.dat
<${DATA} packet2wav_3ch 3 | ./acquire-gps-l5q.py /dev/stdin 69984000 -15191625 >${DEST_DIR}/m-gps-l5q.dat
<${DATA} packet2wav_3ch 3 | ./acquire-galileo-e5ai.py /dev/stdin 69984000 -15191625 >${DEST_DIR}/m-galileo-e5ai.dat
<${DATA} packet2wav_3ch 3 | ./acquire-galileo-e5aq.py /dev/stdin 69984000 -15191625 >${DEST_DIR}/m-galileo-e5aq.dat
<${DATA} packet2wav_3ch 3 | ./acquire-glonass-l3i.py /dev/stdin 69984000 10383375 >${DEST_DIR}/m-glonass-l3i-ch3.dat
<${DATA} packet2wav_3ch 3 | ./acquire-glonass-l3q.py /dev/stdin 69984000 10383375 >${DEST_DIR}/m-glonass-l3q-ch3.dat
<${DATA} packet2wav_3ch 3 | ./acquire-galileo-e5bi.py /dev/stdin 69984000 15498375 >${DEST_DIR}/m-galileo-e5bi-ch3.dat
<${DATA} packet2wav_3ch 3 | ./acquire-galileo-e5bq.py /dev/stdin 69984000 15498375 >${DEST_DIR}/m-galileo-e5bq-ch3.dat
<${DATA} packet2wav_3ch 3 | ./acquire-beidou-b2i.py /dev/stdin 69984000 15498375 >${DEST_DIR}/m-beidou-b2i-ch3.dat
