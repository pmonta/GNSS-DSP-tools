FILE_L1=/mnt/b/gps-2014-04/gps9006-l1.dat

acquire-l1:
	echo gps
	<${FILE_L1} packet2wav | ./acquire-gps-l1.py /dev/stdin 68873142.857 -8662285.714
	echo glonass
	<${FILE_L1} packet2wav | ./acquire-glonass-l1.py /dev/stdin 68873142.857 17917714.286
	echo galileo
	<${FILE_L1} packet2wav | ./acquire-galileo-e1b.py /dev/stdin 68873142.857 -8662285.714
	echo beidou
	<${FILE_L1} packet2wav | ./acquire-beidou-b1i.py /dev/stdin 68873142.857 -22984285.714

FILE_L2=/mnt/b/gps-2014-04/gps9007-l2.dat

acquire-l2:
	echo gps-l2cm
	<${FILE_L2} packet2wav | ./acquire-gps-l2cm.py /dev/stdin 68873142.857 -12116571.429
	echo glonass
	<${FILE_L2} packet2wav | ./acquire-glonass-l2.py /dev/stdin 68873142.857 6283428.571
