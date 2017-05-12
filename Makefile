# Examples of acquisition and tracking

DATA=data/gnss-20170427-L1L2L5.pcap
DEST_DIR=gnss-20170427-L1L2L5

all: acquire track

acquire: ${DATA}
	mkdir -p ${DEST_DIR}
	acquire-all.sh ${DATA} ${DEST_DIR}

track: ${DATA}
	mkdir -p ${DEST_DIR}
	track-all-gnss-2017-0427-L1L2L5.sh ${DATA} ${DEST_DIR}

# Download the sky-recording waveform

${DATA}:
	mkdir -p data
	wget -O ${DATA} https://rf-waveforms.s3.amazonaws.com/gnss-20170427-L1L2L5.pcap
