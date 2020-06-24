#!/bin/bash

USAGE="$0 [config file] [log tarball]"

if [ $# -ne 2 ]
then
	echo -e $USAGE
	exit
fi

if [ ! -e $1 ]
then
	echo "Configureation file doesn't exist, exit."
	exit
else
	CONFIG=$1
fi

if [ ! -e $2 ]
then
	echo "Tarball doesn't exist, exit."
	exit
else
	TAR=$2
fi

##################################################################
### 	Extract information from the configuration file
##################################################################
source $CONFIG

###parsing BIN
BIN_NAME=`basename $LO_BIN`

### parsing number of processes
if [[ $LO_PROC == *..* ]]
then
	PROC=`eval echo {$LO_PROC}`
else
	PROC=($LO_PROC)
fi

### parsing number of block sizes
LO_MAT_B=($LO_MAT_B)


TMP_DIR=tmp_${BIN_NAME}
TAR_NAME=`basename ${TAR}`
if [ -e ${TMP_DIR} ]
then
	rm -R ${TMP_DIR}
fi
mkdir ${TMP_DIR}
cp -f ${TAR} ${TMP_DIR}
tar -xf ${TMP_DIR}/${TAR_NAME} -C ${TMP_DIR}

DAT_FILE=tmp_${BIN_NAME}.dat
GP_FILE=tmp_${BIN_NAME}.gp
DAT_PATH=${TMP_DIR}/${DAT_FILE}
GP_PATH=${TMP_DIR}/${GP_FILE}
touch ${DAT_PATH}
touch ${GP_PATH}

for BS in ${LO_MAT_B[@]}
do
	for proc in ${PROC[@]}
	do
		STATS=`ls ${TMP_DIR}|grep "${BIN_NAME}_P${proc}_.*_${BS}_.*\.log"`
		TIME=`awk '/time:/ {print $2}' ${TMP_DIR}/${STATS}` 
		TIME=`echo "${TIME}/1000"|bc -l`
		echo " $proc $TIME" >> ${DAT_PATH}
	done
	echo -e "\n\n" >> ${DAT_PATH}
done

GP="set title \"Performance: ${BIN_NAME}_${LO_INFO}\"\n"
GP+="set xlabel \"Threads\"\n"
GP+="set ylabel \"Execution Time (ms)\"\n"
GP+="set key below\n"
GP+="set term pdfcairo\n"
GP+="set output \"${BIN_NAME}.pdf\"\n"
GP+="plot "

idx=0
for BS in ${LO_MAT_B[@]}
do
	GP+="\"${DAT_PATH}\" index $idx using 2:xtic(1) with linespoints pointtype 7 title \"b=$BS\","
	((idx++))
done
GP+="\n"
echo -e ${GP} > ${GP_PATH}
###Lose the tailing ","
sed -i "s/\(.*\),$/\1/g" ${GP_PATH}
gnuplot ${GP_PATH}

SPEED_FILE=tmp_speed_${BIN_NAME}.dat
SPEED_PATH=${TMP_DIR}/${SPEED_FILE}
touch ${SPEED_PATH}
GP_FILE=tmp_speed_${BIN_NAME}.gp
GP_PATH=${TMP_DIR}/${GP_FILE}
touch ${GP_PATH}

for BS in ${LO_MAT_B[@]}
do
	num=0
	for proc in ${PROC[@]}
	do
		STATS=`ls ${TMP_DIR}|grep "${BIN_NAME}_P${proc}_.*_${BS}_.*\.log"`
		TIME=`awk '/time:/ {print $2}' ${TMP_DIR}/${STATS}` 
		if [ "$num" -eq 0 ]
		then
			BASE=$TIME
		fi
		((num++))
		TIME=`echo "$BASE/$TIME" | bc -l` 
		echo " $proc $TIME" >> ${SPEED_PATH}
	done
	echo -e "\n\n" >> ${SPEED_PATH}
done

GP="set title \"Speed-up: ${BIN_NAME}_${LO_INFO}\"\n"
GP+="set xlabel \"Threads\"\n"
GP+="set yrange [1:]\n"
GP+="set ylabel \"Speedup\"\n"
GP+="set key below\n"
GP+="set term pdfcairo\n"
GP+="set output \"${BIN_NAME}_speedup.pdf\"\n"
GP+="plot "

idx=0
for BS in ${LO_MAT_B[@]}
do
	GP+="\"${SPEED_PATH}\" index $idx using 2:xtic(1) with linespoints pointtype 7 title \"b=$BS\","
	((idx++))
done
GP+="\n"
echo -e ${GP} > ${GP_PATH}
###Lose the tailing ","
sed -i "s/\(.*\),$/\1/g" ${GP_PATH}
gnuplot ${GP_PATH}
exit
rm -R $TMP_DIR

#		RELERR=`awk '/relerr:/ {print $2}' ${TMP_DIR}/${STATS}`
#		if [ -z "$RELERR" ] 
#		then
#			echo "No RELERR found at ${STATS}"
#			exit
#		fi
#		RELERR=`echo ${RELERR} | sed -e 's/[eE]+*/\\*10\\^/'`
#		THRESH=`echo ${THRESH} | sed -e 's/[eE]+*/\\*10\\^/'`
#		RELERR=`echo "-${RELERR}+${THRESH}"|bc -l`
#		if [[ $RELERR == -* ]]
#		then
#			echo "Suspicious result!!!"
#		fi

