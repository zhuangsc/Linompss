#!/bin/bash

ITER=$1
CGITER=$2
MODE=$3

if [ "$#" -ne 3 ]
then
	echo "USAGE: $0 [ITREF iter] [CG iter] [visualize]"
	exit 1
fi

VERSIONS=("CGPROF")
#VERSIONS=("CGPROF" "CGMOD1")
WARMUPS=(0 5 10 20 30)
DISCARDS=(0 1 3 5 10 20 30 40 50 "dynamic")

DIM=32000
BS=500
RHS=1
RESOLUTION=$DIM

rm -fv *.pdf *.dat *.gp error.log

INVALID=0
for ver in ${VERSIONS[@]}
do
	for warmup in ${WARMUPS[@]}
	do
		PREFIX0="${ver}_w${warmup}_r${RESOLUTION}"
		GP=${PREFIX0}.gp
		RES_DAT=${PREFIX0}_residual.dat
		ACC_DAT=${PREFIX0}_acceptance.dat
		ORT_DAT=${PREFIX0}_orth.dat

		if [ -e "$GP" ]
		then
			rm -fv $GP
		fi
		if [ -e "$RES_DAT" ]
		then
			rm -fv $RES_DAT
		fi
		if [ -e "$ACC_DAT" ]
		then
			rm -fv $ACC_DAT
		fi
		if [ -e "$ORT_DAT" ]
		then
			rm -fv $ORT_DAT
		fi

		for discard in ${DISCARDS[@]}
		do
			TMP=tmp.dat
			PREF=${ver}_w${warmup}_r${discard}
			DIR=${PREF}
			LOG=${DIR}/${PREF}.log
			echo "Processing: $LOG"
			SANITY=`ls -l $LOG | awk '{print $5}'`
			if [ "$SANITY" -eq 0 ]
			then
				echo "$LOG is invalid!" >> error.log
				((INVALID++))
			fi
			iter=0
			for i in `seq 1 $ITER`
			do
				next=$(($iter+1))
				NR_ITER=`grep -n "^ITREF i: $iter\>" $LOG | awk 'BEGIN{FS=":"} {print $1}'`
				NR_NEXT=`grep -n "^ITREF i: $next\>" $LOG | awk 'BEGIN{FS=":"} {print $1}'`
				if [ -n "$NR_NEXT" ] 
				then
					awk -v nr_iter=$NR_ITER -v nr_next=$NR_NEXT 'NR>=nr_iter && NR<=nr_next {print $0}' $LOG > $TMP
				else
					awk -v nr_iter=$NR_ITER 'NR>=nr_iter {print $0}' $LOG > $TMP
				fi

				### RES_DAT
				awk '/CGAS/ {print $8}' $TMP >> ${RES_DAT}
#				MISSING=`grep "CGAS" $TMP|wc -l`
#				GAP=`echo "$CGITER-$MISSING"|bc`
#				for gap in `seq 1 $GAP`
#				do
#					echo "IGNORE" >> ${RES_DAT}
#				done

				## ACC_DAT
				for j in `seq 1 ${warmup}`
				do
					echo "1.0" >> ${ACC_DAT}
				done
				awk '/acc\/tol/ {print $2}' $TMP >> ${ACC_DAT}
#				MISSING=`grep "acc\/tol" $TMP|wc -l`
#				GAP=`echo "$CGITER-$MISSING"|bc`
#				for gap in `seq 1 $GAP`
#				do
#					echo "IGNORE" >> ${ACC_DAT}
#				done

				### ORT_DAT
				awk '/orth/ {print $2}' $TMP >> ${ORT_DAT}
#				MISSING=`grep "orth" $TMP|wc -l`
#				GAP=`echo "$CGITER-$MISSING"|bc`
#				for gap in `seq 1 $GAP`
#				do
#					echo "IGNORE" >> ${ORT_DAT}
#				done

				((iter++))
			done
			echo -e "\n" >> ${RES_DAT}
			echo -e "\n" >> ${ACC_DAT}
			echo -e "\n" >> ${ORT_DAT}
		done
	done
done

rm -fv $TMP

BASELINE_DIR="baseline/CGPROF_baseline.log"
SANITY=`ls -l ${BASELINE_DIR} | awk '{print $5}'`
if [ "$SANITY" -eq 0 ]
then
	echo "${BASELINE_DIR} is invalid!" >> error.log
	((INVALID++))
fi

if [ "$INVALID" -gt 0 ]
then
	echo "WARNING: Invalid logs found, Aborting"
	exit 3
fi

if [ "$MODE" -eq 0 ]
then
	echo "Skip visualization"
	rm -fv *.pdf *.dat *.gp
	exit 2
fi

ITER2=`echo "$ITER^2"|bc -l`
DISCARDS_INFO_RESIDUAL=("discard-0" "discard-1" "discard-3" "discard-5" "discard-10" "discard-20" "discard-30" "discard-40" "discard-50" "dynamic" "baseline")
DISCARDS_INFO_ACCEPTANCE=("discard-0" "discard-1" "discard-3" "discard-5" "discard-10" "discard-20" "discard-30" "discard-40" "discard-50" "dynamic" "baseline")

GP_COLOR="/home/bsc28/bsc28687/git_rep/gnuplot-colorbrewer/qualitative/Accent.plt"
GP_FORMAT=jpeg
GP_EXT=jpeg
##################################################################
##                     Residual plot
##################################################################
for ver in ${VERSIONS[@]}
do
	for warmup in ${WARMUPS[@]}
	do
		PREFIX0="${ver}_w${warmup}_r${RESOLUTION}"
		RES_GP=${PREFIX0}_residual.gp
		RES_DAT=${PREFIX0}_residual.dat
		NR_ITER=`grep -n "^ITREF i: 0\>" ${BASELINE_DIR} | awk 'BEGIN{FS=":"} {print $1}'`
		awk -v nr_iter=$NR_ITER 'NR>=nr_iter {print $0}' ${BASELINE_DIR}|awk '/CGAS/ {print $8}' >> ${RES_DAT}

		GP="set title \"Residual ITER: ${ITER} WARMUP: ${warmup}\"\n"
		GP+="set xlabel \"Iterations\"\n"
		GP+="set xrange [0:$ITER2]\n"
		GP+="set xtics 40 nomirror\n"
		GP+="set mxtics 2\n"
		GP+="set ylabel \"log_{10} ||r||\"\n"
		GP+="set yrange [-18:4]\n"
		GP+="set ytics 4\n"
		GP+="set mytics 2\n"
		GP+="set format x \"\"\n"
		GP+="set key below\n"
		GP+="set term ${GP_FORMAT}\n"
		GP+="set output \"${ver}_i${ITER}_w${warmup}_residual.${GP_EXT}\"\n\n"

		GP+="load '${GP_COLOR}'\n"

		GP+="set datafile missing \"IGNORE\"\n"
		GP+="plot "

		idx=0
		for info in ${DISCARDS_INFO_RESIDUAL[@]}
		do
			i=$(($idx+1))
			GP+="\"${RES_DAT}\" index $idx using (log10(\$1)) with lines ls $i title \"${info}\","
			((idx++))
		done
		GP+="\n"

		echo -e ${GP} > ${RES_GP}
		sed -i "s/\(.*\),$/\1/g" ${RES_GP}
		gnuplot $RES_GP
	done
done

##################################################################
##                     Acceptance plot
##################################################################
for ver in ${VERSIONS[@]}
do
	for warmup in ${WARMUPS[@]}
	do
		PREFIX0="${ver}_w${warmup}_r${RESOLUTION}"
		ACC_GP=${PREFIX0}_acceptance.gp
		ACC_DAT=${PREFIX0}_acceptance.dat
		for a in $ITER
		do
			for b in $CGITER
			do
				echo "1.0000" >> ${ACC_DAT}
			done
		done
#		awk '/acc\/tol/ {print $2}' ${BASELINE_DIR} | sed "s/.*/1.0000/g" >> ${ACC_DAT}

		GP="set title \"Acceptance percentage ITER: ${ITER} WARMUP: ${warmup}\"\n"
		GP+="set xlabel \"Iterations\"\n"
		GP+="set xrange [0:$ITER2]\n"
		GP+="set xtics 40 nomirror\n"
		GP+="set mxtics 2\n"
		GP+="set ylabel \"ACCEPTED/TOTAL\"\n"
		GP+="set yrange [1.5:2]\n"
		GP+="set format x \"\"\n"
		GP+="set key below\n"
		GP+="set term ${GP_FORMAT}\n"
		GP+="set output \"${ver}_i${ITER}_w${warmup}_acceptance.${GP_EXT}\"\n\n"

		GP+="load '${GP_COLOR}'\n"

		GP+="plot "

		idx=0
		for info in ${DISCARDS_INFO_ACCEPTANCE[@]}
		do
			i=$(($idx+1))
			GP+="\"${ACC_DAT}\" index $idx using (log10((\$1)*100)) with lines ls $i title \"${info}\","
			((idx++))
		done
		GP+="\n"

		echo -e ${GP} > ${ACC_GP}
		sed -i "s/\(.*\),$/\1/g" ${ACC_GP}
		gnuplot $ACC_GP
	done
done

##################################################################
##                     Orthogonality plot
##################################################################
for ver in ${VERSIONS[@]}
do
	for warmup in ${WARMUPS[@]}
	do
		PREFIX0="${ver}_w${warmup}_r${RESOLUTION}"
		ORT_GP=${PREFIX0}_orth.gp
		ORT_DAT=${PREFIX0}_orth.dat
		NR_ITER=`grep -n "^ITREF i: 0\>" ${BASELINE_DIR} | awk 'BEGIN{FS=":"} {print $1}'`
		awk -v nr_iter=$NR_ITER 'NR>=nr_iter {print $0}' ${BASELINE_DIR}|awk '/orth/ {print $2}' >> ${ORT_DAT}

		GP="set title \"P orthogonality ITER: ${ITER} WARMUP: ${warmup}\"\n"
		GP+="set xlabel \"Iterations\"\n"
		GP+="set xrange [0:$ITER2]\n"
		GP+="set xtics 40 nomirror\n"
		GP+="set mxtics 2\n"
		GP+="set ylabel \"log_{10} p_{n+1}^T*A*p_{n}\"\n"
#		GP+="set yrange [-18:4]\n"
		GP+="set ytics 4\n"
		GP+="set mytics 2\n"
		GP+="set format x \"\"\n"
		GP+="set key below\n"
		GP+="set term ${GP_FORMAT}\n"
		GP+="set output \"${ver}_i${ITER}_w${warmup}_orth.${GP_EXT}\"\n\n"

		GP+="load '${GP_COLOR}'\n"

		GP+="plot "

		idx=0
		for info in ${DISCARDS_INFO_RESIDUAL[@]}
		do
			i=$(($idx+1))
			GP+="\"${ORT_DAT}\" index $idx using (log10(\$1)) with lines ls $i title \"${info}\","
			((idx++))
		done
		GP+="\n"

		echo -e ${GP} > ${ORT_GP}
		sed -i "s/\(.*\),$/\1/g" ${ORT_GP}
		gnuplot $ORT_GP
	done
done

