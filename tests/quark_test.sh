#! /bin/bash

###########################################
## Parse the configuration file
## and gather the output log files
##########################################

##################################################################
####	Modify the variables below to set up the script properly
##################################################################

OMPSS_PROC="NX_PES"

###################################End###############################

USAGE="$0 [config file]"

if [ "$#" -ne 1 ]
then
	echo -e $USAGE
	exit 1
fi

if [ ! -e "$1" ]
then
	echo "Configuration file not found, exit"
	exit 2
fi

##################################################################
### 	Extract information from the configuration file
##################################################################
source $1
if [ -z "$LO_BIN" ]
then
	echo "LO_BIN empty, cannot proceed"
	exit 3
fi

BIN_NAME=`basename $LO_BIN`

### Extract node infomation
INFO=`uname -n -m -o`
INFO="`date` $INFO"

DATE=`date +%F`

### parsing MACHINE
if [[ $LO_MACHINE != "native" && $LO_MACHINE != "MN" ]]
then
	echo "MACHINE should set to either \"native\" or \"MN\""
	exit
fi

if [[ $LO_MACHINE == "MN" ]] && [ -z $LO_WALLTIME ]
then
	echo "Wall time needed on MN3"
	exit
fi

### parsing number of processes
if [[ $LO_PROC == *..* ]]
then
	LO_PROC=`eval echo {$LO_PROC}`
else
	LO_PROC=($LO_PROC)
fi

### parsing number of block sizes
LO_MAT_B=($LO_MAT_B)

#############################################################
### Generating the parameter sequence per application
#############################################################

if [[ ${BIN_NAME} == *gemm_quark ]]
then
	PARAM_SEQ=(CORES LO_SCHED LO_RUNTIME LO_FUNC LO_TRANSA LO_TRANSB LO_MAT_M LO_MAT_N LO_MAT_K LO_ALPHA LO_LDA LO_LDB LO_BETA LO_LDC BS LO_REP)
elif [[ ${BIN_NAME} == *potrf_quark ]]
then
	PARAM_SEQ=(CORES LO_SCHED LO_RUNTIME LO_FUNC LO_MAT_N LO_LDA BS LO_REP)
elif [[ ${BIN_NAME} == *syrk_quark ]]
then
	PARAM_SEQ=(CORES LO_SCHED LO_RUNTIME LO_FUNC LO_ALPHA LO_BETA LO_MAT_N LO_MAT_K LO_LDA LO_LDC BS LO_REP)
else
	echo "Invalid program"
	exit 1
fi

ii=0
for CORES in ${LO_PROC[@]}
do
	for BS in ${LO_MAT_B[@]}
	do
		for param in ${PARAM_SEQ[@]}
		do
			if [ -z "${!param}" ]
			then
				echo "${BIN_NAME}: Parameter(s) missing or incorrect. ${param}"
				exit 3
			fi

			PARAM=${!param}

			if [[ $PARAM == yes ]] || [[ $PARAM == csr ]] || [[ $PARAM == upper ]] || [[ $PARAM == right ]]
			then
				PARAM=1
			fi

			if [[ $PARAM == no ]] || [[ $PARAM == csc ]] || [[ $PARAM == lower ]] || [[ $PARAM == left ]]
			then
				PARAM=0
			fi

			PARAM=`echo "$PARAM" | tr -s " "`
			PARAM_ARRAY="$PARAM_ARRAY $PARAM"
		done
	CMD_POOL[$ii]=$PARAM_ARRAY
	((ii++))
	PARAM_ARRAY=""
	done
done

#echo "${CMD_POOL[@]}"
#############################################################
### Set up prefixes for the shell script
#############################################################

INFO_PREFIX="#########################################\n"
INFO_PREFIX+="### Automatically generated script\n"
INFO_PREFIX+="### For $LO_MACHINE\n"
INFO_PREFIX+="#########################################\n"

NATIVE_PREFIX=$INFO_PREFIX
NATIVE_PREFIX+="#!/bin/bash\n\n\n"
NATIVE_PREFIX+='export NX_ARGS="--smp-workers=1 $NX_ARGS"'

MN_PREFIX=$INFO_PREFIX
MN_PREFIX+="#!/bin/bash\n"
MN_PREFIX+="#\n"
MN_PREFIX+="#BSUB -n 4\n"
MN_PREFIX+="#BSUB -oo output_${BIN_NAME}_%J.out\n"
MN_PREFIX+="#BSUB -eo output_${BIN_NAME}_%J.err\n"
MN_PREFIX+="#BSUB -J ${BIN_NAME}\n"
MN_PREFIX+="#BSUB -W $LO_WALLTIME\n"
MN_PREFIX+="#BSUB -x\n"
MN_PREFIX+="ulimit -c unlimited"
MN_PREFIX+="\n\n"
MN_PREFIX+='export NX_ARGS="--smp-workers=1 $NX_ARGS"'

NATIVE_JOB="${BIN_NAME}.sh"
MN_JOB="${BIN_NAME}.sh"

#############################################################
### Writing temporary shell script
############################################################
TMP_JOB=tmp.sh
if [ -e $TMP_JOB ]
then
	rm -f $TMP_JOB
fi
for param in "${CMD_POOL[@]}"
do
	echo "$LO_BIN $param" >> $TMP_JOB
	##Remove the absolute path
	pp=`echo $param | tr " " "_"|sed 's/\.*\///g'`
	##Name the actual log file
	log_file=${BIN_NAME}_P${pp}.log
	echo -e "cat ${LO_STATS} >> $log_file" >> $TMP_JOB
	if [ -n "$LO_LOG" ]
	then
		echo -e "cat ${LO_LOG} >> $log_file" >> $TMP_JOB
	fi
	echo "rm -f $LO_LOG $LO_STATS" >> $TMP_JOB
done
echo "tar czf ${DATE}_${BIN_NAME}-logs.tar.gz ${BIN_NAME}*.log" >> $TMP_JOB
echo "rm ${BIN_NAME}*.log" >> $TMP_JOB
#############################################################
### Constructing the actual shell script
############################################################
MACHINE_ERROR="Unrecognizable machine name, exit"
if [[ $LO_MACHINE == native ]]
then
	echo -e $NATIVE_PREFIX > $NATIVE_JOB
	cat $TMP_JOB >> $NATIVE_JOB
	chmod +x $NATIVE_JOB
	echo "Transfer control to the job script"
	./$NATIVE_JOB
elif [[ $LO_MACHINE == MN ]]
then
	echo -e $MN_PREFIX	> $MN_JOB
	cat $TMP_JOB >> $MN_JOB
	chmod +x $NATIVE_JOB
	echo "submitting job..."
	bsub < $MN_JOB
else
	echo $MACHINE_ERROR
	exit
fi
rm -f $TMP_JOB

