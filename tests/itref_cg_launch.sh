#!/bin/bash
DIM=32000
BS=500
RHS=1
RES=1E-16
REP=1
WTIME="00:30"
RESOLUTION=$DIM

ITER_IR=40
ITER_CG=40
WARMUP=(0 5 10 20 30)
#DISCARD=("dynamic")
DISCARD=(0 1 3 5 10 20 30 40 50 "dynamic")
VER=(CGPROF)
#VER=(CGPROF CGMOD1)

ZBIN=/home/bsc28/bsc28687/apps/linompss_test/build/bin/itref

A_DIR=/gpfs/scratch/bsc28/bsc28687/ibm_async/fac1/A.dat
RHS_DIR=/gpfs/scratch/bsc28/bsc28687/ibm_async/fac1/RHS.dat
CFG=./linompss.cfg

#########################
#####    mkdirs
########################
BASE=`pwd`
for ver in ${VER[@]}
do
	for warmup in ${WARMUP[@]}
	do
		for discard in ${DISCARD[@]}
		do
			DIR=${ver}_w${warmup}_r${discard}
			if [ -e $DIR ]
			then
				rm -Rfv $DIR
			fi
			mkdir $DIR
			cd $DIR
			if [ "${discard}" == "dynamic" ]
			then
				BIN=${ZBIN}
				ASYNC_CMD="$BIN $DIM $BS $RHS $RHS $ITER_IR $ITER_CG $RES 0 1 $RESOLUTION $REP $warmup 1989 $RHS_DIR $A_DIR"
			elif [ "${discard}" == "0" ]
			then
				BIN=${ZBIN}
				ASYNC_CMD="$BIN $DIM $BS $RHS $RHS $ITER_IR $ITER_CG $RES 0 1 $RESOLUTION $REP $warmup 0 $RHS_DIR $A_DIR"
			else
				BIN=${ZBIN}
				ASYNC_CMD="$BIN $DIM $BS $RHS $RHS $ITER_IR $ITER_CG $RES 0 1 $RESOLUTION $REP $warmup $discard $RHS_DIR $A_DIR"
			fi
			#########################
			#####    setup script
			########################
			MN="#! /bin/bash\n"
			MN+="#BSUB -n 1\n"
			MN+="#BSUB -oo ${ver}_%J.out\n"
			MN+="#BSUB -eo ${ver}_%J.err\n"
			MN+="#BSUB -J ${ver}_${warmup}_${discard}\n"
			MN+="#BSUB -W ${WTIME}\n"
			MN+="#BSUB -x\n"
			MN+="\n\n"
			MN+="export LINOMPSS_CFG=$CFG\n"
			MN+="export NX_ARGS=\"--summary --schedule=bf --schedule-priority=yes\"\n"
			MN+="#./itref n bm bn s refit [it] [prec] [lookahead] [async] [profile] [rep] [warmup] [RHS] [A]\n"
			MN+="sed -i \"s/\\(solver=\\).*/\1${ver}/\" $CFG\n"
			MN+="${ASYNC_CMD} >> ${ver}_w${warmup}_r${discard}.log\n"
			MN+="mv itref_0.log ${ver}_w${warmup}_r${discard}_res.log\n"
			########################
			#####  Finish
			########################
			########################
			####  setup cfg
			########################
			MN1="ITREF: \n"
			MN1+="{\n"
			MN1+="	A="covar"\n"
			MN1+="	Apar=[2.0,2.0]\n"
			MN1+="	solver=CGMOD1\n"
			MN1+="};\n\n"
			MN1+="CG: \n"
			MN1+="{\n"
			MN1+="	A="covar"\n"
			MN1+="	Apar=[2.0,2.0]\n"
			MN1+="	dump=("A", "rhs")\n"
			MN1+="	log="no"\n"
			MN1+="};\n"
			########################
			#####  Finish
			########################
			echo -e $MN > linompss.sh
			chmod +x *.sh
			echo -e $MN1 > linompss.cfg
			bsub < linompss.sh
			cd $BASE
		done
	done
done

#########################
##### BASELINE
#########################
BASELINE="baseline"
mkdir $BASELINE
cd $BASELINE
BIN=${ZBIN}
ASYNC_CMD="$BIN $DIM $BS $RHS $RHS $ITER_IR $ITER_CG $RES 0 0 $RESOLUTION $REP 1 0 $RHS_DIR $A_DIR"
#########################
#####    setup script
########################
MN="#! /bin/bash\n"
MN+="#BSUB -n 1\n"
MN+="#BSUB -oo ${ver}_%J.out\n"
MN+="#BSUB -eo ${ver}_%J.err\n"
MN+="#BSUB -J ${ver}_baseline\n"
MN+="#BSUB -W ${WTIME}\n"
MN+="#BSUB -x\n"
MN+="\n\n"
MN+="export LINOMPSS_CFG=$CFG\n"
MN+="export NX_ARGS=\"--summary --schedule=bf --schedule-priority=yes\"\n"
MN+="#./itref n bm bn s refit [it] [prec] [lookahead] [async] [profile] [rep] [warmup] [RHS] [A]\n"
MN+="sed -i \"s/\\(solver=\\).*/\1${ver}/\" $CFG\n"
MN+="${ASYNC_CMD} >> ${ver}_baseline.log\n"
MN+="mv itref_0.log ${ver}_baseline_res.log\n"
########################
#####  Finish
########################
########################
####  setup cfg
########################
MN1="ITREF: \n"
MN1+="{\n"
MN1+="	A="covar"\n"
MN1+="	Apar=[2.0,2.0]\n"
MN1+="	solver=CGMOD1\n"
MN1+="};\n\n"
MN1+="CG: \n"
MN1+="{\n"
MN1+="	A="covar"\n"
MN1+="	Apar=[2.0,2.0]\n"
MN1+="	dump=("A", "rhs")\n"
MN1+="	log="no"\n"
MN1+="};\n"
########################
#####  Finish
########################
echo -e $MN > linompss.sh
chmod +x *.sh
echo -e $MN1 > linompss.cfg
bsub < linompss.sh
cd $BASE
