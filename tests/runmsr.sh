#!/bin/sh
#runs ONCVPSP with the command-line argument <prefix> and the graphics
#which review the results
#uses the scalar-relativistic all-electron atom calculation

PREFIX=/home/drh/metapsp-1.0.1

INFILE=$1.dat

OUTFILE=$1_m.out

GNUFILE=$1_m.scr

PLOTFILE=$1_m.plot

TEMP=$$.tmp

echo 'SR' >SR

$PREFIX/src/oncvpspm.x <$INFILE >$OUTFILE  #Edit if your executable is
                                            #in another directory

rm SR

awk 'BEGIN{out=0};/GNUSCRIPT/{out=0}; {if(out == 1) {print}};\
	/DATA FOR PLOTTING/{out=1}' $OUTFILE >$PLOTFILE

awk 'BEGIN{out=0};/END_GNU/{out=0}; {if(out == 1) {print}};\
	/GNUSCRIPT/{out=1}' $OUTFILE >$TEMP

#sed -e 1,1000s/t1/$PLOTFILE/ $TEMP >$GNUFILE

sed -e s/t1/$PLOTFILE/ $TEMP | sed -e s/t2/$1/ >$GNUFILE

if [ "$2" != "-np" ]
then
	gnuplot $GNUFILE
fi

 rm  $GNUFILE $TEMP $PLOTFILE
