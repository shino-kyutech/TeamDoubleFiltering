#!/bin/bash
if [ $# -ne 3 ] ; then
echo "Usage> expand_filenames fc6_ .ftr"
echo "1. prefix      = prefix"
echo "2. range       = 00_01, 00_09, or 00_99"
echo "3. sufix       = sufix"
exit 1
fi

prefix=$1; shift
range=$1; shift
sufix=$1; shift

if [ ${range} == "00" ] ; then
s=${prefix}${range}${sufix}
elif [ ${range} == "00_01" ] ; then
s=""
for f in {00..01} ; do
s="$s ${prefix}$f${sufix}"
done
elif [ ${range} == "00_09" ] ; then
s=""
for f in {00..09} ; do
s="$s ${prefix}$f${sufix}"
done
elif [ ${range} == "10_19" ] ; then
s=""
for f in {10..19} ; do
s="$s ${prefix}$f${sufix}"
done
elif [ ${range} == "20_29" ] ; then
s=""
for f in {20..29} ; do
s="$s ${prefix}$f${sufix}"
done
elif [ ${range} == "30_39" ] ; then
s=""
for f in {30..39} ; do
s="$s ${prefix}$f${sufix}"
done
elif [ ${range} == "40_49" ] ; then
s=""
for f in {40..49} ; do
s="$s ${prefix}$f${sufix}"
done
elif [ ${range} == "50_59" ] ; then
s=""
for f in {50..59} ; do
s="$s ${prefix}$f${sufix}"
done
elif [ ${range} == "60_69" ] ; then
s=""
for f in {60..69} ; do
s="$s ${prefix}$f${sufix}"
done
elif [ ${range} == "70_79" ] ; then
s=""
for f in {70..79} ; do
s="$s ${prefix}$f${sufix}"
done
elif [ ${range} == "80_89" ] ; then
s=""
for f in {80..89} ; do
s="$s ${prefix}$f${sufix}"
done
elif [ ${range} == "90_99" ] ; then
s=""
for f in {90..99} ; do
s="$s ${prefix}$f${sufix}"
done
elif [ ${range} == "00_99" ] ; then
s=""
for f in {00..99} ; do
s="$s ${prefix}$f${sufix}"
done
elif [ ${range} == "00_96" ] ; then
s=""
for f in {00..96} ; do
s="$s ${prefix}$f${sufix}"
done
else
from=$(echo ${range} | tr _ '\n' | head -n 1)
to=$(echo ${range} | tr _ '\n' | tail -n 1)
#echo from = $from
#echo to = $to
    s=""
    for f in $(seq ${from} ${to}) ; do
    if [ $f -lt 10 ] ; then
    s="$s ${prefix}0$f${sufix}"
    else
    s="$s ${prefix}$f${sufix}"
    fi
    done
fi

echo $s
