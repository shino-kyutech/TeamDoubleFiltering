#!/bin/bash
if [ $# -ne 5 ] ; then
echo "Usage> expand_filenames fc6_ 0 19 .ftr 1"
echo "1. prefix      = prefix"
echo "2. from        = 0, 10, or 90"
echo "3. to          = 9, 19, or 99"
echo "4. sufix       = sufix"
echo "5. mode    = 0 -> range, 1 -> filenames"
exit 1
fi

prefix=$1; shift
from=$1; shift
to=$1; shift
sufix=$1; shift
mode=$1; shift

if [ ${mode} -eq 0 ] ; then

if [ ${from} -lt 0 ] ; then
    range=${to}
else
    if [ ${from} -lt 10 ] ; then
    range=0${from}
    else
    range=${from}
    fi

    if [ ${to} -lt 10 ] ; then
    range=${range}_0${to}
    else
    range=${range}_${to}
    fi
fi
echo ${range}
exit

fi

if [ ${from} -lt 0 ] ; then
    s="${prefix}${to}${sufix}"
else
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
exit

# ここから下は，以前の名残で，現在は未使用
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
fi

echo $s
