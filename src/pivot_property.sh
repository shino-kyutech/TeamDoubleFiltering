#!/bin/bash
if [ ! -e $2 ] ; then
echo $2 does not exist!
exit 0
fi

if [ $1 == "-w" ] ; then
head -n 2 $2 | tail -n 1 | tr "," " " | wc -w
exit 0
fi

if [ $1 == "-d" ] ; then
head -n 3 $2 | tail -n 1 | tr "," " " | wc -w
exit 0
fi

if [ $1 == "-p" ] ; then
head -n 1 $2 | tr -d " " | tr "," "\n" | head -n 1
exit 0
fi

if [ $1 == "-n" ] ; then
head -n 1 $2 | tr -d " " | tr "," "\n" | head -n 2 | tail -n 1
exit 0
fi

if [ $1 == "-r" ] ; then
r=$(head -n 1 $2 |  tr "," " " | wc -w)
if [ $r -gt 2 ] ; then
echo "ROTATION"
else
echo "NO_ROTATION"
exit 0
fi
fi