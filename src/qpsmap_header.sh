#!/bin/bash
# QPSMAP ファイルのヘッダー情報を得る

sm=$1; shift
op=$1; shift

if [ $op == "DIM" ] ; then
skip=0
elif [ $op == "BIT" ] ; then
skip=4
elif [ $op == "NUM" ] ; then
skip=8
fi

dd if=${sm} bs=1 skip=${skip} count=4 2>/dev/null | hexdump -e '1/4 "%d\n"'
