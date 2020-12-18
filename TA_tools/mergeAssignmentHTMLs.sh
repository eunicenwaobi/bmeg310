#!/bin/bash

first=1
for file in "$@"
do
	if [ $first -ne 0 ]
	then
		from="0"
	else
		from=`grep -A 1 -n "<body>" "$file"  | head -n 1 | sed 's/:.*$//g'`
	fi
	first=0;
	to=`grep -n  "</body>" "$file" | sed 's/:.*$//g'`
	#echo Getting file $file from $from to $to
	cat "$file" | awk "(NR > $from && NR < $to){print}"
done
cat "$file" | awk "(NR >= $to){print}"

