#!/bin/bash

echo for each line in $1 \(user	repo	assignmentFile\), cloning all repos for $2 date \(yyyy-mm-dd\), and putting assignments \(renamed to user name.ext\) in $3

startDir=`pwd`
mkdir -p $3
while read line; do
	splitID=($line) 
	curUN=${splitID[0]}
	curRepo=${splitID[1]}
	curBranch=${splitID[2]}
	curPath=`echo ${splitID[3]} | sed 's/%20/ /g'`
	echo user: $curUN
	echo repo: $curRepo
	echo branch: $curBranch
	echo file: $curPath
	#make dir for user
	mkdir $curUN
	cd $curUN

	#clone repo
	git clone git@github.com:$curUN/$curRepo.git
	#set date for repo
	cd $curRepo
	tempVar=`git rev-parse --until=$2`
	echo git rev-parse date: $tempVar
	echo git rev-list -1 $tempVar $curBranch
	curCommit=`git rev-list -1 $tempVar $curBranch`
	echo last commit before $2: $curCommit
	git reset --hard $curCommit

	#move file to desired folder
	newFN=$3/$curUN"."`echo $curPath | sed 's/.*\.\([^.]\+\)$/\1/g'`

	echo cp "$curPath" "$startDir/$newFN"
	cp "$curPath" "$startDir/$newFN"
	cd $startDir
done < $1
