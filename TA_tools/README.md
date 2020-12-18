# TA_tools
Tools for downloading student assignments in bash 


You will first need to set up SSH Key with this github account (use your own email):
https://jdblischak.github.io/2014-09-18-chicago/novice/git/05-sshkeys.html#:~:text=Login%20to%20github.com%20and,hit%20Add%20key%20to%20save



In the example below, `allHTMLURLs.A1.txt` is a file with one line per assignment, where the line is the full path to the `blob` HTML submission of the student's assignment
An example `allHTMLURLs.A1.txt` is included.

Example usage:
```
TA_tools/preprocess.sh allHTMLURLs.A1.txt allA1Paths.txt

mkdir A1_all; cd A1_all;
#../allA1Paths.txt is the file that contains the information needed to download, rename, and move the assignments to the desired location. This was created above.
#2020-10-23 is the date before which assignments will be considered (e.g. only updates before this date will be considered)
#../A1_compiled is where the files will end up
../TA_tools/downloadAllAssignments.sh ../allA1Paths.txt 2020-10-23 ../A1_compiled/
```
Now all the html assignments should be in `../A1_compiled`

If some people have submitted the assignment as multiple HTML files, you can use `mergeAssignmentHTMLs.sh`. An example of how to use it is included here: `exampleMergeHTMLs.sh`. Note that the order here matters. Make sure the components are compiled in the correct order. 
