#!/bin/bash
echo processing URLs in $1 and storing as "user\trepo\tbranch\tpath" in $2
cat $1 | sed 's/.*github.com\/\([^/]\+\)\/\([^/]\+\)\/blob\/\(main\|master\)\//\1\t\2\t\3\t/g' >  $2
