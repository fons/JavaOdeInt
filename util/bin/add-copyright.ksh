#!/usr/bin/env ksh

shcop="/tmp/shell_cop.txt"

if [[ -e $shcop ]];then
    rm $shcop
fi

while read line
do
    echo "# $line" >> $shcop
done < `pwd`/../../copyright/copyright.txt


cfile_cop="/tmp/cfile_cop.txt"
if [[ -e $cfile_cop ]]; then
    rm $cfile_cop
fi

echo "/* " >> $cfile_cop
while read line
do
    echo "* $line" >> $cfile_cop
done < `pwd`/../../copyright/copyright.txt
echo "*/ " >> $cfile_cop

path=$1


for file in `find $path -name "*.h"`
do
    has_notice=`grep -i copyright $file`
    if [[ -z $has_notice ]]; then
        fn=`basename $file`
        cat $cfile_cop > /tmp/$fn
        cat $file >> /tmp/$fn
        cp /tmp/$fn $file
    fi
done
