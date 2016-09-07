#!/usr/bin/env ksh

make_let()
{
    
    package=$1
    action=$2
    PACK_SRC="./${package}/src"
    if [[ ! -e $PACK_SRC ]];then
        return 0
    fi 

    if [[ ! -e $PACK_SRC/Makefile ]];then
        return 0
    fi 

    cmd="cd $PACK_SRC && make $action && cd -"
    
    echo $cmd
    
    eval $cmd
    
    if  [[ $? -ne 0 ]]; then
        
        echo "C TARGET MAKE FAILURE make failed for package $package"
        
        exit 1
        
    fi

    
}

action=$1
for f in `ls`
do
    make_let $f $action
done

