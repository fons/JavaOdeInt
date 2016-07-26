#!/bin/bash
#

action=$1

echo "action $action"

BASE=`dirname $0`
TARGET="$BASE/../../../target"
RESOURCE="$BASE/../main/resources"


CODEPACK_SRC="$BASE/../../../codeint/codepack/src"
CODEPACK_INCL="$BASE/../../../codeint/codepack/include"
CODEPACK_TARGET_INCL="$BASE/../../src/main/native/codeint/codepack/include"

CVODE_SRC="$BASE/../../../codeint/cdvode/src"
CVODE_INCL="$BASE/../../../codeint/cdvode/include"
CVODE_TARGET_INCL="$BASE/../../src/main/native/codeint/cdvode/include"

realclean()
{
    cmd="rm -rf /tmp/resource_*"
    echo $cmd
    eval $cmd
    
    cmd="rm -rf /tmp/codepack_*"
    echo $cmd
    eval $cmd
    
    cmd="rm -rf /tmp/cdvode_*"
    echo $cmd
    eval $cmd

}

clean()
{
    TS=`date +20%y%m%d_%H%M%S`

    if [[ -d $RESOURCE ]]; then
        cmd="mv $RESOURCE /tmp/resource_${TS}"
        echo $cmd
        eval $cmd
    fi

    if [[ -d $CODEPACK_TARGET_INCL ]]; then
        cmd="mkdir -p /tmp/codepack_${TS}/include"
        echo $cmd
        eval $cmd

        cmd="cp $CODEPACK_TARGET_INCL/* /tmp/codepack_${TS}/include"
        echo $cmd
        eval $cmd

        cmd="rm $CODEPACK_TARGET_INCL/* "
        echo $cmd
        eval $cmd
        
    fi

    if [[ -d $CVODE_TARGET_INCL ]]; then
        cmd="mkdir -p /tmp/cdvode_${TS}/include"
        echo $cmd
        eval $cmd

        cmd="cp $CVODE_TARGET_INCL/* /tmp/cdvode_${TS}/include"
        echo $cmd
        eval $cmd

        cmd="rm -f $CVODE_TARGET_INCL/* "
        echo $cmd
        eval $cmd
        
    fi


}

copy()
{
    cmd="cp -R $TARGET/ $RESOURCE/"
    echo $cmd
    eval $cmd

    cmd="cp $CODEPACK_INCL/*.h $CODEPACK_TARGET_INCL/."
    echo $cmd
    eval $cmd

    cmd="cp $CVODE_INCL/*.h $CVODE_TARGET_INCL/."
    echo $cmd
    eval $cmd
}

make_all()
{
    cmd="cd $CODEPACK_SRC && make && cd -"
    echo $cmd
    eval $cmd

    cmd="cd $CVODE_SRC && make && cd -"
    echo $cmd
    eval $cmd
    
}

if [[ $action = "realclean" ]]; then
    clean
    realclean
    exit 0
fi
if [[ $action = "clean" ]]; then
    clean
    exit 0
fi

make_all
copy


