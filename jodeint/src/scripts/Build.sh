#!/bin/bash
#


realclean_let()
{
    package=$1
    
    cmd="rm -rf /tmp/${package}_*"
    echo $cmd
    eval $cmd
}


clean_let()
{
    package=$1
    PACK_TARGET_INCL="$BASE/../../src/main/native/codeint/${package}/"

    if [[ -d $PACK_TARGET_INCL ]]; then
        cmd="mkdir -p /tmp/${package}_${TS}/include"
        echo $cmd
        eval $cmd

        cmd="cp -R $PACK_TARGET_INCL/* /tmp/${package}_${TS}/include"
        echo $cmd
        eval $cmd

        cmd="rm $BASE/../../src/main/native/codeint/${package}/include/*"
        echo $cmd
        eval $cmd

        cmd="rm -rf $BASE/../../src/main/native/codeint/${package}/include/"
        echo $cmd
        eval $cmd

        cmd="rm -rf $BASE/../../src/main/native/codeint/${package}/"
        echo $cmd
        eval $cmd
    fi

}

clean()
{
    action=$1

    if [[ -z $action ]]; then
        return
    fi

    if [[ "$action" != "clean" && "$action" != "realclean" ]];then
        return
    fi

    if [[ -d $RESOURCE ]]; then
        cmd="mv $RESOURCE /tmp/resource_${TS}"
        echo $cmd
        eval $cmd
    fi

    #
    # this will clean the root package dir if it exists uder ./src/main/.
    # by default the sources are generated under ./target/generated-sources/
    #
    PACKAGE_E=$(echo $ROOTPACKAGE | sed 's/\./\//g' )
    if [[ -n $PACKAGE_E ]]; then
        PTJ="$BASE/../../src/main/java/$PACKAGE_E"
        echo $PACKAGE_E
        echo $PTJ
        if [[ -d $PTJ ]]; then
            cmd="rm -rf $PTJ"
            echo $cmd
            eval $cmd
        fi
    fi

    #
    # clean up the trash in /tmp/
    #

    if [[ $action = "realclean" ]]; then    
        cmd="rm -fr /tmp/resource_*"
        echo $cmd
        eval $cmd
        
        cmd="rm -f $JNAERATOR_TARGET/prev*"
        echo $cmd
        eval $cmd

        if [[ -d $JNAERATOR_TARGET ]]; then
            cmd="rmdir $JNAERATOR_TARGET"
            echo $cmd
            eval $cmd
        fi
        
    fi
    
    for f in `ls $BASE/../../../codeint/`
    do
        if [[ "$f" != "test" ]]; then
            realclean_let $f
            clean_let $f
        fi
    done


}

copy_let()
{
    package=$1
    PACK_INCL="$BASE/../../../codeint/${package}/include"
    PACK_TARGET_INCL="$BASE/../../src/main/native/codeint/${package}/include"

    if [[ ! -d $PACK_TARGET_INCL ]];
    then
        cmd="mkdir -p $PACK_TARGET_INCL"
        echo $cmd
        eval $cmd
    fi

    cmd="cp $PACK_INCL/*.h $PACK_TARGET_INCL/."
    echo $cmd
    eval $cmd

}

clean_f77_let()
{
    package=$1
    PACK_SRC="$BASE/../../../fodeint/${package}/src"
    cmd="cd $PACK_SRC && make clean && cd -"
    echo $cmd
    eval $cmd
    if  [[ $? -ne 0 ]]; then
        echo "F77 TARGET CLEAN FAILURE make failed for package $package"
        exit 1
    fi

}

make_f77_let()
{
    package=$1
    PACK_SRC="$BASE/../../../fodeint/${package}/src"
    cmd="cd $PACK_SRC && make && cd -"
    echo $cmd
    eval $cmd
    if  [[ $? -ne 0 ]]; then
        echo "F77 TARGET MAKE FAILURE make failed for package $package"
        exit 1
    fi

}

clean_test_let()
{
    echo `pwd`
    
    package=$1
    PACK_SRC="$BASE/../../../codeint/test/${package}/"
    cmd="cd $PACK_SRC && make clean && cd -"
    echo $cmd
    eval $cmd
    if  [[ $? -ne 0 ]]; then
        echo "TEST TARGET CLEAN FAILURE make failed for package $package"
        exit 1
    fi

}
make_test_let()
{
    echo `pwd`
    
    package=$1
    PACK_SRC="$BASE/../../../codeint/test/${package}/"
    cmd="cd $PACK_SRC && make && cd -"
    echo $cmd
    eval $cmd
    if  [[ $? -ne 0 ]]; then
        echo "TEST TARGET MAKE FAILURE make failed for package $package"
        exit 1
    fi

}

clean_let()
{
    package=$1
    PACK_SRC="$BASE/../../../codeint/${package}/src"
    cmd="cd $PACK_SRC && make clean && cd -"
    echo $cmd
    eval $cmd
    if  [[ $? -ne 0 ]]; then
        echo "C TARGET MAKE FAILURE make failed for package $package"
        exit 1
    fi
    
    
}


make_let()
{
    package=$1
    PACK_SRC="$BASE/../../../codeint/${package}/src"
    cmd="cd $PACK_SRC && make && cd -"
    echo $cmd
    eval $cmd
    if  [[ $? -ne 0 ]]; then
        echo "C TARGET MAKE FAILURE make failed for package $package"
        exit 1
    fi
}

jnaerator_let()
{

    package=$1
    JNAERATOR_CONFIG_PACKAGE="${NEW_JNAERATOR_CONFIG}_${package}"    

    echo "// -----  generated on ${TS} by `whoami` for $package  ---------------" > $JNAERATOR_CONFIG_PACKAGE
    echo "      " >> $JNAERATOR_CONFIG_PACKAGE
    echo "-rootPackage $ROOTPACKAGE" >> $JNAERATOR_CONFIG_PACKAGE
    echo "      " >> $JNAERATOR_CONFIG_PACKAGE
    echo "-library ${package}" >> $JNAERATOR_CONFIG_PACKAGE
    echo "-I\$(DIR)/../native/codeint/${package}/include" >> $JNAERATOR_CONFIG_PACKAGE
    echo "\$(DIR)/../native/codeint/${package}/include/${package}.h" >>  $JNAERATOR_CONFIG_PACKAGE
    echo " " >> $JNAERATOR_CONFIG_PACKAGE
    echo " " >> $JNAERATOR_CONFIG_PACKAGE


    if [[ ! -d $JNAERATOR_TARGET ]]; then
        mkdir -p $JNAERATOR_TARGET
    fi

    if [[ -f  $JNAERATOR/config.jnaerator_${package} ]]; then
        cmd="mv  $JNAERATOR/config.jnaerator_${package} $JNAERATOR_TARGET/prev_${TS}_config.jnaerator_${package}"
        echo $cmd
        eval $cmd
    fi    

    if [[ -f  $JNAERATOR_CONFIG_PACKAGE ]]; then
        cmd="mv  $JNAERATOR_CONFIG_PACKAGE $JNAERATOR/config.jnaerator_${package}"
        echo $cmd
        eval $cmd
    fi    


    cmd="ln -f $JNAERATOR/config.jnaerator_${package} $JNAERATOR/config.jnaerator" 
    echo $cmd
    eval $cmd

    cmd="mvn jnaerator:generate"
    echo $cmd
    eval $cmd
    if [[ $? -ne 0 ]]; then
        exit 1
    fi
}

clean_all_f77()
{

    for f in `ls $BASE/../../../fodeint/`
    do
        if [[ "$f" != "test" ]]; then
           clean_f77_let $f
        fi 
    done
}
make_all_f77()
{

    for f in `ls $BASE/../../../fodeint/`
    do
        if [[ "$f" != "test" ]]; then
            make_f77_let $f
        fi 
    done
}

clean_all_test()
{

    for f in `ls $BASE/../../../codeint/test`
    do
        clean_test_let $f
    done
}

make_all_test()
{

    for f in `ls $BASE/../../../codeint/test`
    do
        make_test_let $f
    done
}

clean_all_c()
{


    for f in `ls $BASE/../../../codeint/`
    do
        if [[ "$f" != "test" ]]; then
            clean_let $f
        fi 
    done
}

make_all_c()
{

    if [[ ! -e $JNAERATOR ]]; then
	cmd="mkdir -p $JNAERATOR"
	echo $cmd
	eval $cmd
    fi

    make_let "codeintdeps"
    for f in `ls $BASE/../../../codeint/`
    do
        if [[ "$f" != "test" ]]; then
            make_let $f
            if [[ "$f" != "codeintdeps" ]]; then
                echo " ==> ${JNAERATOR_CONFIG}_${f}"
                if [[ ! -e "${JNAERATOR_CONFIG}_${f}" ]]; then
                    copy_let $f
                    jnaerator_let $f
                fi
            fi
        fi 
    done
}


clean_all()
{
    clean_all_f77

    clean_all_test

    clean_all_c
    
    cmd="rm -rf $TARGET/"
    echo $cmd
    eval $cmd

    cmd="rm -rf $JNAERATOR/config.jnaerator*"
    echo $cmd
    eval $cmd

}

make_all()
{
    make_all_f77
    make_all_c
    make_all_test

    if [[  ! -e  "$RESOURCE" ]]; then
        cmd="mkdir -p $RESOURCE"
        echo $cmd
        eval $cmd
    fi
    
    if [[  ! -e  "$RESOURCE/lib/" ]]; then
        cmd="cp -R $TARGET/lib $RESOURCE/lib/"
        echo $cmd
        eval $cmd
    fi

}
action=$1

if [[ -n $action ]]; then 
    echo "action $action"
fi

BASE=`dirname $0`
TARGET="$BASE/../../../target"
RESOURCE="$BASE/../main/resources"
TS=`date +20%y%m%d_%H%M%S`

JNAERATOR="$BASE/../main/jnaerator/"
JNAERATOR_TARGET="$TARGET/jnaerator/"

NEW_JNAERATOR_CONFIG="$BASE/../main/jnaerator/new_jnaerator.config"
JNAERATOR_CONFIG="$BASE/../main/jnaerator/config.jnaerator"
ROOTPACKAGE="com.kabouterlabs.jodeint"

if [[ $action = "realclean" || $action = "clean" ]]; then
    clean $action
    clean_all
    exit 0
fi
 make_all



