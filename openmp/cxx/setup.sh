if [ "$PE_ENV" != 'GNU' ]; then
    echo ">>>>>>>>>>>>>>> configuring the GNU programming environment"
    pe=`echo $PE_ENV | awk '{print tolower($0)}'`
    if [ "$PE_ENV" == '' ]; then
        echo ">>>>>>>>>>>>>>> no PE loaded: loading PrgEnv-gnu directly"
        module load PrgEnv-gnu
    else
        echo ">>>>>>>>>>>>>>> swapping out PrgEnv-${pe} for PrgEnv-gnu"
        module swap PrgEnv-${pe} PrgEnv-gnu
    fi
fi
