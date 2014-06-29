#!/bin/bash

# set up modules environment
source /etc/bash.bashrc
source /opt/modules/default/init/bash

# visit requires the gnu programming environment
# the following should take the steps required to set it up
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


# generate the image
module load visit
echo "======== running visit to generate image ========"
visit -nowin -cli -s ./phi.py &> log
if [ $? -ne 0 ]; then
    echo "==================="
    echo "ERROR: visit failed"
    cat log
    exit
fi
grep Saved log

# get output file name
fname=`awk '/VisIt: Message - Saved/ {print $5}' log`

# draw picture
# only do so if x is running!
if env | grep -q ^DISPLAY=
then
    module load image-magick
    echo "======== drawing $fname ========"
    display $fname
else
    echo "can't draw image: you need to have x-windows running"
fi
