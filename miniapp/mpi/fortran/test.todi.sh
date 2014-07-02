xdim=256
echo "======================================="
echo "every second core"
echo "======================================="
echo aprun -n 1 ./main $xdim $xdim 20 0.005
aprun -n 1 ./main $xdim $xdim 20 0.005 | grep "simulation took"
echo aprun -n 2 ./main $xdim $xdim 20 0.005
aprun -n 2 ./main $xdim $xdim 20 0.005 | grep "simulation took"
echo aprun -n 3 ./main $xdim $xdim 20 0.005
aprun -n 3 ./main $xdim $xdim 20 0.005 | grep "simulation took"
echo aprun -n 4 ./main $xdim $xdim 20 0.005
aprun -n 4 ./main $xdim $xdim 20 0.005 | grep "simulation took"
echo aprun -n 8 ./main $xdim $xdim 20 0.005
aprun -n 8 ./main $xdim $xdim 20 0.005 | grep "simulation took"

