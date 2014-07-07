export MPICH_RDMA_ENABLED_CUDA=1
#unset MPICH_RDMA_ENABLED_CUDA

dim=512
nprocx=1
nprocy=1
nt=50
echo "============================================="
echo "weak scaling $dim*$dim for $nt steps"
echo "============================================="
for i in 1 2 3 4 5
do
    xdim=$[$nprocx * $dim]
    ydim=$[$nprocy * $dim]
    nproc=$[$nprocx * $nprocy]
    steps_per_second=`aprun -n $nproc  -N 1 ./main $xdim $ydim $nt 0.0025 | grep "per second" | awk '{printf("%9.1f", $5)}'`
    echo "CG iterations per second = $steps_per_second :: $nproc nodes ($xdim*$ydim)"
    if [ "$nprocx" -lt "$nprocy" ]
    then
        nprocx=$nprocy
    else
        nprocy=$[$nprocy * 2]
    fi
done

dim=1024
nt=100
echo "============================================="
echo "strong scaling $dim*$dim for $nt steps"
echo "============================================="
nproc=1
for i in 1 2 3 4 5
do
    steps_per_second=`aprun -n $nproc  -N 1 ./main $dim $dim $nt 0.0025 | grep "per second" | awk '{printf("%9.1f", $5)}'`
    echo "CG iterations per second = $steps_per_second :: $nproc nodes ($dim*$dim)"
    nproc=$[$nproc * 2];
done

