
pydir=/hosts/sisyphus/home/jamadagni/Klein/OZSolver/SALR

m=100
n=50
ptype='v2'

mkdir m$m-n$n-$ptype
cd m$m-n$n-$ptype

for eps in 3.0
  do
  mkdir eps-$eps
  cd eps-$eps
  
  
  for cs in 0.50
    do
    mkdir cs-$cs
    cd cs-$cs
    python $pydir/IsothermCalculation.py -m $m -n $n --eps $eps --Cs $cs --ptype $ptype

    python $pydir/IsothermCalculation.py -m $m -n $n --eps $eps --Cs $cs --ptype $ptype --B2 

    cd ../    
  done
  cd ..
done



