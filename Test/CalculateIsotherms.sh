echo "Usage: CalculateSpinodal.sh m n eps_crit" 
m=$1 #50
n=$2 #12
eps_crit=$3 #1.05  #from the Tc-estimate.dat file in the Identify-Tc/ folder

branch_vap="vap"
#branch_liq="liq"

mkdir m$m-n$n
cd m$m-n$n
python ../IsothermCalculation.py $m $n $eps_crit $branch_vap
#python ../IsothermCalculation.py $m $n $eps_crit $branch_liq 
cd ../

