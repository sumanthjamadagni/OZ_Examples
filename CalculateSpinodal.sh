
for T in 0.80 0.90 1.0 1.10 1.20 1.30 
  do
  echo T = $T
  python IsothermCalculation.py $T T$T-Vap.dat vap
done 

for T in 0.80 0.90 1.0 1.10 1.20 1.30 
  do
  echo T = $T
  python IsothermCalculation.py $T T$T-Liq.dat liq
done 

