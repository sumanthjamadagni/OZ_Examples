rm vap*_Pkappa.dat -rf  #These are created by the python script. Deleting old ones. 
rm liq*_Pkappa.dat -rf  #These are created by the python script. Deleting old ones. 

flist_liq=`(ls liq* )`
flist_vap=`(ls vap* )`

##Liquid side
rm -rf LiqSpinodal.dat 
for f in $flist_liq
  do
  echo file = $f
  rm -rf SpinodalPoint.dat
  python ../FitSpinodal_SingleIsotherm.py $f liq 

  if [ -e SpinodalPoint.dat ] 
      then
      cat SpinodalPoint.dat >> LiqSpinodal.dat
  fi

done 

#Vapor side
rm -rf VapSpinodal.dat 
for f in $flist_vap
  do
  echo file = $f
  rm -rf SpinodalPoint.dat
  python ../FitSpinodal_SingleIsotherm.py $f vap 

  if [ -e SpinodalPoint.dat ] 
      then
      cat SpinodalPoint.dat >> VapSpinodal.dat
  fi

done 
