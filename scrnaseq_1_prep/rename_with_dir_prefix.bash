#!/bin/bash

#space separated list of dirs
pathname=("/home/rhooads/sra/aligned")
echo "dir list is: "${pathname[@]}

tmp=("/home/rhooads/sra/tmp")
echo "temp folder is: "$tmp

for d in ${pathname[@]}
do
cd $d
pwd


#only need if recursing through dirs (i.e. Sample_control-1)
for i in $(ls -d $d/*/)
do
cd ${i%%/}
echo  ${i%%/}
#SUFFIX=$(basename $i)
#echo $SUFFIX

for file in $(ls)
do
SUFFIX=$(basename $(ps -aef | grep -o `pwd`))
echo $SUFFIX
echo $file
mv $file ${SUFFIX}_$file
done

#end of initial 2 location loops
echo "leaving: "$(ps -aef | grep -o `pwd`)
cd ../
done

done


