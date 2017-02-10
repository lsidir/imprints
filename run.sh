#!/bin/sh


for db in airtraffic bestdr7 cnet  nozhup sf100 sinterklaas
do
	echo "$db"
	cd /export/scratch2/lsidir/finger_dbfarm/$db
	sh job$db > ${db}all.res
	echo "done $db"
done
