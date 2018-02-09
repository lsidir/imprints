# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

#!/bin/sh

for db in airtraffic bestdr7 cnet  nozhup sf100 sinterklaas
do
	echo "$db"
	cd /net/milan/export/scratch2/lsidir/finger_dbfarm/$db
	sh job$db > ${db}all.res
	echo "done $db"
done
