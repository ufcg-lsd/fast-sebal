#!/bin/bash

export LC_ALL="C"
TIME_BETWEEN_COMMANDS=1
echo TIMESTAMP, PID, USER, SYSTEM, GUEST, CPU, COMMAND
while [ -e /proc/$1 ]; do
 /usr/bin/pidstat -p $1 -u $TIME_BETWEEN_COMMANDS 1 | grep -e '[0-9]\{2\}:[0-9]\{2\}:[0-9]\{2\}[ \t]*[0-9]' | awk -v date="$( date +"%s" )" '{print date", "$3", "$4", "$5", "$6", "$7", "$9; }' 2> /dev/null
done
