#!/bin/bash
IFS=$'\n' inputs=($(pbsnodes | grep '^c0' | grep -E -v 'c03|c05s10'))
n_inputs=${#inputs[@]}
broken_nodes_file="broken_nodes.txt"
rm -f $broken_nodes_file
for (( i=1; i<=$n_inputs; i++ )); do
    input=${inputs[$i-1]}
    disk_info=$(ssh $input 'hostname && df -h | grep /dev/sda1')
    if [[ "$disk_info" =~ ( |\') ]]
    then
       echo $disk_info
    else
        echo "node $disk_info is BROKEN"
        echo $disk_info >> $broken_nodes_file
    fi
    echo "    in queue '"$(pbsnodes $input | grep 'properties' | cut -d = -f 2 | sed -e 's/^ *//g' -e 's/ *$//g')"'"
done
[ -s $broken_nodes_file ] && printf "\n>>> Have a look at $broken_nodes_file for a list of broken nodes\n"
