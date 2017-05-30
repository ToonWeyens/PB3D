#!/bin/bash
################################################################################
# Executes command with a timeout
# Params:
#   $1 timeout in seconds
#   $2 command
# Returns 1 if timed out 0 otherwise
timeout() {

    time=$1

    # start the command in a subshell to avoid problem with pipes
    # (spawn accepts one command)
    command="/bin/sh -c \"$2\""

    expect -c "set echo \"-noecho\"; set timeout $time; spawn -noecho $command; expect timeout { exit 1 } eof { exit 0 }"

    if [ $? = 1 ] ; then
        echo "Timeout after ${time} seconds"
    fi

}

IFS=$'\n' inputs=($(pbsnodes | grep '^c0' | grep -E -v 'c03|c05s10'))
n_inputs=${#inputs[@]}
broken_nodes_file="broken_nodes.txt"
rm -f $broken_nodes_file
for (( i=1; i<=$n_inputs; i++ )); do
    input=${inputs[$i-1]}
    echo "node $input in queue '"$(pbsnodes $input | grep 'properties' | cut -d = -f 2 | sed -e 's/^ *//g' -e 's/ *$//g')"'"
    disk_info=$(timeout 5 "ssh $input 'hostname && df -h | grep /dev/sda1'")
    if [[ ${disk_info} != *"Timeout"* && "${disk_info}" == *\%* ]]
    then
        echo $disk_info
    else
        echo " is BROKEN"
        echo $input >> $broken_nodes_file
    fi
done
[ -s $broken_nodes_file ] && printf "\n>>> Have a look at $broken_nodes_file for a list of broken nodes:\n"
[ -s $broken_nodes_file ] && cat $broken_nodes_file
