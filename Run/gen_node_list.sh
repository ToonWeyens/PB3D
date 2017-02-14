#!/bin/bash

usage() {
  echo "usage: $0 [-m <max_nodes>] <queue> [disabled_node]+"
  echo ""
  echo "This program reads the output of the pbsnodes command and creates a list of all free nodes"
  echo ""
  echo "Optional arguments"
  echo "  -m max_nodes output at most max_nodes in the list"
  echo "  -[h?]        show this help"
}
function join { local IFS="$1"; shift; echo "$*"; }

while getopts m: o
do	case "$o" in
        m)      num_nodes="$OPTARG";;
	[?h])	usage; exit 1;;
	esac
done
shift $((OPTIND-1))

if [ "$#" -lt 1 ]; then
  echo "Too few arguments given"
  echo ""
  usage; exit 1
fi
queue=$1
shift 1

# Keep a list of eligible nodes
nodes=()

# Works for torque, adjust syntax for other systems?
IFS=$'%';for node in `pbsnodes | sed 's/^$/%/g'`; do # Loop over each node
  # Test if the node has the right queue, is free and has no jobs
  if ! grep --quiet "properties = ${queue}$" <<< $node; then continue; fi
  if ! grep --quiet "state = free$" <<< $node; then continue; fi
  if grep --quiet "jobs = " <<< $node; then continue; fi

  node_name=`grep -o "^[a-zA-Z0-9]*$" <<< $node`
  # Exclude it if it matches the list of disabled nodes
  if grep --quiet "\b$node_name\b" <<< "$@"; then continue; fi

  # Get number of processors on this node
  [[ $node =~ "np = ([0-9]+)" ]]
  np="${BASH_REMATCH[1]}"

  nodes+=("$node_name:ppn=$np")
done

# Restore IFS
IFS=$' \t\n'
# truncate length of array to num_nodes if it exists
if [ ! -z "$num_nodes" ]; then
  nodes=( "${nodes[@]:0:num_nodes}" )
fi

# output result joined with +
join + ${nodes[@]}
