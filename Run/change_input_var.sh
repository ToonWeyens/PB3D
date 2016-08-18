#!/bin/bash
# to change an input variable of a PB3D or POST input file.
if [ "$#" -ne 3 ]; then
    echo -e "\nUsage:\n$0 INPUT_FILE VAR_NAME VAR_VAL \n" 
    exit 1
fi
# find the line of var
var_line=$(grep -nr -m 1 $2 $1)
echo var_line=$var_line
# find the line number of var
var_line_nr=$(echo $var_line | cut -d : -f 1)
echo var_line_nr=$var_line_nr
# replace whole line by new line
sed -i "${var_line_nr}s/.*/    $2 = $3/" "$1"

# Note: in an old version, only integer variables were treated and they were replaced inline, using
# # find the old value of var
# var_old=$(echo $var_line | cut -d : -f 2 | tr -dc '0-9')
# echo var_old=$var_old
# # replace old value of var by new value
# sed -i "${var_line_nr}s/${var_old}/${3}/" "$1"
