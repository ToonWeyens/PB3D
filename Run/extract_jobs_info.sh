#!/bin/bash
# Extracts the Eigenvalues of jobs for an input array
# run this function as "[[ $# -ne A ]] && err_arg A"
err_arg() { 
    [[ $# -ge 1 ]] && nr_args_req=$1 || nr_args_req=some
    echo "ERROR: need $nr_args_req argument(s)"
    exit 1
}

[[ $# -ne 3 ]] && echo "usage: $0 [ARRAY INPUT] [RUN FOLDER] [OUTPUT FILE]"
[[ $# -ne 3 ]] && err_arg 3

array_input=$1
run_folder=${2%/}
output_file=$3
echo "Extracting Eigenvalue information from '$run_folder', which should correspond to '$array_input'"
echo
echo "The output is saved to '$output_file'"
echo

# write info in output file
echo "Eigenvalues from '$run_folder', corresponding to '$array_input':" > $output_file
grep "^[^#;]" "$array_input" | sed -e 's/^/: /' | cat -n >> $output_file
sed -i 's/^/# /' $output_file

# get all inputs of array input
IFS=$'\n' inputs=($(grep "^[^#;]" "$array_input"))
n_inputs=${#inputs[@]}

# loop over all inputs to count how many columns needed
nr_columns=0
for (( i=1; i<=$n_inputs; i++ )); do
    nr_commas_loc=$(echo ${inputs[$i-1]} | tr -cd , | wc -c)
    (( nr_commas_loc+1 > nr_columns )) && nr_columns=$((nr_commas_loc+1))
done
column_width_int=10
column_width_EV=26

# loop over all inputs
for (( i=1; i<=$n_inputs; i++ )); do
    input=${inputs[$i-1]}
    echo "input $i/$n_inputs:"
    echo "    input file: '$input'"
    
    # split input line in ',' elements (from http://stackoverflow.com/a/918931)
    IFS=',' read -ra input_elem <<< "$input"
    
    # set up first part of output string: input modifications
    output_string=$(printf "%-${column_width_int}s" $i)
    for (( j=1; j<=$nr_columns; j++ )); do
        # extract value, trim whitespaces
        val=$(echo ${input_elem[$j-1]} | cut --complement -d = -f 1 | sed -e 's/^ *//g' -e 's/ *$//g')
        
        if [[ -z "${val// }" ]]; then
            output_string=$output_string$(printf "%-${column_width_int}s" "..")
        else
            output_string=$output_string$(printf "%-${column_width_int}s" $val)
        fi
    done
    
    # check if error present
    if [[ -s $run_folder/$i/PB3D.err ]]; then
        echo "    ERROR FOUND:"
        sed -e 's/^/        /' $run_folder/$i/PB3D.err | cat
        echo
        sed -e 's/^/        /' $run_folder/$i/PB3D.out | tail
        output_string=$output_string$(printf "%-${column_width_EV}s" "ERR" "ERR")
    else
        # set up second part of output string: extracted Eigenvalue
        EV_files=$(ls $run_folder/$i/PB3D_out_EV_R*.txt 2> /dev/null)
        EV_file=$(echo $EV_files | rev | cut -d" " -f1 | rev)
        # check if something found: replace occurences of space by nothing
        if [[ -z "${EV_file// }" ]]; then
            echo "    no EV file found"
            for (( k=1; k<=${#EVs[@]}; k++ )); do
                IFS=' ' read -ra EV <<< "${EVs[$k-1]}"
                output_string=$output_string$(printf "%-${column_width_EV}s" ".." "..")
            done
        else
            echo "    EV file:    '$EV_file'"
            IFS=$'\n' EVs=($(grep "^[^#;]" "$EV_file"))
            echo "    writing ${#EVs[@]} Eigenvalue(s)"
            for (( k=1; k<=${#EVs[@]}; k++ )); do
                IFS=' ' read -ra EV <<< "${EVs[$k-1]}"
                output_string=$output_string$(printf "%-${column_width_EV}s" ${EV[1]} ${EV[2]})
            done
        fi
    fi
    echo $output_string >> $output_file
    echo
done
