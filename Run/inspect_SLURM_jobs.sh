#!/bin/bash
# Finds the directories of all the jobs currently running
# In every directory, it gives the tail of PB3D.out
# The user can decide to cancel the job
user=toon
IFS=$'\n'
jobs=($(squeue --format "%Z %A" -t R -u $user))

for (( i=0; i<${#jobs[@]}; i++ )); do
    # first element is directory, rest is job id
    job_dir=$(echo ${jobs[i]} | cut -d " " -f 1 )
    job_id=$(echo ${jobs[i]} | cut --complement -d " " -f 1 )
    
    if [[ "$job_dir" != "WORK_DIR" ]]; then
        echo
        echo found job $job_id in directory $job_dir with tail ;
        echo 
        tail $job_dir/PB3D.out
        read -p "Do you want to cancel this job (y/N)?" choice
        case "$choice" in 
            y|Y ) 
                read -p ">>> Are you completely sure? (y/N)?" choice
                case "$choice" in 
                    y|Y ) 
                        scancel $job_id
                        echo 
                        echo ">>> canceled job" $job_id "in directory" $job_dir
                    ;;
                    * ) 
                    ;;
                esac
            ;;
            * ) 
            ;;
        esac
    fi
done

echo
echo "queue:"
echo
squeue
