#!/usr/bin/bash


######## Change these #######
REMOTE="remote_name"
REMOTE_USER=$USER
REMOTE_DIR="/path/on/remote/server"
JOB_PREFIX="AMPA_fst"
MAXEPOCH=10

# Max number of connection failures before giving up and exiting
MAXFAILS=100
#############################


# wait times (seconds)
# time before retrying after connection failures
connection_retry=10
# If a newly submitted job is finished before this, assume things went badly and exit
new_epoch_fail=1800
# If the simulation is running, wait this long until rechecking
sim_running=600


function time_echo () {
    echo "[$(date +%T)]" "$@"
}

function time_echo_previous_line () {
    echo -e "\e[1A\e[K[$(date +%T)]" "$@"
}


function get_queue_njobs () {
    if sque_output=$(ssh ${REMOTE_USER}@${REMOTE} "squeue -u ${REMOTE_USER} -o '%A %j' "); then
        njobs=$(echo ${sque_output} | grep ${JOB_PREFIX} | wc -l) 
        return 0
    else
        rc=$?
        # The return code for a failed connection is 255. Any other nonzero code most likely means the command failed on the remote.
        if [ $rc -eq 255 ]; then 
            return 255
        else
            >&2 time_echo "ssh command failed with return code $rc. This most likely means the remote command failed instead of the connection. Stopping..."
            exit 1
        fi
    fi
}


function send_sbatch_job () {
    # the -t parameter is needed for module-command to be present.
    # Something about login and non-login shells.       I think...
    if ssh -t ${REMOTE_USER}@${REMOTE} "cd ${REMOTE_DIR}/epoch$(printf '%02d' $i ); sbatch sbatch_launch.sh" ; then
        return 0
    else
        rc=$?
        # The return code for a failed connection is 255. Any other nonzero code most likely means the command failed on the remote.
        if [ $rc -eq 255 ]; then 
            return 255
        else
            >&2 time_echo "ssh command failed with return code $rc. This most likely means the remote command failed instead of the connection. Stopping..."
            exit 1
        fi
    fi
}



function fst_rsync_cmd () {
    if fst $1 ; then
        return 0
    else
        rc=$?
        # The return code for a failed connection is 255. Any other nonzero code most likely means the command failed on the remote.
        if [ $rc -eq 255 ]; then 
            return 255
        else
            >&2 time_echo "fst $1 command failed with return code $rc. This most likely the problem is not the connection. Stopping..."
            exit 1
        fi
    fi
}

nfails=0
nsleep=0

while true ; do
    # if initial epoch does not exist, initialise it
    if ! [ -d epoch01 ] ; then
        nsleep=0

        time_echo "Running 'fst init'"
        # Make initial epoch
        if ! fst init; then
            >&2 time_echo "'fst init' failed, stopping..."
            exit 1
        fi

    else

        # retry connecting until works
        while ! get_queue_njobs ; do
            nsleep=0
            # keep track of number of failed attempts
            nfails=$((nfails+1))
            if ! [ -z $nfails ] && [[ $nfails -gt $MAXFAILS ]]; then
                >&2 time_echo "Connecting failed too many times, exiting"
                exit 1
            fi
            time_echo "getting squeue failed, retrying in $connection_retry s"
            sleep $connection_retry
        done

        # njobs is a global var set by get_queue_njobs
        if ! [[ $njobs -eq 0 ]]; then
            if [[ $nsleep -eq 0 ]]; then
                time_echo "Found job matching prefix, sleeping for $sim_running s"
            elif [[ $nsleep -eq 1 ]]; then
                time_echo "Found job matching prefix, sleeping again."
            elif [[ $nsleep -eq 2 ]]; then
                time_echo_previous_line "Found job matching prefix, sleeping for the 3rd time."
            else
                time_echo_previous_line  "Found job matching prefix, sleeping for the "$((nsleep+1))"th time."
            fi
            sleep $sim_running
            nsleep=$((nsleep+1))
            continue
        fi

        nsleep=0

        time_echo "No matching jobs found, making a new epoch"

        # Download existing data
        while ! fst_rsync_cmd pull ; do
            nfails=$((nfails+1))
            if ! [ -z $nfails ] && [[ $nfails -gt $MAXFAILS ]]; then
                >&2 time_echo "Connecting failed too many times, exiting"
                exit 1
            fi
            time_echo "'fst pull' failed, retrying in $connection_retry s"
            sleep $connection_retry
        done

        # Check if enough epochs are done
        if ! [ -z $MAXEPOCH ] && [ -d $(printf "epoch%02d\n" $MAXEPOCH) ]; then
            time_echo "MAXEPOCH folder ($(printf "epoch%02d\n" $MAXEPOCH)) found and simulation no longer running, I reckon we're done here"
            time_echo "Making final plots"
            fst choose --choose-only
            time_echo "All done, see you again soon"
            exit
        fi

        time_echo "Running 'fst choose'"
        # Make new epoch
        if ! fst choose; then
            >&2 time_echo "'fst choose' failed, stopping..."
            exit 1
        fi

    fi

    # The rest is done for both initialisation and main loop

    time_echo "Pushing changes"
    while ! fst_rsync_cmd push ; do
        nfails=$((nfails+1))
        if ! [ -z $nfails ] && [[ $nfails -gt $MAXFAILS ]]; then
            time_echo "Connecting failed too many times, exiting"
            exit 1
        fi
        time_echo "'fst push' failed, retrying in $connection_retry s"
        sleep $connection_retry
    done

    # Find the highest epoch number found and assume it is the new epoch to run
    i=0
    while [ -d epoch$(printf "%02d" $((i+1)) ) ]; do
        i=$((i+1))
    done

    time_echo "The newest epoch is $i"

    # Try putting in queue until it works
    while ! send_sbatch_job ; do
        nfails=$((nfails+1))
        if ! [ -z $nfails ] && [[ $nfails -gt $MAXFAILS ]]; then
            >&2 time_echo "Connecting failed too many times, exiting"
            exit 1
        fi
    done

    time_echo "Letting job run for $new_epoch_fail s"

    # Sleeping for 1 min just to make sure the job is firmly in the queue before checking
    sleep $new_epoch_fail

    time_echo "Checking if job is still in queue"
    # If the job is already done, something must be wrong.
    if ! get_queue_njobs; then
        time_echo "Failed to check if squeue output, we'll just continue looping and hope everything is fine"
        continue
    fi

    if [[ $njobs -eq 0 ]]; then
        >&2 time_echo "Job finished in $new_epoch_fail s, something is wrong"
        exit 1
    fi

done
