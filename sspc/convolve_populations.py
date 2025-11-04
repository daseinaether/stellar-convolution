"""
Main file to handle the convolution of populations
"""

import copy
import json
import multiprocessing
import os
import pickle
import traceback
import warnings
from functools import partial

import h5py
import setproctitle

from syntheticstellarpopconvolve.convolve_on_the_fly import convolve_on_the_fly
from syntheticstellarpopconvolve.convolve_pre_calculated_data import (
    convolve_pre_calculated_data,
)
from syntheticstellarpopconvolve.general_functions import (
    JsonCustomEncoder,
    create_job_dict,
    create_time_bin_info_dict,
    generate_data_dict,
    generate_group_name,
    get_tmp_dir,
    has_unit,
    maybe_strip_scaled_dimensionless,
)


def handle_storing_convolution_results(  # DH0001
    config, grp, convolution_results, bin_center
):
    """
    Function to manage the storing of the convolution results
    """

    #########
    # Handle storing convolution results
    if isinstance(convolution_results, list):
        for convolution_result in convolution_results:

            # Create group
            current_time_bin_grp = grp.create_group(
                "convolution_results/{}/{}".format(
                    convolution_result["name"], str(bin_center)
                )
            )

            ############
            # handle storing entries and units
            config["logger"].debug(
                "Storing convolution results {} of bin-center {}".format(
                    convolution_result["name"], str(bin_center)
                )
            )

            #
            store_convolution_result_entries(
                config=config,
                current_time_bin_group=current_time_bin_grp,
                convolution_result=convolution_result,
            )

    else:
        # Create group
        current_time_bin_grp = grp.create_group(
            "convolution_results/{}".format(str(bin_center))
        )

        ############
        # handle storing entries and units
        config["logger"].debug(
            "Storing convolution results of bin-center {}".format(str(bin_center))
        )

        #
        store_convolution_result_entries(
            config=config,
            current_time_bin_group=current_time_bin_grp,
            convolution_result=convolution_results,
        )


def store_convolution_result_entries(  # DH0001
    config, current_time_bin_group, convolution_result
):
    """
    Function to handle storing an entry of the convolution_result
    """

    units_to_store = {}

    # loop over the entries
    for entry in convolution_result.keys():
        # skip name field
        if entry == "name":  # pass
            continue

        #
        config["logger"].error(f"Storing {entry}")

        # unpack data
        entry_data = convolution_result[entry]

        # handle fake units (dimensionless but scaled)
        entry_data = maybe_strip_scaled_dimensionless(entry_data)

        # handle storing data with units
        if has_unit(entry_data):
            current_time_bin_group.create_dataset(entry, data=entry_data.value)
            units_to_store[entry] = entry_data.unit
        # handle storing data without units
        else:
            current_time_bin_group.create_dataset(entry, data=entry_data)

    ###########
    # store units
    current_time_bin_group.attrs["units"] = json.dumps(
        units_to_store, cls=JsonCustomEncoder
    )


def pre_convolution(config, convolution_instruction, sfr_dict):  # DH0001
    """
    Function to handle things before a convolution.
    - Prepare output group structures in hdf5 as far as possible.
    - Stores SFR information in the output group.
    - Creates temporary directories.
    """

    ########
    # get groupname
    groupname, elements = generate_group_name(
        convolution_instruction=convolution_instruction, sfr_dict=sfr_dict
    )

    ########
    # Apply correct structure in hdf5 file
    with h5py.File(config["output_filename"], "a") as output_hdf5file:
        ########
        # Create output data group
        config["logger"].debug("Creating output data groups '{}'".format(groupname))

        #
        if "output_data" not in output_hdf5file.keys():
            output_hdf5file.create_group("output_data")

        # Create further structure of data group
        for depth in range(len(elements)):
            group = "/".join(elements[: depth + 1])
            if group not in output_hdf5file["output_data"]:
                output_hdf5file["output_data"].create_group(group)

        ########
        # store SFR dict
        if "name" in sfr_dict:
            group_ = "output_data/{}".format(sfr_dict["name"])
        else:
            group_ = "output_data"

        config["logger"].debug(
            "Storing SFR dict in attribute of group '{}'".format(group_)
        )

        #
        output_hdf5file[group_].attrs["SFR_info"] = json.dumps(
            sfr_dict, cls=JsonCustomEncoder
        )

    ########
    # create tmp dir
    tmp_dir = get_tmp_dir(
        config=config,
        convolution_instruction=convolution_instruction,
        sfr_dict=sfr_dict,
    )
    os.makedirs(tmp_dir, exist_ok=True)


def post_convolution(config, convolution_instruction, sfr_dict):  # DH0001
    """
    Function to handle post-convolution.

    Mostly stores tmp pickle files that contain the data
    """

    #################
    # Put pickle data in the hdf5 file
    tmp_dir = get_tmp_dir(
        config=config,
        convolution_instruction=convolution_instruction,
        sfr_dict=sfr_dict,
    )

    ########
    # Write results to output file
    if not config["write_to_hdf5"]:
        return

    # Get groupname
    groupname, _ = generate_group_name(
        convolution_instruction=convolution_instruction, sfr_dict=sfr_dict
    )
    full_groupname = "output_data/" + groupname

    with h5py.File(config["output_filename"], "a") as output_hdf5file:
        config["logger"].debug("Writing results to {}".format(full_groupname))

        # Readout group
        grp = output_hdf5file[full_groupname]

        ###########
        # loop over all pickle files that contain data
        content_dir = os.listdir(tmp_dir)

        sorted_content_dir = sorted(
            content_dir,
            key=lambda x: float(".".join(x.split(".")[:-1]).split(" ")[0]),
        )
        for pickle_file in sorted_content_dir:
            #########
            # check if file is actually pickle file
            if not pickle_file.endswith(".p"):
                continue

            #########
            # Load pickled data
            full_path = os.path.join(tmp_dir, pickle_file)
            with open(full_path, "rb") as picklefile:
                payload = pickle.load(picklefile)

            ##########
            # Unpack
            if "convolution_results" in payload.keys():
                convolution_results = payload["convolution_results"]
            else:  # TODO: do we want to raise an error or just continue?
                raise ValueError("No convolution result present in the data")

            #########
            # Handle storing convolution results
            handle_storing_convolution_results(
                config=config,
                grp=grp,
                convolution_results=convolution_results,
                bin_center=payload["bin_center"],
            )

            # remove the pickled file
            if config["remove_pickle_files"]:
                os.remove(full_path)


def handle_convolution_choice(  # DH0001
    config,
    job_dict,
    sfr_dict,
    convolution_instruction,
    data_dict,
    persistent_data=None,
    previous_convolution_results=None,
):
    """
    Function to handle the convolution choice
    """

    #
    time_bin_info_dict = job_dict["time_bin_info_dict"]

    ##################
    # Handle choice of convolution type.
    #   here we just handle whether the convolution uses pre-calculated data or not.
    #
    if convolution_instruction["convolution_type"] == "on-the-fly":
        warnings.warn("On-the-fly convolution is currently not fully tested")

        if convolution_instruction["contains_binned_data"]:
            raise ValueError(
                "Convolving binned data with on-the-fly convolution is currently not supported."
            )

        ##########
        #
        convolution_results = convolve_on_the_fly(
            config=config,
            sfr_dict=sfr_dict,
            time_bin_info_dict=time_bin_info_dict,
            convolution_instruction=convolution_instruction,
            #
            persistent_data=persistent_data,
            previous_convolution_results=previous_convolution_results,
        )
    else:
        convolution_results = convolve_pre_calculated_data(
            config=config,
            sfr_dict=sfr_dict,
            data_dict=data_dict,
            time_bin_info_dict=time_bin_info_dict,
            convolution_instruction=convolution_instruction,
            #
            persistent_data=persistent_data,
            previous_convolution_results=previous_convolution_results,
        )

    return convolution_results


def convolution_job_worker(job_queue, error_queue, worker_ID, config):  # DH0001
    """
    Function that handles running the job
    """

    setproctitle.setproctitle(
        "convolution multiprocessing worker process {}".format(worker_ID)
    )

    # Get items from the job_queue
    for job_dict in iter(job_queue.get, "STOP"):
        #########
        # Stopping or working
        if job_dict == "STOP":
            return None

        # TODO: most of the parts below are shared with the sequential convolution. Abstract

        # Unpack info
        convolution_instruction = job_dict["convolution_instruction"]
        data_dict = job_dict["data_dict"]
        time_bin_info_dict = job_dict["time_bin_info_dict"]
        sfr_dict = job_dict["sfr_dict"]

        job_dict["worker_ID"] = worker_ID

        ##########
        # Set up output dict
        payload = {}

        ##############
        # run convolution
        try:
            convolution_results = handle_convolution_choice(
                config=config,
                job_dict=job_dict,
                sfr_dict=sfr_dict,
                convolution_instruction=convolution_instruction,
                data_dict=data_dict,
            )

            ##############
            # Construct dictionary that is stored in the pickle files
            payload["bin_center"] = time_bin_info_dict["bin_center"]
            payload["convolution_instruction"] = convolution_instruction
            payload = {
                **payload,
                "convolution_results": convolution_results["convolution_results"],
            }

            ##############
            # Store info
            with open(
                os.path.join(
                    job_dict["output_dir"], "{}.p".format(payload["bin_center"])
                ),
                "wb",
            ) as f:
                pickle.dump(payload, f)

        ##############
        # handle errors
        except Exception as e:
            error_queue.put(
                (
                    "exception",
                    e,
                    "".join(traceback.format_tb(e.__traceback__)),
                    worker_ID,
                )
            )


def create_bin_iterator(config, convolution_instruction, sfr_dict):  # DH0001
    """
    Function to create the bin iterator data
    """

    ######
    # Determine bins to loop over
    # - backward conv loops over convolution bins
    # - forward conv loops over sfr bins

    ###
    # Support checks

    # integrate convolution does not support forward convolution (yet) TODO: not too difficult to support.
    if (convolution_instruction["convolution_type"] == "integrate") and (
        convolution_instruction["convolution_direction"] != "backward"
    ):
        raise ValueError(
            "Choice of convolution-method {} for convolution_type=`convolution_by_integration` is not supported. Only convolution_direction=`backward` is supported.".format(
                convolution_instruction["convolution_direction"]
            )
        )

    # on-the-fly convolution does not support backward convolution
    if (convolution_instruction["convolution_type"] == "on-the-fly") and (
        convolution_instruction["convolution_direction"] != "forward"
    ):
        raise ValueError(
            "Choice of convolution-method {} for convolution_type=`on-the-fly` is not supported. Only convolution_direction=`forward` is supported.".format(
                convolution_instruction["convolution_direction"]
            )
        )

    #
    if convolution_instruction["convolution_direction"] == "backward":
        bin_type = "convolution time"
        zipped_bin_data = zip(
            config["convolution_time_bin_centers"],
            config["convolution_time_bin_sizes"],
            config["convolution_time_bin_edges"][:-1],
        )
    elif convolution_instruction["convolution_direction"] == "forward":
        bin_type = "star formation time"
        zipped_bin_data = zip(
            sfr_dict["time_bin_centers"],
            sfr_dict["time_bin_sizes"],
            sfr_dict["time_bin_edges"][:-1],
        )
    else:
        raise ValueError(
            "`convolution_direction` {} not supported".format(
                convolution_instruction["convolution_direction"]
            )
        )

    # flip if we want to reverse convolution direction. Related to persistant data and previous results
    if convolution_instruction["reverse_convolution"]:
        zipped_bin_data = zipped_bin_data[::-1]

    return zipped_bin_data, bin_type


def convolution_queue_filler(  # DH0001
    job_queue,
    num_cores,
    config,
    sfr_dict,
    convolution_instruction,
    data_dict,
    processes,
):
    """
    Function to handle filling the queue for the multiprocessing

    When the convolution instruction is a sampling-based convolution,
    we use forward convolution, which loops over starformation bins
    rather than convolution bins
    """

    # Set up bin iterator data
    zipped_bin_data, bin_type = create_bin_iterator(
        config=config,
        convolution_instruction=convolution_instruction,
        sfr_dict=sfr_dict,
    )

    ######
    # Fill the queue with centres
    for bin_number, (
        bin_center,
        bin_size,
        bin_edge_lower,
    ) in enumerate(zipped_bin_data):

        # store current bin info, which is different in different cases.
        time_bin_info_dict = create_time_bin_info_dict(
            config=config,
            convolution_instruction=convolution_instruction,
            bin_number=bin_number,
            bin_center=bin_center,
            bin_edge_lower=bin_edge_lower,
            bin_size=bin_size,
            bin_type=bin_type,
        )

        # Set up job dict
        job_dict = create_job_dict(
            config=config,
            sfr_dict=sfr_dict,
            data_dict=data_dict,
            convolution_instruction=convolution_instruction,
            time_bin_info_dict=time_bin_info_dict,
            bin_number=bin_number,
        )

        #
        config["logger"].debug("job {} in the queue".format(job_dict["job_number"]))

        # Put job in queue
        job_queue.put(job_dict)

    # Signal stop to workers
    config["logger"].debug("Sending job termination signals")
    for _ in range(num_cores):
        job_queue.put("STOP")


def handle_multiprocessing_convolution(  # DH0001
    config, convolution_instruction, sfr_dict
):
    """
    Main function to handle convolution by multiprocessing

    This allows several cores to handle convolutions at the same time, but in
    this case the user cannot store persistent information and use the results
    of the previous convolution.
    """

    ###################
    # Set up data_dict: dictionary that contains the arrays or ensembles that are required for the convolution.
    config, data_dict, convolution_instruction = generate_data_dict(
        config=config, convolution_instruction=convolution_instruction
    )

    ###################
    # Run the convolution through multiprocessing

    # Set process name
    setproctitle.setproctitle("Convolution parent process")

    # Set up the manager object that can share info between processes
    manager = multiprocessing.Manager()
    job_queue = manager.Queue(config["max_job_queue_size"])
    error_queue = manager.Queue(config["max_job_queue_size"])

    # Create process instances
    processes = []
    for worker_ID in range(config["num_cores"]):
        processes.append(
            # Process(
            multiprocessing.Process(
                target=convolution_job_worker,
                args=(job_queue, error_queue, worker_ID, config),
            )
        )

    # Activate the processes
    for p in processes:
        p.start()

    # Start the system_queue and process
    convolution_queue_filler(
        job_queue=job_queue,
        num_cores=config["num_cores"],
        config=config,
        sfr_dict=sfr_dict,
        convolution_instruction=convolution_instruction,
        data_dict=data_dict,
        processes=processes,
    )

    # Join the processes to wrap up
    for p in processes:
        p.join()

    # Pass errors
    if not error_queue.empty():
        result_type, result_value, tb_string, worker_id = error_queue.get()
        if result_type == "exception":
            amended_args = tuple(
                [f"{result_value.args[0]}\n{str(tb_string)}", *result_value.args[1:]]
            )
            result_value.args = amended_args
            raise result_value


def handle_sequential_convolution(config, convolution_instruction, sfr_dict):  # DH0001
    """
    Main function to handle sequential convolution.

    This handles the convolution steps in sequence, but also allows the user
    to provide persistent information and use results of the previous
    convolution step
    """

    ###################
    # Set up data_dict: dictionary that contains the arrays or ensembles that
    # are required for the convolution.
    config, data_dict, convolution_instruction = generate_data_dict(
        config=config, convolution_instruction=convolution_instruction
    )

    # Set up bin iterator data
    zipped_bin_data, bin_type = create_bin_iterator(
        config=config,
        convolution_instruction=convolution_instruction,
        sfr_dict=sfr_dict,
    )

    # #############
    #
    persistent_data = {}
    previous_convolution_results = None

    # #############
    # TODO: this loop is shared with the queue filler. abstract
    # loop over bins
    for bin_number, (
        bin_center,
        bin_size,
        bin_edge_lower,
    ) in enumerate(zipped_bin_data):

        # store current bin info, which is different in different cases.
        time_bin_info_dict = create_time_bin_info_dict(
            config=config,
            convolution_instruction=convolution_instruction,
            bin_number=bin_number,
            bin_center=bin_center,
            bin_edge_lower=bin_edge_lower,
            bin_size=bin_size,
            bin_type=bin_type,
        )

        # Set up job dict
        job_dict = create_job_dict(
            config=config,
            sfr_dict=sfr_dict,
            data_dict=data_dict,
            convolution_instruction=convolution_instruction,
            time_bin_info_dict=time_bin_info_dict,
            bin_number=bin_number,
        )

        # #############
        # run convolution
        convolution_results = handle_convolution_choice(
            config=config,
            job_dict=job_dict,
            sfr_dict=sfr_dict,
            convolution_instruction=convolution_instruction,
            data_dict=data_dict,
            #
            persistent_data=persistent_data,
            previous_convolution_results=previous_convolution_results,
        )

        # Store previous results
        previous_convolution_results = copy.deepcopy(convolution_results)

        # add persistent data to the convolution_results that is stored
        if persistent_data is not None:
            if isinstance(persistent_data, dict):
                convolution_results["convolution_results"] = {
                    **convolution_results["convolution_results"],
                    **persistent_data,
                }
            else:
                raise ValueError(
                    "persistent_data ({}) should be a dictionary".format(
                        persistent_data
                    )
                )

        # #############
        # store information
        # TODO: below is copied quite roughly from the multiprocessing version. Should be cleaned

        # Get groupname
        groupname, _ = generate_group_name(
            convolution_instruction=convolution_instruction, sfr_dict=sfr_dict
        )
        full_groupname = "output_data/" + groupname

        #########
        # Handle storing convolution results
        with h5py.File(config["output_filename"], "a") as output_hdf5file:
            config["logger"].debug("Writing results to {}".format(full_groupname))

            # Readout group
            grp = output_hdf5file[full_groupname]

            # Handle storing
            handle_storing_convolution_results(
                config=config,
                grp=grp,
                bin_center=bin_center,
                convolution_results=convolution_results,
            )


def handle_sequential_or_multiprocessing_convolution(  # DH0001
    config, convolution_instruction, sfr_dict
):
    """ """

    if config["multiprocessing"] is True:
        handle_multiprocessing_convolution(
            config=config,
            convolution_instruction=convolution_instruction,
            sfr_dict=sfr_dict,
        )
    else:
        handle_sequential_convolution(
            config=config,
            convolution_instruction=convolution_instruction,
            sfr_dict=sfr_dict,
        )


def handle_convolution_steps(config, convolution_instruction, sfr_dict):  # DH0001
    """
    Function to handle the pre-convolution, convolution, and post-convolution steps for a particular set of SFR dict and convolution_instruction
    """

    # pre-convolution
    pre_convolution(
        config=config,
        convolution_instruction=convolution_instruction,
        sfr_dict=sfr_dict,
    )

    # actual convolution
    handle_sequential_or_multiprocessing_convolution(
        config=config,
        convolution_instruction=convolution_instruction,
        sfr_dict=sfr_dict,
    )

    #
    post_convolution(
        config=config,
        convolution_instruction=convolution_instruction,
        sfr_dict=sfr_dict,
    )


def convolve_populations(config):  # DH0001
    """
    Main function to handle the convolution of populations
    """

    #######
    # Check if we need to provide info for the SFR loop of not
    actual_sfr_dict_loop = False
    if isinstance(config["SFR_info"], dict):
        sfr_dicts = [config["SFR_info"]]
    else:
        actual_sfr_dict_loop = True
        sfr_dicts = config["SFR_info"]

    ########
    # Loop over all sfr dicts
    for sfr_dict_number, sfr_dict in enumerate(sfr_dicts):

        # provide info for sfr loop if necessary
        if actual_sfr_dict_loop:
            config["logger"].debug(
                "Handling SFR {} (number {}) ".format(sfr_dict["name"], sfr_dict_number)
            )

        ########
        # Convolution
        for convolution_instruction in config["convolution_instructions"]:

            ########
            # check if we chunk
            bound_handle_convolution_steps = partial(
                handle_convolution_steps,
                config=config,
                sfr_dict=sfr_dict,
            )
            if convolution_instruction["chunked_readout"]:

                # extract total number of chunks we should go over.
                total_chunk_number = convolution_instruction["chunk_total"]

                # loop over chunk
                for chunk in range(total_chunk_number):
                    convolution_instruction["chunk_number"] = chunk
                    bound_handle_convolution_steps(
                        convolution_instruction=convolution_instruction
                    )
            else:
                bound_handle_convolution_steps(
                    convolution_instruction=convolution_instruction
                )
