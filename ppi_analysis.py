import csv
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import os
import sys
import argparse
import itertools as it

files = [
    "names.csv",
    "ppi.csv",
    "newppi.csv",
    "fullcorr.csv",
    "Affinity Capture-Luminescence.csv",
    "Affinity Capture-MS.csv",
    "Affinity Capture-RNA.csv",
    "Affinity Capture-Western.csv",
    "Biochemical Activity.csv",
    "Co-crystal Structure.csv",
    "Co-fractionation.csv",
    "Co-localization.csv",
    "Co-purification.csv",
    "Dosage Growth Defect.csv",
    "Dosage Lethality.csv",
    "Dosage Rescue.csv",
    "FRET.csv",
    "Far Western.csv",
    "Negative Genetic.csv",
    "PCA.csv",
    "Phenotypic Enhancement.csv",
    "Phenotypic Suppression.csv",
    "Positive Genetic.csv",
    "Protein-RNA.csv",
    "Protein-peptide.csv",
    "Proximity Label-MS.csv",
    "Reconstituted Complex.csv",
    "Synthetic Growth Defect.csv",
    "Synthetic Haploinsufficiency.csv",
    "Synthetic Lethality.csv",
    "Synthetic Rescue.csv",
    "Two-hybrid.csv",
]


def create_arrays(filename):
    """
    returns input file as a list or np array
    """
    with files_loaded.get_lock():
        files_loaded.value += 1
    print(f"\rloading file {files_loaded.value} of {len(files)}...", end="")
    if filename == "names.csv":
        return list(csv.reader(open(filename, "r"), delimiter=","))[0]
    elif filename.islower():
        return np.array(
            list(csv.reader(open(filename, "r"), delimiter=",")), dtype=np.float
        )
    else:
        return np.array(
            list(
                csv.reader(
                    open(f"biogrid_arrays{os.path.sep}{filename}", "r"), delimiter=","
                )
            ),
            dtype=np.int8,
        )


def biogrid_verification(test_array):
    """
    returns a list of the percentage of the input matrix that overlaps with
    1. each BioGrid method's sparse matrix (28 total)
    2. the union of all these matrices
    3. inverse of this union (1's and 0's flipped)
    INPUT
    test_array must be 2708x2078 np array
    OUTPUT
    results is a 30 item list containing overlap percentages
    """
    # get the sum of the input array; our matrices are symmetric
    # so their sums must be halved to get the correct count
    sum_test_array = np.sum(test_array) / 2
    results = []
    # iterate through each biogrid matrix, the union of all matrices, and the
    # inverse of the union
    for array in it.chain(biogrid.values(), (all_biogrid_data, not_biogrid_data)):
        # overlap of each matrix and the input matrix
        intersection = array & test_array
        sum_intersection = np.sum(intersection) / 2
        # overlap between matrices is recorded as a percentage of the sum of
        # the input matrix
        results.append((sum_intersection / sum_test_array) * 100)
    return results


def edge_count(cutoff):
    """
    creates new PPI matrix with the cutoff given as an argument
    and returns the number of edges in the new matrix
    INPUT
    cutoff could be a single value or an np array of values, both work with the syntax
    OUTPUT
    single number representing the matrix edge count
    """
    # create new PPI matrix with the current cutoff
    new_matrix = (fullcorr >= cutoff).astype(np.int8)
    return np.sum(new_matrix)


def cutoff_test(cutoff):
    """
    creates new PPI matrix with the cutoff given as an argument
    and returns the results of biogrid verification on the new matrix
    INPUT
    cutoff could be a single value or an np array of values, both work with the syntax
    OUTPUT
    returns the 30 item list from verification function
    """
    # create new PPI matrix with the current cutoff
    new_matrix = (fullcorr >= cutoff).astype(np.int8)
    # perform biogrid verification on the new PPI matrix
    return biogrid_verification(new_matrix)


def correlation_testing():
    """
    generates a range of correlation cutoffs, performs testing on each
    cutoff and plots the results
    """
    # test edge prediction of the full correlation matrix at varying cutoffs
    cutoff_values = np.arange(0.01, 1.01, 0.01)
    pool = mp.Pool(num_cores)
    print("correlation coefficient")
    print("generating new networks and testing them...")

    # get edge counts for new networks produced from cutoffs.
    # edge_counts will be a 100x1 matrix containing the total number of
    # edges for each of the 100 cutoff value networks tested.
    edgecount_results = pool.map(edge_count, cutoff_values)
    edge_counts = np.array(edgecount_results, dtype=np.int32)
    edge_counts = edge_counts.transpose()

    # get PPI overlap for new networks produced from cutoffs.
    # overlap will be a 30x100 matrix containing the overlap percentages
    # for each of the 100 new networks with the 30 biogrid verification methods
    cutoff_results = pool.map(cutoff_test, cutoff_values)
    overlap = np.array(cutoff_results, dtype=np.float)
    overlap = overlap.transpose()

    pool.close()
    pool.join()
    print("finished testing!")
    plot_results(edge_counts, overlap, "correlation coefficient")


def plot_results(edge_counts, results, test="unknown"):
    """
    plot the results of biogrid verification against the number of network edges.
    test -> text description of the evaluation metric used to generate the new network
    """
    # if you're just doing a comparison summary, make a few bits of the
    # data accessible and return from the function
    if compare == True:
        global overall_results
        new_summary = []
        new_summary.append(test)
        new_summary.append(edge_counts)
        new_summary.append(results[-2])
        new_summary.append(results[-1])
        overall_results.append(new_summary)
        return

    print("plotting the results...")
    made_dir = False
    for x in range(results.shape[0]):
        plt.plot(edge_counts, results[x], "-")
        if x == results.shape[0] - 2:
            plt.title(f"PPIs verified by any BioGrid method\n{test}")
            img_title = "verified_0any"
        elif x == results.shape[0] - 1:
            plt.title(f"PPIs unverified\n{test}")
            img_title = "unverified"
        else:
            plt.title(f"Verified by {biogrid_methods[x]}\n{test}")
            img_title = "verified_" + biogrid_methods[x].replace(" ", "_")
        plt.xlabel("edge count")
        plt.ylabel("correct edge %")
        plt.ylim([0, 100])
        plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
        if (
            is_cmd
        ):  # if the program is running from the command-line, save plots to directory
            if not made_dir:
                timestamp = datetime.now().strftime(f"{test.split()[0]}_%y%m%d%H%M%S")
                dir_name = f"image_results{os.path.sep}{timestamp}"
                os.makedirs(dir_name)
                made_dir = True
            plt.savefig(f"{dir_name}{os.path.sep}{img_title}")
            plt.close()
        else:  # interactive program only shows plots on screen
            plt.show()
    print("finished plotting!")
    print()


def plot_summary():
    """
    plot the results of biogrid verification for verification by any method as
    well as unverified results. this is done for all network generation methods,
    allowing a direct comparison of their efficacy.
    overall_results -> list of lists; each list contains the following 4 items:
        1. Name of testing method
        2. NP array of the edge counts for each network generated from a parameter value
        3. NP array of biogrid verification results for PPIs verified by any biogrid method
        4. NP array of biogrid verification results for PPIs unverified by all biogrid methods
        5. Plot color - verified by any biogrid data
        6. Plot color - unverified by biogrid data
    """
    print("plotting the results...")
    plt.figure(figsize=(12, 6))
    for info in overall_results:
        legend_label = "Strategy "
        if info[0][0] == "c":
            legend_label += "1"
        elif info[0][0] == "t":
            legend_label += "2"
        elif info[0][0] == "o":
            legend_label += "3"
        elif info[0][0] == "g":
            legend_label += "4"
        elif info[0][0] == "l":
            legend_label += "5"
        plt.plot(
            info[1],
            info[2],
            color=info[4],
            linestyle="-",
            label=f"{legend_label} - Verified",
        )
        plt.plot(
            info[1],
            info[3],
            color=info[5],
            linestyle="-",
            label=f"{legend_label} - Unverified",
        )
    plt.xlabel("edge count")
    plt.xlim([0, 120_000])
    plt.ylabel("correct edge %")
    plt.ylim([0, 100])
    plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    plt.legend(loc="center right")
    if (
        is_cmd
    ):  # if the program is running from the command-line, save summary plot to directory
        timestamp = datetime.now().strftime(f"summary_%y%m%d%H%M%S")
        dir_name = "image_results"
        plt.savefig(f"{dir_name}{os.path.sep}{timestamp}")
        plt.close()
    else:  # interactive program only shows plots on screen
        plt.show()
    print("finished plotting!")
    print()


def generate_coefficients(degree_array):
    """
    takes a range of row indices as input and for each row, calculates
    the correlation cutoff for each node that will yield the desired degree,
    and returns
    INPUT
    degree_array is a 1x2708 np array containing the degrees to take for each node
    OUTPUT
    results is a 2708 item list containing the appropriate correlation coefficients
    to use for each node in order to get the desired degree
    """
    results = []
    for row in range(fullcorr.shape[0]):
        # copy the corresponding row of the full correlation matrix
        # and sort in ascending order
        sorted_row = np.sort(fullcorr[row].copy())
        # use the row's degree in the original PPI matrix to
        # determine the cutoff that, when applied to the full
        # correlation matrix, will yield a row of the same degree
        try:
            cutoff = sorted_row[-degree_array[row]]
        except:
            cutoff = sorted_row[0]
        results.append(cutoff)
    return results


def degree_testing(mode):
    """
    "original" mode -> keep degrees the same as the original PPI matrix
    "top" mode -> select top N most correlated edges for each node
    "global" mode -> adjust the number of edges based on global average
    "local" mode -> adjust the number of edges based on local average

    degree_array -> indices for each node to select the desired number of
    degrees for the given mode
    """
    print(mode)
    print("generating new networks and testing them...")
    degree_list = []
    if mode == "top":
        test_text = "top degrees"
        degree_array = np.zeros((fullcorr.shape[0],), dtype=np.int16)
        iter_list = np.arange(1, 131)
        # after the for loop, degree_list will contain 130 lists, each with 2708 items.
        # biogrid overlap for top 1-130 degrees, for each of the 2708 nodes
        for _ in iter_list:
            degree_array += 1
            degree_list.append(generate_coefficients(degree_array))
    else:
        # "original", "global", "local" modes
        # create a boolean matrix from the original PPI matrix
        orig_ppi = oldppi > 0
        # calculate the degree of each row
        degree_array = np.sum(orig_ppi, axis=1, dtype=np.int16)
        if mode == "original":
            test_text = "original number of degrees"
            iter_list = np.arange(0.1, 10.1, 0.1)
            for adjustment in iter_list:
                tmp_degree_array = (
                    adjustment * degree_array.astype(np.float64)
                ).astype(np.int16)
                # make sure the values of the tmp_degree_array cannot be less than 1
                tmp_degree_array = np.clip(tmp_degree_array, 1, None)
                # after this for loop, degree_list will contain 100 lists, one for each
                # adjustment. each list will contain 2708 items
                degree_list.append(generate_coefficients(tmp_degree_array))
        else:
            # adjust original PPI degrees based on a global or local average
            if mode == "global":
                iter_list = np.arange(0.1, 10.1, 0.1)
                test_text = "global average adjustment"
                # average is a vector here
                average = np.ones((fullcorr.shape[0],), dtype=np.float64)
                average *= np.mean(degree_array)
            elif mode == "local":
                iter_list = np.arange(0.1, 10.1, 0.1)
                test_text = "local average adjustment"
                # average is a vector here
                average = np.zeros((fullcorr.shape[0],), dtype=np.float64)
                for x in range(fullcorr.shape[0]):
                    if degree_array[x] == 0:
                        average[x] = 0
                    else:
                        sum = np.sum(degree_array[orig_ppi[x]])
                        average[x] = sum / degree_array[x]
            for adjustment in iter_list:
                tmp_degree_array = degree_array + np.around(
                    adjustment * (average - degree_array.astype(np.float64))
                ).astype(np.int16)
                # make sure the values of the tmp_degree_array cannot be less than 1
                tmp_degree_array = np.clip(tmp_degree_array, 1, None)
                # after this for loop, degree_list will contain 100 lists, one for each
                # adjustment. each list will contain 2708 items
                degree_list.append(generate_coefficients(tmp_degree_array))
    # degree_np_array is a large array filled with correlation coefficients that are used to generate
    # the desired node degrees.
    degree_np_array = np.array(degree_list, dtype=np.float)
    # multi-process worker pool
    pool = mp.Pool(num_cores)

    # get edge counts for new networks produced from cutoffs
    edgecount_results = pool.map(edge_count, degree_np_array)
    edge_counts = np.array(edgecount_results, dtype=np.int32)
    edge_counts = edge_counts.transpose()

    # get PPI overlap for new networks produced from cutoffs
    cutoff_results = pool.map(cutoff_test, degree_np_array)
    overlap = np.array(cutoff_results, dtype=np.float)
    overlap = overlap.transpose()

    pool.close()
    pool.join()
    print("finished testing!")
    plot_results(edge_counts, overlap, test_text)


def compare_all():
    """
    generates a comparison summary for different network generation strategies
    """
    global compare, overall_results
    compare = True
    overall_results = []
    colors = [
        ["#ff0000", "#ff7070"],  # "correlation" - red
        ["#18b300", "#57b548"],  # "top" - green
        ["#ff9500", "#ffbc5e"],  # "original" - orange
        ["#004cff", "#6b97ff"],  # "global" - blue
        ["#7300ff", "#a963ff"],  # "local" - purple
    ]
    correlation_testing()
    degree_testing("top")
    degree_testing("original")
    degree_testing("global")
    degree_testing("local")
    
    for idx, pair in enumerate(colors):
        overall_results[idx].append(pair[0])
        overall_results[idx].append(pair[1])
    plot_summary()


def data_load():
    """
    loads the data needed to perform the subsequent analyses
    """
    global files_loaded, num_cores, names, oldppi, newppi, fullcorr, biogrid, biogrid_methods, all_biogrid_data, not_biogrid_data
    files_loaded = mp.Value("B", 0)
    # import all data as lists or np arrays using a multi-process worker pool
    num_cores = mp.cpu_count()
    pool = mp.Pool(num_cores)
    # names -> every protein name in the original PPI network
    # oldppi -> matrix of PPI strength for the original PPI network
    # newppi -> matrix of highest correlated PPIs, keeping the same number of
    # edges as the original network
    # fullcorr -> full correlation matrix created from random walk algorithm
    # biogrid_arrays -> for each method, matrix of PPIs verified by that method
    names, oldppi, newppi, fullcorr, *biogrid_arrays = pool.map(create_arrays, files)
    # create the biogrid as a dictionary of np arrays
    biogrid_methods = [x.replace(".csv", "") for x in files[4:]]
    biogrid = dict(zip(biogrid_methods, biogrid_arrays))

    pool.close()
    pool.join()

    # fill the main diagonal of the full correlation matrix with zero
    np.fill_diagonal(fullcorr, 0)

    # create the data needed for biogrid verification
    # 1.) union of all PPIs in the biogrid data into one array
    all_biogrid_data = np.zeros(fullcorr.shape, dtype=np.int8)
    for array in biogrid.values():
        all_biogrid_data |= array
    # 2.) create an inversion of this matrix (for measuring the amount of
    # unverified data)
    not_biogrid_data = all_biogrid_data ^ 1  # XOR with 1 to invert the matrix
    print("\rfinished loading data!  ")


if __name__ == "__main__":
    global compare
    compare = False
    if sys.stdin.isatty():
        # PROGRAM RUNNING FROM THE COMMAND-LINE
        is_cmd = True
        parser = argparse.ArgumentParser(
            description="Create new PPI networks based on various criteria and validate \
            against biological data."
        )
        parser.add_argument(
            "-a",
            help="Perform all testing and generate a brief comparison of the results.",
            action="store_true",
        )
        parser.add_argument(
            "-c",
            help="Create new networks by applying various correlation coefficient cutoff values for the network's \
                edges.",
            action="store_true",
        )
        parser.add_argument(
            "-o",
            help="Create new networks with the most correlated network edges such that each node has a degree that \
                is equal to its degree in the original PPI network multiplied by an adjustment factor.",
            action="store_true",
        )
        parser.add_argument(
            "-t",
            help="Create new networks with the most correlated network edges such that each node has a degree that \
                is equal to NUMBER_OF_DEGREES.",
            action="store_true",
        )
        parser.add_argument(
            "-g",
            help="Create new networks with the most correlated network edges such that each node has a degree that \
                is equal to the global average multiplied by an adjustment factor.",
            action="store_true",
        )
        parser.add_argument(
            "-l",
            help="Create new networks with the most correlated network edges such that each node has a degree that \
                is equal to the local average multiplied by an adjustment factor.",
            action="store_true",
        )
        args = parser.parse_args()
        data_load()
        if args.a:
            compare_all()
        if args.c:
            correlation_testing()
        if args.o:
            degree_testing("original")
        if args.t:
            degree_testing("top")
        if args.g:
            degree_testing("global")
        if args.l:
            degree_testing("local")
    else:
        # PROGRAM RUNNING INTERACTIVELY
        is_cmd = False
        data_load()
        print("Perform testing by executing any of the following functions:")
        print("\tcorrelation_testing()")
        print("\tdegree_testing(mode)")
        print('\t\tmode : { "original", "top", "global", "local" }')
        print("Perform all testing with a brief comparison of the results:")
        print("\tcompare_all()")
