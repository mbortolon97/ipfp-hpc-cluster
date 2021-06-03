import os

id = []
number_of_processes = []
process_per_node = []

for basename in os.path.listdir("/home/matteo/uniTN/hpc/results"):
    file_path = os.path.join("/home/matteo/uniTN/hpc/results", basename)
    root_ext = os.path.splitext(basename)
    if not root_ext[-1].startswith("o"):
        continue
    fields = {'id': root_ext[-1][1:]}
    with open(file_path) as f:
        line = f.readline()
        while line:
            line = f.readline()
            
            if line.startswith("Number of processes: "):
                fields['number_of_processes'] = line[len("Number of processes: "):]
            if line.startswith("Process per node: "):
                fields['process_per_node'] = line[len("Process per node: "):]
            
            if line.startswith("Average time per hour: "):
                fields['average_time_per_hour'] = line[len("Average time per hour: "):]
            if line.startswith("Average saving time per hour: "):
                fields['average_saving_time_per_hour'] = line[len("Average saving time per hour: "):]

            if line.startswith("io_time: "):
                fields['io_time'] = line[len("io_time: "):]
            if line.startswith("permutation_time: "):
                fields['permutation_time'] = line[len("permutation_time: "):]
            if line.startswith("distribution_operations_time: "):
                fields['distribution_operations_time'] = line[len("distribution_operations_time: "):]
            
    if len(fields) != 8:
        print("Invalid file: ", file_path)

