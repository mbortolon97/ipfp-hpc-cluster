import os
import pandas as pd

dir_path = "/home/matteo/uniTN/hpc/results_only_rows"

id = []
number_of_processes = []
process_per_node = []

data = []
for basename in os.listdir(dir_path):
    file_path = os.path.join(dir_path, basename)
    root_ext = os.path.splitext(basename)
    if not root_ext[-1].startswith(".o"):
        continue
    fields = {'id': root_ext[-1][2:]}
    with open(file_path) as f:
        line = "Hello"
        while line:
            line = f.readline()
            
            if line.startswith("Number of processes: "):
                fields['number_of_processes'] = line[len("Number of processes: "):].strip()
            if line.startswith("Process per node: "):
                fields['process_per_node'] = line[len("Process per node: "):].strip()
            
            if line.startswith("Average time per hour: "):
                fields['average_time_per_hour'] = line[len("Average time per hour: "):].strip()
            if line.startswith("Average saving time per hour: "):
                fields['average_saving_time_per_hour'] = line[len("Average saving time per hour: "):].strip()
            if line.startswith("Average aggregation distribution time per hour: "):
                fields['average_aggregation_distribution_time_per_hour'] = line[len("Average aggregation distribution time per hour: "):].strip()
            
            if line.startswith("io_time: "):
                fields['io_time'] = line[len("io_time: "):].strip()
            if line.startswith("permutation_time: "):
                fields['permutation_time'] = line[len("permutation_time: "):].strip()
            if line.startswith("distribution_operations_time: "):
                fields['distribution_operations_time'] = line[len("distribution_operations_time: "):].strip()
            if line.startswith("Hardware clock tick: "):
                fields['hardware_clock_tick'] = line[len("Hardware clock tick: "):].strip()
    
    if len(fields) < 8:
        print("Invalid file: ", file_path)
        continue

    data.append(fields)

df = pd.DataFrame(data)
df.to_csv("result-summary.csv", index=False)

