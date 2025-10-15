import os
import sys


def extract_times(file_path, start_str, result):
    with open(file_path, "r") as file:
        for line in file:
            if start_str in line:
                start_index = line.find(start_str) + len(start_str)
                end_index = line.find("s", start_index)

                substring = line[start_index:end_index]

                try:
                    valueAtTime = float(substring.strip().replace(",", "."))
                    result.append(valueAtTime)
                except ValueError:
                    print("CANNOT CONVERT SUBSTRING to float", substring.strip())

if len(sys.argv) != 2:
    print(
        "Usage: python extract_times.py <folder_path>\n"
        + "<folder_path> the path to the folder containing the log files\n"
    )
    sys.exit(1)

folder_path = sys.argv[1]

files = os.listdir(folder_path)
start_str = "INFO: compute() took"
output_path = "times.txt"

# Go through each file in the folder
times = []
for file_name in files:
    file_path = os.path.join(folder_path, file_name)
    extract_times(file_path, start_str, times)

if times:
    with open(output_path, "w") as file:
        average = sum(times) / len(times)
        file.write(f"Average: {average}s\n")
        for i, time in enumerate(times):
            file.write(f"Run {i}: {time}s\n")
    
    if len(times) != len(files):
        print(f"WARNING: Extracted {len(times)} run times, but there are {len(files)} files in total!")
    print(f"Average: {average}s of {len(times)} time(s): {times}")
