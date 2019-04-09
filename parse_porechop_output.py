
import re
import sys

path_to_poreoutput = sys.argv[1]

porefile = open(path_to_poreoutput, "r")

read_adapters = {'readid':[],
                 'runid':[],
                 'sampleid':[],
                 'read':[],
                 'ch':[],
                 'start_time':[],
		 'adapters':[]}

for line in porefile:
	if "reads had adapters trimmed" in line:
		break
	if re.match('\w{8}\-\w{4}\-', line):
		read_info = line.split(" ")
		read_adapters['readid'].append(read_info[0])
		for i in range(1, len(read_info)):
			split_info = read_info[i].split('=')
			read_adapters[split_info[0]].append(split_info[1].rstrip())
		adapters = []
		which_end = "start"
		line = porefile.readline().strip()
		while line:
			if line.startswith("end"):
				which_end = "end"
			elif line.startswith("start") == False:
				adapter_info = line.split(", ")
				this_adapter = {'which_end':which_end,
						'adapter_name':adapter_info[0]}
				for i in range(1, len(adapter_info)):
					split_info = re.split('=|: ', adapter_info[i])
					this_adapter[split_info[0]] = split_info[1]
				adapters.append(this_adapter)
			line = porefile.readline().strip()
		read_adapters['adapters'].append(adapters)

output_file = open(path_to_poreoutput.replace("output", "table"), "w")


for i in range(len(read_adapters['readid'])):
	if len(read_adapters['adapters'][i]) > 0:
		break


output_file.write("\t".join(list(read_adapters.keys()) + list(read_adapters['adapters'][i][0].keys())) + "\n")
for i in range(len(read_adapters['readid'])):
	for j in range(len(read_adapters['adapters'][i])):
		string_output = []
		for read_key in read_adapters.keys():
			if read_key != 'adapters':
				string_output.append(read_adapters[read_key][i])
		output_file.write("\t".join(string_output + list(read_adapters['adapters'][i][j].values())) + "\n")

output_file.close()


