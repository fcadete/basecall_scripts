
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
		 'adapters':[],
                 'start_barcodes':[],
		 'end_barcodes':[],
 		 'final_barcode_call':[]}

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
			if line.startswith("Barcodes:"):
				line = porefile.readline().strip()
				line = line.replace("start barcodes: ", "").strip()
				read_adapters['start_barcodes'].append(line.split(", "))
				line = porefile.readline().strip()
				line = line.replace("end barcodes: ", "").strip()
				read_adapters['end_barcodes'].append(line.split(", "))
				line = porefile.readline()
				line = porefile.readline()
				line = porefile.readline().strip()
				read_adapters['final_barcode_call'].append(line.split(":")[1].strip())
			elif line.startswith("end"):
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

output_file = open(path_to_poreoutput.replace("output", "finalcall_table"), "w")

output_file.write("\t".join(['readid', 'runid', 'sampleid', 'read', 'ch', 'start_time', 'final_barcode_call']) + "\n")
for i in range(len(read_adapters['readid'])):
	string_output = []
	for read_key in read_adapters.keys():
		if read_key not in ['adapters', 'start_barcodes', 'end_barcodes']:
			string_output.append(read_adapters[read_key][i])
	#	elif read_key == 'adapters':
	#		adapters_list = []
	#		for j in range(0, len(read_adapters[read_key][i])):
	#			adapters_list.append(";".join(list(read_adapters[read_key][i][0].values())))
	#		string_output.append(":".join(adapters_list))
	#	else:
	#		string_output.append(";".join(read_adapters[read_key][i]))
	output_file.write("\t".join(string_output) + "\n")

output_file.close()


for i in range(len(read_adapters['readid'])):
	if len(read_adapters['adapters'][i]) > 0:
		break

output_file = open(path_to_poreoutput.replace("output", "adapters_table"), "w")

output_file.write("\t".join(['readid', 'runid', 'sampleid', 'read'] + list(read_adapters['adapters'][i][0].keys())) + "\n")
for i in range(len(read_adapters['readid'])):
	for j in range(len(read_adapters['adapters'][i])):
		string_output = []
		for read_key in ['readid', 'runid', 'sampleid', 'read']:
			string_output.append(read_adapters[read_key][i])
		output_file.write("\t".join(string_output + list(read_adapters['adapters'][i][j].values())) + "\n")

output_file.close()

output_file = open(path_to_poreoutput.replace("output", "barcodes_table"), "w")

output_file.write("\t".join(['readid', 'runid', 'sampleid', 'read', 'which_end', 'barcode', 'score']) + "\n")
for i in range(len(read_adapters['readid'])):
	for j in range(len(read_adapters['start_barcodes'][i])):
		string_output = []
		for read_key in ['readid', 'runid', 'sampleid', 'read']:
			string_output.append(read_adapters[read_key][i])
		string_output.append('start')
		
		output_file.write("\t".join(string_output +
						read_adapters['start_barcodes'][i][j].replace("(", "").replace("%)", "").split(" ")) +
					"\n")

output_file.close()


