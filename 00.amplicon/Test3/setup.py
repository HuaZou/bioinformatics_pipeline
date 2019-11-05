\#!/usr/bin/env python

import argparse as ap 
import sys, os, re 

def parse_argument(args):
	parser = ap.ArgumentParser(description='amplicon sequencing pipeline')
	parser.add_argument(
		'-p', '--path', metavar='<sample path file>', type=str,
		help='sample path file (SampleID|fq1|fq2)', required=True)
	parser.add_argument(
		'-c', '--config', metavar='<configure file>', type=str, 
		help='provide a configure file including needed database and parameters for each setp',
		required=True)
	parser.add_argument(
		'-o', '--outdir', metavar='<output directory>', type=str,
		help='output directory path. Conatins the results and scripts.', required=True)

	return parser.parse_args()


###############################################
# init result directories
# Result:
#      Raw_data_status
#	   Trim_adaptor (FWD|REV)
#	   Filter_lowQC
#	   DADA2
#	   Taxonomy
#      Format_output
#      script
##############################################
def init_result_directory(directory):
	dir_list = {}
	total  = "/".join([os.getcwd(), directory])
	fqscan = "/".join([total, "00.scan_fastq"])
	trim   = "/".join([total, "01.trim_adaptor"])
	filt   = "/".join([total, "02.filter_reads"])
	dada2  = "/".join([total, "03.dada2_table"])
	tax    = "/".join([total, "04.taxonomy"])
	output = "/".join([total, "05.format_out"])
	script = "/".join([total, "script"])

	directory_list = [total, fqscan, trim, filt, \
			dada2, tax, output, script]
	name_list = ["total", "fqscan", "trim", "filt",\
			"dada2", "tax", "output", "script"]

	for name, dir in zip(name_list, directory_list):
		create_dir(dir)
		dir_list[name] = dir

	return dir_list


def create_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)
    return path


###################################################
# scripts to run 16s pipeline
###################################################
bin_dir = '/ldfssz1/ST_META/P17Z10200N0188_BGIXIAOPANG_CKY/zouhua/Practice/02.amplicon/Workflow/bin/'
fqscan_script = bin_dir + 'FqScan.py'
trim_script   = bin_dir + 'Trim.py'
filt_script   = bin_dir + 'Filter.R'
dada2_script  = bin_dir + 'Dada2.R'

###################################################
# get configure parameters for next script
###################################################
def get_config_parameters(configure_file):
	config_list = {}
	with open(configure_file, 'r') as cfg_f:
		for line in cfg_f:
			content = line.strip()
			if not content.startswith('#'):
				sample = re.split(r'\s+=\s+', content)
				config_list[sample[0]] = sample[1]
	
	return config_list


##################################################
# obtain the content of sample path file
##################################################
def obtain_fq_path(sample_file):
	sample_list = {}
	with open(sample_file, 'r') as sample_f:
		for line in sample_f:
			line = line.strip()
			samples = re.split(r'\s+', line)
			sample_list[samples[0]] = "\t".join([samples[1], samples[2]])
	
	return sample_list

#####################################################
# step1: generate scripts in check the fastq's status
#####################################################
def check_fastq(sample_dict, dir_dict, cfg_dict):
	check_fastq_dir   = dir_dict['fqscan']
	script_individual = dir_dict['script']
	script_total      = dir_dict['total']
	threads           = cfg_dict['threads'] 

	total_fastqc_file = script_total + "/" + "Run.S1_fastqc.sh"
	total_f = open(total_fastqc_file, 'w')
	for key, value in sample_dict.items():
		sampleid = key
		fastqs = re.split(r'\t', value)
		fq1 = fastqs[0]
		fq2 = fastqs[1]

		shell_fqc = " ".join(["python", fqscan_script, "-r1", fq1, "-r2", fq2, \
					'-t', threads, '-o', check_fastq_dir])
		ind_file = script_individual + "/" + sampleid + ".fastqscan.sh"
		with open(ind_file, 'w') as ind_f:
			ind_f.write(shell_fqc + "\n")
		
		total_f.write(shell_fqc + "\n")
	
	shell_mqc = " ".join(["multiqc", "--title", "fastqc", \
				"--module", "fastqc", check_fastq_dir, "--outdir",  check_fastq_dir])
	total_f.write(shell_mqc + "\n")

	return 'OK'

######################################################
# step2: generate scripts in trimming adapter sequence
######################################################
def cut_adapter(sample_dict, dir_dict, cfg_dict):
	cut_adapter_dir   = dir_dict['trim']
	script_individual = dir_dict['script']
	script_total      = dir_dict['total']
	ADAPTER_FWD       = cfg_dict['ADAPTER_FWD']
	ADAPTER_REV       = cfg_dict['ADAPTER_REV']
	match_bases       = cfg_dict['match_bases']
	error_rate        = cfg_dict['error_rate']
	min_length        = cfg_dict['min_length']
	 
	total_cut_file = script_total + "/" + "Run.S2_cutadapt.sh"
	total_f = open(total_cut_file, 'w')
	for key, value in sample_dict.items():
		sampleid = key
		fastqs = re.split(r'\t', value)
		fq1 = fastqs[0]
		fq2 = fastqs[1]

		shell_cut = " ".join(["python", trim_script, "-i", sampleid, "-r1", fq1, "-r2", fq2, \
					"-f", ADAPTER_FWD, "-r", ADAPTER_REV, "-b", match_bases, \
					"-e", error_rate, "-m", min_length, "-o", cut_adapter_dir])
		ind_file = script_individual + "/" + sampleid + ".cutadapter.sh"
		with open(ind_file, 'w') as ind_f:
			ind_f.write(shell_cut + "\n")
		
		total_f.write(shell_cut + "\n")

	total_f.close()

	return 'OK'

######################################################
# step3: filter low quality reads by dada2
######################################################
def filter_reads(dir_dict, cfg_dict):
	cut_adapter_dir   = dir_dict['trim']
	filt_read_dir     = dir_dict['filt']
	script_total      = dir_dict['total']
	truncLen_1        = cfg_dict['truncLen_1']
	truncLen_2        = cfg_dict['truncLen_2']
	maxEE             = cfg_dict['maxEE']
	truncQ            = cfg_dict['truncQ']
	maxN              = cfg_dict['maxN']

	filtpathF = cut_adapter_dir + "/FWD"
	filtpathR = cut_adapter_dir + "/REV"

	total_filt_file = script_total + "/" + "Run.S3_filter.sh"
	total_filt_f = open(total_filt_file, 'w')
	shell_filt = " ".join(["Rscript", filt_script, "-f", filtpathF, "-r", filtpathR, \
				"-t1", truncLen_1, "-t2", truncLen_2, \
				"-me", maxEE, "-tq", truncQ, "-mn", maxN, "-o", filt_read_dir])
	
	total_filt_f.write(shell_filt + "\n")
	total_filt_f.close()

	return "OK"

######################################################
# step4: infer sequence variance and assign taxonomy
######################################################
def Infer_feature_tax(dir_dict, cfg_dict):
	filt_read_dir = dir_dict['filt']
	dada_dir	  = dir_dict['dada2']
	tax_dir       = dir_dict['tax']		
	script_total  = dir_dict['total']
	database      = cfg_dict['database']


	filtpathF = filt_read_dir + "/FWD"
	filtpathR = filt_read_dir + "/REV"

	total_dada_file = script_total + "/" + "Run.S4_OTU.sh"
	total_dada_f = open(total_dada_file, 'w')
	shell_dada = " ".join(["Rscript", dada2_script, "-f", filtpathF, "-r", filtpathR, \
				"-o1", dada_dir, "-d", database, "-o2", tax_dir])
	total_dada_f.write(shell_dada + "\n")
	total_dada_f.close()

	return "OK"

######################################################
# step5: format output types and summary pipeline
######################################################

def joint_step(dir_dict):
	summary_step_file = dir_dict['total'] + "/Run.Total.sh" 
	with open(summary_step_file, 'w') as sum_f:
		sum_f.write("sh Run.S1_fastqc.sh\n")
		sum_f.write("sh Run.S2_cutadapt.sh\n")
		sum_f.write("sh Run.S3_filter.sh\n")
		sum_f.write("sh Run.S4_OTU.sh\n")



######################################################
# main: run each step
######################################################
def main():
	args = parse_argument(sys.argv)

	sample_file = args.path
	config_file = args.config
	out_dire    = args.outdir

	# samples path
	sample_lst = obtain_fq_path(sample_file)
	# output directory 
	dir_lst = init_result_directory(out_dire)
	# configure parameters
	config_lst = get_config_parameters(config_file)
	
	# step1 
	if check_fastq(sample_lst, dir_lst, config_lst):
		print("Step1: Scanning fastq succeed")
	else:
		print("Step1: Scanning fastq failed")

	# step2 
	if cut_adapter(sample_lst, dir_lst, config_lst):
		print("Step2: Cutting adapter succeed")
	else:
		print("Step2: Cutting adapter failed")
	
	# step3 
	if filter_reads(dir_lst, config_lst):
		print("Step3: filtering reads succeed")
	else:
		print("Step3: filtering reads failed")

	# step4 
	if Infer_feature_tax(dir_lst, config_lst):
		print("Step4: Running dada2 succeed")
	else:
		print("Step4: Running dada2 failed")

	joint_step(dir_lst)

if __name__ == '__main__':
	main()
