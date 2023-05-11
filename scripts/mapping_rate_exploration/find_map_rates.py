#!/usr/bin/env python3

__author__ = "Anna Pardo"

# import modules
import argparse
import os
import pandas as pd
import json
import statistics

# write a function to generate a dataframe with filepaths and mapping rates from Salmon
def make_df(salmon_dir):
        # salmon_dir is the full path to the directory containing Salmon outputs

        # initialize empty lists
	percent_mapped = []
	sample_info = []
	salmon_file = []

        # iterate through subdirectories of salmon_dir
	for folder in os.listdir(salmon_dir):
                # this is added for data downloaded from SRA
		if folder.startswith("SRR"):
                        # make a variable containing the path to the meta_info.json file
			json_file = salmon_dir+folder+"/aux_info"+"/meta_info.json"
			# open meta_info.json, pull out path to quant.sf file, sample name, and mapping rate, and append to previously initialized lists
			with open(json_file, "r") as infile:
				meta_info = json.load(infile)
				map_percent = meta_info["percent_mapped"]
				dir_split = json_file.strip().split("/")
				name_info = dir_split[8]
				sample_info.append(name_info)
				percent_mapped.append(map_percent)
				salmon_file.append(salmon_dir+folder+"/quant.sf")
        # make dataframe and return dataframe
	df = pd.DataFrame({"File":salmon_file,"Sample":sample_info,"Percent_Mapped":percent_mapped})
	return df

def main():
	parser = argparse.ArgumentParser(description="Parse args")
        # command line arguments
	parser.add_argument("--directory","-d",type=str,help="full path to directory containing salmon output files")
	parser.add_argument("--outfilename","-o",type=str,help="full path to desired output file for mapping rate json")
	args = parser.parse_args()
	directory = str(args.directory)
	outfilename = str(args.outfilename)

	df = make_df(directory)

	# print the summary statistics of mapping rates
	print("maximum mapping rate:",max(df["Percent_Mapped"]))
	print("minimum mapping rate:",min(df["Percent_Mapped"]))
	print("mean mapping rate:",sum(df["Percent_Mapped"])/len(df.index))
	print("median mapping rate:",statistics.median(df["Percent_Mapped"]))

	# save the dataframe to a file
	df.to_csv(outfilename,sep=",",index=False,header=True)

if __name__ == "__main__":
	main()
