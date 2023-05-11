#!/usr/bin/env python3

__author__ = "Anna Pardo"

import argparse
import os
import pandas as pd
import json
import statistics

def make_df(salmon_dir):
	percent_mapped = []
	sample_info = []
	salmon_file = []
	for folder in os.listdir(salmon_dir):
		if folder.startswith("SRR"):
			json_file = salmon_dir+folder+"/aux_info"+"/meta_info.json"
			#print(json_file)
			with open(json_file, "r") as infile:
				meta_info = json.load(infile)
				map_percent = meta_info["percent_mapped"]
				dir_split = json_file.strip().split("/")
				#print(dir_split)
				name_info = dir_split[8]
				print(name_info)
				#if name_info.endswith("_1"):
				#	name_info = name_info.split("_")[0]
				sample_info.append(name_info)
				percent_mapped.append(map_percent)
				salmon_file.append(salmon_dir+folder+"/quant.sf")
	df = pd.DataFrame({"File":salmon_file,"Sample":sample_info,"Percent_Mapped":percent_mapped})
	return df

def main():
	parser = argparse.ArgumentParser(description="Parse args")
	parser.add_argument("--directory","-d",type=str,help="full path to directory containing salmon output files")
	parser.add_argument("--outfilename","-o",type=str,help="full path to desired output file for mapping rate json")
	args = parser.parse_args()
	directory = str(args.directory)
	outfilename = str(args.outfilename)

	df = make_df(directory)
	print(df.head())

	# print the statistics
	print("maximum mapping rate:",max(df["Percent_Mapped"]))
	print("minimum mapping rate:",min(df["Percent_Mapped"]))
	print("mean mapping rate:",sum(df["Percent_Mapped"])/len(df.index))
	print("median mapping rate:",statistics.median(df["Percent_Mapped"]))

	# save the dataframe to a file
	df.to_csv(outfilename,sep=",",index=False,header=True)

if __name__ == "__main__":
	main()
