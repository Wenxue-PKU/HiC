#!/bin/bash
# Prepare Juicer's restriction files for all "sites.txt" file

### Usage
# sh /home/hidekit/Library/HiC/Package/JUICER/Prepare_JUICER_restriction_files_for_all_data.sh


for f in $(ls *sites.txt)
do
	FILE_convert=${f/.txt/.juicer.txt}
	if [ ! -e $FILE_convert ]; then
		echo "Convert $f"
		perl /home/hidekit/Library/HiC/Package/JUICER/Convert_restriction_file.pl -i $f -o $FILE_convert
	else
		echo "$FILE_convert is already exists"
	fi
done
