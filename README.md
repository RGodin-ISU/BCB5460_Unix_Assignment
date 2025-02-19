# UNIX Assignment 

## Data Inspection

### Attributes of `fang_et_al_genotypes`

```
# Changed directory into raw data folder.
cd ./raw_data

# I checked how large the file size is in human-readable format.
du -h fang_et_al_genotypes.txt

# I checked how many lines, words, and number of characters there were in the file.
wc fang_et_al_genotypes.txt

# I inspected the first two lines of the file to see what they looked like.
head -n 2 fang_et_al_genotypes.txt

# I checked the last two lines of the file to see if they included any metadata.
tail -n 2 fang_et_al_genotypes.txt

# I checked how many columns there were in the file, assuming it is tab-delimitered.
awk -F "\t" '{print NF; exit}' fang_et_al_genotypes.txt

# Since there were so many columns, I only wanted to see what the first 20 columns
# looked like for the first 20 lines in a human readable format.
cut -f 1-20 fang_et_al_genotypes.txt | column -t | head -n 20

# I checked if there were any lines which didn't have missing data (?/?) using grep.
grep "^?/?" fang_et_al_genotypes.txt

# I checked what the encoding of the file type was:
file fang_et_al_genotypes.txt

# I checked if there were any duplicate gene types
cut -f 2 fang_et_al_genotypes.txt | sort | uniq -c | sort -n | tail -n 3

# I checked what the most common group types were
cut -f 3 fang_et_al_genotypes.txt | sort | uniq -c | sort -nr | head -n 3
```

By inspecting this file I learned that:

* There are 2,783 lines in the file, 2744038 words, and 11051939 characters, with a file size of 6.3 MB. Thus, there are not many lines, but there appear to be a lot of words and characters relatively speaking. This is due to there being a large number of columns on each line, one for each single nucleotide polymorphisms (SNP) in the sample along with some metadata.
* Since only the first line of the file contained metadata with the head camp, and since the tail doesn't show any metadata, it appears there is only metadata for the first row. This implies there are 2,782 samples in the SNP dataset.
* Inspecting the number of columns, I got 986. Since the first three columns are metadata (sample ID, JG_OTU, and group), this appears to imply each sample has info on 983 SNPs.
* There is missing data for some of the SNPs with "?/?" appearing. At least, I assume this represents missing data. I couldn't find any sample which did not have "?/?" for at least one SNP using grep. Thus, I assume this is pretty common.
* The file is ASCII text, with very long lines.
* There were no duplicate gene types, but there were duplicate groups with the largest being the ZMMLR group at 1256 members.

### Attributes of `snp_position.txt`

```
# Changed directory into raw data folder.
cd ./raw_data

# I checked how large the file size is in human-readable format.
du -h snp_position.txt

# I checked how many lines, words, and number of characters there were in the file.
wc snp_position.txt

# I inspected the first two lines of the file to see what they looked like.
head -n 2 snp_position.txt

# I checked the last two lines of the file to see if they included any metadata.
tail -n 2 snp_position.txt

# I checked how many columns there were in the file, assuming it is tab-delimitered.
awk -F "\t" '{print NF; exit}' snp_position.txt

# Since the column names were so long, I only wanted to see what the first 9 columns
# looked like for the first 15 lines in a human readable format.
cut -f 1-9 snp_position.txt | column -t | head -n 15

# I also checked the end columns to see anything I missed
cut -f 9- snp_position.txt | column -t | head -n 15

# I checked what the encoding of the file type was:
file snp_position.txt

# I checked what the most common chromosomes were
cut -f 3 snp_position.txt | sort | uniq -c | sort -nr | head -n 3

# I checked what the most common genes were
cut -f 9 snp_position.txt | sort | uniq -c | sort -nr | head -n 3
```

By inspecting this file I learned that:

* There are 984 lines in the file, 13198 words, and 82763 characters, with a file size of 4.3 KB. Thus, the file is considerably smaller than the fang_et_al_genotypes.txt file as expected.
* The file doesn't appear to have any metadata with the only the first line providing the column names, with there appearing to be 15 columns. Thus, there appear to be 983 individual SNPs documented in the file. For each SNP, there appears to be 
	* From looking at the first 15 lines in the file (only the first 9 columns) it did not appear that there was any missing data. I would search for it, but since it is not obvious how it would be documented, I can't easily grep for it.
* The file is ASCII text, without very long lines.
* The most common chromosomes containing a SNP were chromosomes 1 (155 SNPs), 2 (127 SNPs) and 5 (122 SNPs). The most common genes containing SNPs were zmm28 (11 SNPs), PZA03450 (9 SNPs), and zag1 (8 SNPs). 

## Data Processing

### Maize Data

```
# Prepare maize directory files
mkdir -p ./maize/{interim,final} && mkdir ./maize/final/{increasing_pos_chromosomes,decreasing_pos_chromosomes}

# Extract maize samples from fang_et_al_genotypes.txt file and sorts for joining
awk 'NR == 1 || $3 ~ /ZMMIL|ZMMLR|ZMMMR/' fang_et_al_genotypes.txt | cut -f 4- | awk -f transpose.awk | sort -k1,1 > maize/interim/transposed_maize_genotypes.txt

# Extract the correct columns from the snp_position.txt file and sort for joining.
tail -n +2 snp_position.txt | sort -k3,3n | awk 'BEGIN {OFS="\t"} { print $1, $3, $4 }'  | sort -k1,1 > maize/interim/sorted_snp_position.txt

# Join the files together, allowing them to have unpairable lines
join -1 1 -2 1 -a 1 -a 2 maize/interim/sorted_snp_position.txt maize/interim/transposed_maize_genotypes.txt > maize/interim/joined_maize_snp_data.txt

# Sort based on increasing position values and generate 10 files
sort -k3,3n maize/interim/joined_maize_snp_data.txt | awk '$2 >= 1 && $2 <= 10 && $3 ~ /^[0-9]+$/ { print $0 > ("maize/final/increasing_pos_chromosomes/maize_chromosome_" $2 ".txt")}'

# Sort based on decreasing position values, replace ? with - for missing data, and 
# generate 10 chromosome files.
sort -k3,3nr maize/interim/joined_maize_snp_data.txt | sed -E 's/\?/-/g' | awk '$2 >= 1 && $2 <= 10 && $3 ~ /^[0-9]+$/ { print $0 > ("maize/final/decreasing_pos_chromosomes/maize_chromosome_" $2 ".txt")}'

# Generate unknown positions file
awk '$2 ~ /unknown/ || $3 ~ /unknown/  { print $0 >  "maize/final/maize_chromosome_unknown.txt"}' maize/interim/joined_maize_snp_data.txt

# Generate multiple positions file
awk '$2 ~ /multiple/ || $3 ~ /multiple/  { print $0 >  "maize/final/maize_chromosome_multiple.txt"}' maize/interim/joined_maize_snp_data.txt
```

- The code creates a folder for the processed maize files including a "final" folder where the final processed data files are stored as well as an "interim" folder where intermediate files that are needed for processing are stored.
- The code then extracts the data corresponding to maize samples from fang_et_al_genotypes.txt and transposes them with the provided awk script so that they can be joined with the SNP_ID, chromosome, and position columns that are extracted from snp_position.txt.
	- Note that the joined file is space delimitered, not tab delimitered so some commands are modified accordingly.
- Once the files are joined into the interim joined_maize_snp_data.txt file, the file is sorted and processed into individual files for each chromosome which are sorted based on increasing or decreasing position. These individual files are then saved in the proper subfolder "increasing_pos_chromosome" or "decreasing_pos_chromosome".
	- For the increasing position chromosome files, missing data is already coded by "?", so all that is needed is to sort by the third position column and extract the data corresponding to a certain chromosome (excluding both unknown and multiple files).
	- For the decreasing position chromosome files, after it is sorted the missing data is replaced by "-" globally using the sed command before being processed into the individual chromosome files.
- Once the chromosome files are generated, all the SNPs with unknown chromosomes (2nd column) or positions (3rd column) are added to the "maize_chromosome_unknown.txt" using awk. Similarly, all the SNPS with multiple chromosomes or positions are filed into the "maize_chromosome_multiple.txt" file.

### Teosinte Data

```
# Prepare teosinte directory files
mkdir -p ./teosinte/{interim,final} && mkdir ./teosinte/final/{increasing_pos_chromosomes,decreasing_pos_chromosomes}

# Extract teosinte samples from fang_et_al_genotypes.txt file and sorts for joining
awk 'NR == 1 || $3 ~ /ZMPBA|ZMPIL|ZMPJA/' fang_et_al_genotypes.txt | cut -f 4- | awk -f transpose.awk | sort -k1,1 > teosinte/interim/transposed_teosinte_genotypes.txt

# Extract the correct columns from the snp_position.txt file and sort for joining.
tail -n +2 snp_position.txt | sort -k3,3n | awk 'BEGIN {OFS="\t"} { print $1, $3, $4 }'  | sort -k1,1 > teosinte/interim/sorted_snp_position.txt

# Join the files together, allowing them to have unpairable lines
join -1 1 -2 1 -a 1 -a 2 teosinte/interim/sorted_snp_position.txt teosinte/interim/transposed_teosinte_genotypes.txt > teosinte/interim/joined_teosinte_snp_data.txt

# Sort based on increasing position values and generate 10 files
sort -k3,3n teosinte/interim/joined_teosinte_snp_data.txt | awk '$2 >= 1 && $2 <= 10 && $3 ~ /^[0-9]+$/ { print $0 > ("teosinte/final/increasing_pos_chromosomes/teosinte_chromosome_" $2 ".txt")}'

# Sort based on decreasing position values, replace ? with - for missing data, and 
# generate 10 chromosome files.
sort -k3,3nr teosinte/interim/joined_teosinte_snp_data.txt | sed -E 's/\?/-/g' | awk '$2 >= 1 && $2 <= 10 && $3 ~ /^[0-9]+$/ { print $0 > ("teosinte/final/decreasing_pos_chromosomes/teosinte_chromosome_" $2 ".txt")}'

# Generate unknown positions file
awk '$2 ~ /unknown/ || $3 ~ /unknown/  { print $0 >  "teosinte/final/teosinte_chromosome_unknown.txt"}' teosinte/interim/joined_teosinte_snp_data.txt

# Generate multiple positions file
awk '$2 ~ /multiple/ || $3 ~ /multiple/  { print $0 >  "teosinte/final/teosinte_chromosome_multiple.txt"}' teosinte/interim/joined_teosinte_snp_data.txt
```

- The code functions exactly as the maize processing script, with file names being modified to replace "maize" with "teosinte."