#!/bin/bash


# Help menu

Help()
{

## Display Help
   echo -ne "\nThe extract_amplicons.sh script allows extracting from a set of (reference) sequences the fragment(s) amplified by a primer set provided by the user. The mention \"_amplicon1\", \"_amplicon2\", etc. is added at the end of the amplimer sequence identifiers so that amplimers extracted from the same (reference) sequence can be distinguished."
   echo
   echo -ne "\nSyntax: extract_amplicons.sh [--fastafile|--fwdprimer|--revprimer|--mismatchperc|--minlength|--maxlength|--outdirectory|--help]"
   echo
   echo "Options:"
   echo -ne "\t--fastafile\tPath to the reference sequence fasta file.\n"
   echo -ne "\t--fwdprimer\tNucleotide sequence of the forward primer (5'->3'). Degenerated nucleotides are allowed.\n"
   echo -ne "\t--revprimer\tNucleotide sequence of the reverse primer (5'->3'). Degenerated nucleotides are allowed.\n"
   echo -ne "\t--mismatchperc\tMaximum number of mismatches allowed between primer sequence and homologous region in reference sequence. [Options: integer]\n"
   echo -ne "\t--minlength\tMinimum length allowed for the fragment amplified by the primers. [Options: integer]\n"
   echo -ne "\t--maxlength\tMaximum length allowed for the fragment amplified by the primers. [Options: integer]\n"
   echo -ne "\t--outdirectory\tPath to the output directory where extracted amplicons will be stored.\n"
   echo -ne "\t--help\t\tPrint this help."
   echo


## Usage example
   echo -ne "\nUsage example:"
   echo -ne "\n-------------"
   echo
   echo "extract_amplicons.sh \\"
   echo "  --fastafile input_file.fasta \\"
   echo "  --fwdprimer CGTTTGGCTCCACCACTACT \\"
   echo "  --revprimer ACGTTCAGTCTATCTTTCTATCACA \\"
   echo "  --mismatchperc 0 \\"
   echo "  --minlength 140 \\"
   echo "  --maxlength 150 \\"
   echo "  --outdirectory output_directory"
   echo
}

# Transform long options to short ones

for arg in "$@"; do
  shift
  case "$arg" in
    '--help')            set -- "$@" '-h'   ;;
    '--fastafile')       set -- "$@" '-a'   ;;
    '--fwdprimer')       set -- "$@" '-b'   ;;
    '--revprimer')       set -- "$@" '-c'   ;;
    '--mismatchperc')    set -- "$@" '-d'   ;;
    '--minlength')       set -- "$@" '-e'   ;;
    '--maxlength')       set -- "$@" '-f'   ;;
    '--outdirectory')    set -- "$@" '-g'   ;;
    *)                   set -- "$@" "$arg" ;;
  esac
done


# Parse short options

while getopts :ha:b:c:d:e:f:g: flag
do
    case "${flag}" in
        h) Help
        exit;;
        a) fasta_input=${OPTARG};;
        b) forward_primer=${OPTARG};;
        c) reverse_primer=${OPTARG};;
        d) mismatch_percent=${OPTARG};;
        e) minimum_length=${OPTARG};;
        f) maximum_length=${OPTARG};;
        g) output_directory=${OPTARG};;
        \?) # Invalid option
                echo "Error: Invalid option"
                exit;;
    esac
done


# Checking that the input fasta file does not contain in-seq linebreaks

awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" }END { printf "%s", n }' $fasta_input > $fasta_input\_tmp && mv $fasta_input\_tmp $fasta_input


# Creating the primer file

echo -ne "target1\t${forward_primer}\t${reverse_primer}\n" > primer_file


# Searching for putative sequences amplified by the primers

mkdir ${output_directory}

primersearch -seqall $fasta_input -infile primer_file -mismatchpercent $mismatch_percent -outfile ${output_directory}/output_file


# Reformatting files and keeping amplimers between minlength and maxlength

paste \
  <(awk 'NR >=3 {print $0}' ${output_directory}/output_file | awk '{printf "%s", $0} NR % 6 == 0 {print ""} END {print ""}' | cut -d $'\t' -f 2 | cut -d " " -f 2) \
  <(awk 'NR >=3 {print $0}' ${output_directory}/output_file | awk '{printf "%s", $0} NR % 6 == 0 {print ""} END {print ""}' | cut -d $'\t' -f 4 | awk -F "strand at " '{print $2}' | awk -F " with" '{print $1}') \
  <(awk 'NR >=3 {print $0}' ${output_directory}/output_file | awk '{printf "%s", $0} NR % 6 == 0 {print ""} END {print ""}' | cut -d $'\t' -f 5 | awk -F "strand at " '{print $2}' | awk -F " with" '{print $1}' | tr -d "[" | tr -d "]") \
  <(awk 'NR >=3 {print $0}' ${output_directory}/output_file | awk '{printf "%s", $0} NR % 6 == 0 {print ""} END {print ""}' | cut -d $'\t' -f 6 | cut -d ' ' -f 3) | awk -F '\t' -v min="$minimum_length" -v max="$maximum_length" '$4 >= min && $4 <= max {print}' > ${output_directory}/output_file2


# Adding corresponding sequences into the table

while IFS=$'\t' read -r id_seq reste_colonnes; do
    # Searching the sequence in the fasta file
    sequence=$(grep -A 1 "^>$id_seq" $fasta_input | tail -n 1)
    
    # Adding the sequence to the second column of the table
    echo -e "$id_seq\t$sequence\t$reste_colonnes"
done < ${output_directory}/output_file2 > ${output_directory}/output_file3   # Might be a time-consuming step


# Removing the column corresponding to amplicon length and replacing it by the entire sequence length

paste \
  <(awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,$4}' ${output_directory}/output_file3) \
  <(awk 'BEGIN {FS=OFS="\t"} {print length($2)}' ${output_directory}/output_file3) > ${output_directory}/output_file4


# Computing the reverse completement position for the reverse primer

awk 'BEGIN {FS=OFS="\t"} {print $1,$2,$3,$5-$4}' ${output_directory}/output_file4 > ${output_directory}/output_file5


# Extracting amplimers from each fasta file

while IFS=$'\t' read -r id_seq sequence start end; do
    # Extracting the amplicon according to the start and end positions
    amplicon=$(echo "$sequence" | cut -c "$start"-"$end")
    
    # Writing the amplimer sequence in the output fasta file
    echo ">$id_seq" >> ${output_directory}/extracted_amplicons.fasta
    echo "$amplicon" >> ${output_directory}/extracted_amplicons.fasta
done < ${output_directory}/output_file5

rm ${output_directory}/output*
rm primer_file


# Adding in the fasta deflines the amplicon number to be able to discriminate sequences with identical seq ids (ie amplicons extracted from the same parent sequence)

declare -A occurrences

# Parsing the fasta file
while IFS= read -r line; do
    if [[ $line =~ ^\>(.*)$ ]]; then
        id_seq="${BASH_REMATCH[1]}"
        # Checking if this identifier has already been encountered
        if [[ ${occurrences[$id_seq]+_} ]]; then
            # Increment the occurrence counter
            occurrences[$id_seq]=$((occurrences[$id_seq] + 1))
            # Add the suffix '_ampliconN' to the identifier
            id_seq="${id_seq}_amplicon${occurrences[$id_seq]}"
        else
            # The first occurrence of this identifier
            occurrences[$id_seq]=1
            # Adding the suffix '_amplicon1' to the identifier
            id_seq="${id_seq}_amplicon1"
        fi
        # Writing the modified identifier to the output file
        echo ">$id_seq" >> tmp_file
    else
        # Writing the sequence to the output file
        echo "$line" >> tmp_file
    fi
done < ${output_directory}/extracted_amplicons.fasta

mv tmp_file ${output_directory}/extracted_amplicons.fasta

unset fasta_input
unset forward_primer
unset reverse_primer
unset mismatch_percent
unset minimum_length
unset maximum_length
unset output_directory
unset occurences

