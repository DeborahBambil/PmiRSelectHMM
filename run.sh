#!/bin/bash
unrar x Stockholm.rar
unrar x data.rar

mkdir HMM

for hmm_file in Stockholm/*; do
    hmm_name=$(basename "$hmm_file")
    hmmsearch --domtblout "HMM/${hmm_name}" --cpu $(nproc) "$hmm_file" genome
done

mkdir Coordenadas

for arquivo in HMM/*; do
    nome_arquivo=$(basename "$arquivo")
    awk '$1 != "#" {print $1, $20, $21}' "$arquivo" > "Coordenadas/${nome_arquivo%.*}.txt"
done

find Coordenadas -size  0 -print -delete

mkdir PredictedHairpin

for arquivo_info in Coordenadas/*; do
    nome_arquivo_info=$(basename "$arquivo_info")
    while read -r seq_id start end; do
        awk -v seq_id="$seq_id" -v start="$start" -v end="$end" '
            /^>/ {if (id) { print ">" seq_id ":" start "-" end "\n" substr(seq, start, end - start + 1); seq=""; } id = ($0 ~ seq_id); next;}
            id {seq = seq $0}
            END {if (id) print ">" seq_id ":" start "-" end "\n" substr(seq, start, end - start + 1);}
        ' genome >> "PredictedHairpin/${nome_arquivo_info%.*}"
    done < "$arquivo_info"
    
done

find PredictedHairpin -size  0 -print -delete

#REMOVE>300
#Loop to iterate over the files in the 'PredictedHairpin' directory
for file in PredictedHairpin/*; do
    awk '/^>/ { if (seq != "" && length(seq) <= 301) { print header ORS seq } header = $0; seq = "" } !/^>/ { seq = seq $0 } END { if (length(seq) <= 301) { print header ORS seq } }' "$file" > arquivo_filtrado.fasta
    mv arquivo_filtrado.fasta "$file"
done

find PredictedHairpin -size  0 -print -delete

for file in PredictedHairpin/*profile.hmm_saida_domtblout; do
    newname=$(echo "$file" | sed 's/profile\.hmm_saida_domtblout//')
    mv "$file" "$newname"
done

sed -i 's/:/\//' PredictedHairpin/*


#BLASTNMATURETAB
mkdir ID
for file in PredictedHairpin/*; do
    filename=$(basename "$file")
     blastn -query "$file" -db "data/families/$filename" -out "ID/$filename" -evalue 0.01 -outfmt "6 qseqid" -word_size 15
done

find ID -size  0 -print -delete

#ID COM FA ">"
#Loop to iterate over the files in the "ID" directory.
for file in ID/*; do
    sed -i 's/^/>/' "$file"
done

mkdir PredictedCurated
#Loop to iterate over the files in the "PredictedHairpin" directory.
for file in PredictedHairpin/*; do
    filename=$(basename "$file")
    grep -x -F -A 1000 -f "ID/$filename" "$file" > "PredictedCurated/$filename"
done

input_dir="PredictedHairpin"  # Substitua pelo caminho do diretório contendo os arquivos

for file in "$input_dir"/*; do
    if [ -f "$file" ]; then
        sed -i -e '/^>/!y/T/U/' "$file"
    fi
done

input_dir="PredictedHairpin"  # Substitua pelo caminho do diretório contendo os arquivos

for file in "$input_dir"/*; do
    if [ -f "$file" ]; then
        sed -i -e '/^>/!y/t/u/' "$file"
    fi
done


input_dir="PredictedCurated"  # Substitua pelo caminho do diretório contendo os arquivos

for file in "$input_dir"/*; do
    if [ -f "$file" ]; then
        sed -i -e '/^>/!y/T/U/' "$file"
    fi
done


input_dir="PredictedCurated"  # Substitua pelo caminho do diretório contendo os arquivos

for file in "$input_dir"/*; do
    if [ -f "$file" ]; then
        sed -i -e '/^>/!y/t/u/' "$file"
    fi
done

find PredictedCurated -size  0 -print -delete

mkdir PredictedNonIdentical
mkdir PredictedIdentical
mkdir Redundant70
mkdir NonRedundant70
mkdir NonRedundant75
mkdir Redundant75
mkdir NonRedundant80
mkdir Redundant80
mkdir NonRedundant85
mkdir Redundant85
mkdir NonRedundant90
mkdir Redundant90
mkdir NonRedundant95
mkdir Redundant95

#RedundanceIdenticalNonIdentical
for file in PredictedHairpin/*; do
    filename=$(basename "$file")
    skipredundant -feature toggle -sequences "$file" [-datafile matrixf] -mode 1 -threshold 100 -minthreshold 30 -maxthreshold 90 -gapopen 10 -gapextend 5 -outseq "PredictedNonIdentical/$filename" -redundantoutseq "PredictedIdentical/$filename"
done

#REDUNDANCE70
for file in PredictedHairpin/*; do
    filename=$(basename "$file")
    skipredundant -feature toggle -sequences "$file" [-datafile  matrixf] -mode 1 -threshold 70 -minthreshold 30 -maxthreshold 90 -gapopen  10 -gapextend 5 -outseq "NonRedundant70/$filename" -redundantoutseq "Redundant70/$filename"
done

#REDUNDANCE75
for file in PredictedHairpin/*; do
    filename=$(basename "$file")
    skipredundant -feature toggle -sequences "$file" [-datafile  matrixf] -mode 1 -threshold 75 -minthreshold 30 -maxthreshold 90 -gapopen  10 -gapextend 5 -outseq "NonRedundant75/$filename" -redundantoutseq "Redundant75/$filename"
done

#REDUNDANCE80
for file in PredictedHairpin/*; do
    filename=$(basename "$file")
    skipredundant -feature toggle -sequences "$file" [-datafile  matrixf] -mode 1 -threshold 80 -minthreshold 30 -maxthreshold 90 -gapopen  10 -gapextend 5 -outseq "NonRedundant80/$filename" -redundantoutseq "Redundant80/$filename"
done

#REDUNDANCE85
for file in PredictedHairpin/*; do
    filename=$(basename "$file")
    skipredundant -feature toggle -sequences "$file" [-datafile  matrixf] -mode 1 -threshold 85 -minthreshold 30 -maxthreshold 90 -gapopen  10 -gapextend 5 -outseq "NonRedundant85/$filename" -redundantoutseq "Redundant85/$filename"
done

#REDUNDANCE90
for file in PredictedHairpin/*; do
    filename=$(basename "$file")
    skipredundant -feature toggle -sequences "$file" [-datafile  matrixf] -mode 1 -threshold 90 -minthreshold 30 -maxthreshold 90 -gapopen  10 -gapextend 5 -outseq "NonRedundant90/$filename" -redundantoutseq "Redundant90/$filename"
done

#REDUNDANCE95
for file in PredictedHairpin/*; do
    filename=$(basename "$file")
    skipredundant -feature toggle -sequences "$file" [-datafile  matrixf] -mode 1 -threshold 95 -minthreshold 30 -maxthreshold 90 -gapopen  10 -gapextend 5 -outseq "NonRedundant95/$filename" -redundantoutseq "Redundant95/$filename"
done

#Delete gff
rm -f *.gff

find PredictedNonIdentical -size  0 -print -delete
find PredictedIdentical -size  0 -print -delete
find Redundant70 -size  0 -print -delete
find NonRedundant70 -size  0 -print -delete
find NonRedundant75 -size  0 -print -delete
find Redundant75 -size  0 -print -delete
find NonRedundant80 -size  0 -print -delete
find Redundant80 -size  0 -print -delete
find NonRedundant85 -size  0 -print -delete
find Redundant85 -size  0 -print -delete
find NonRedundant90 -size  0 -print -delete
find Redundant90 -size  0 -print -delete
find NonRedundant95 -size  0 -print -delete
find Redundant95 -size  0 -print -delete

#PredictedCurated
mkdir CuratedNonIdentical
mkdir CuratedIdentical
mkdir CuratedNonRedundant70
mkdir CuratedRedundant70
mkdir CuratedNonRedundant75
mkdir CuratedRedundant75
mkdir CuratedNonRedundant80
mkdir CuratedRedundant80
mkdir CuratedNonRedundant85
mkdir CuratedRedundant85
mkdir CuratedNonRedundant90
mkdir CuratedRedundant90
mkdir CuratedNonRedundant95
mkdir CuratedRedundant95


#CompareMatureRedundanceIdenticalNonIdentical
for file in PredictedCurated/*; do
    filename=$(basename "$file")
    skipredundant -feature toggle -sequences "$file" [-datafile matrixf] -mode 1 -threshold 100 -minthreshold 30 -maxthreshold 90 -gapopen 10 -gapextend 5 -outseq "CuratedNonIdentical/$filename" -redundantoutseq "CuratedIdentical/$filename"
done

#CompareMatureREDUNDANCE70
for file in PredictedCurated/*; do
    filename=$(basename "$file")
    skipredundant -feature toggle -sequences "$file" [-datafile matrixf] -mode 1 -threshold 70 -minthreshold 30 -maxthreshold 90 -gapopen 10 -gapextend 5 -outseq "CuratedNonRedundant70/$filename" -redundantoutseq "CuratedRedundant70/$filename"
done

#CompareMatureREDUNDANCE75
for file in PredictedCurated/*; do
    filename=$(basename "$file")
    skipredundant -feature toggle -sequences "$file" [-datafile matrixf] -mode 1 -threshold 75 -minthreshold 30 -maxthreshold 90 -gapopen 10 -gapextend 5 -outseq "CuratedNonRedundant75/$filename" -redundantoutseq "CuratedRedundant75/$filename"
done

#CompareMatureREDUNDANCE80
for file in PredictedCurated/*; do
    filename=$(basename "$file")
    skipredundant -feature toggle -sequences "$file" [-datafile matrixf] -mode 1 -threshold 80 -minthreshold 30 -maxthreshold 90 -gapopen 10 -gapextend 5 -outseq "CuratedNonRedundant80/$filename" -redundantoutseq "CuratedRedundant80/$filename"
done

#CompareMatureREDUNDANCE85

for file in PredictedCurated/*; do
    filename=$(basename "$file")
    skipredundant -feature toggle -sequences "$file" [-datafile matrixf] -mode 1 -threshold 85 -minthreshold 30 -maxthreshold 90 -gapopen 10 -gapextend 5 -outseq "CuratedNonRedundant85/$filename" -redundantoutseq "CuratedRedundant85/$filename"
done

#CompareMatureREDUNDANCE90
for file in PredictedCurated/*; do
    filename=$(basename "$file")
    skipredundant -feature toggle -sequences "$file" [-datafile matrixf] -mode 1 -threshold 90 -minthreshold 30 -maxthreshold 90 -gapopen 10 -gapextend 5 -outseq "CuratedNonRedundant90/$filename" -redundantoutseq "CuratedRedundant90/$filename"
done

#CompareMatureREDUNDANCE95
for file in PredictedCurated/*; do
    filename=$(basename "$file")
    skipredundant -feature toggle -sequences "$file" [-datafile matrixf] -mode 1 -threshold 95 -minthreshold 30 -maxthreshold 90 -gapopen 10 -gapextend 5 -outseq "CuratedNonRedundant95/$filename" -redundantoutseq "CuratedRedundant95/$filename"
done


find PredictedCurated -size  0 -print -delete
find CuratedNonIdentical -size  0 -print -delete
find CuratedIdentical -size  0 -print -delete
find CuratedRedundant70 -size  0 -print -delete
find CuratedNonRedundant70 -size  0 -print -delete
find CuratedNonRedundant75 -size  0 -print -delete
find CuratedRedundant75 -size  0 -print -delete
find CuratedNonRedundant80 -size  0 -print -delete
find CuratedRedundant80 -size  0 -print -delete
find CuratedNonRedundant85 -size  0 -print -delete
find CuratedRedundant85 -size  0 -print -delete
find CuratedNonRedundant90 -size  0 -print -delete
find CuratedRedundant90 -size  0 -print -delete
find CuratedNonRedundant95 -size  0 -print -delete
find CuratedRedundant95 -size  0 -print -delete
#deletefiles

#Delete gff
rm -f *.gff
