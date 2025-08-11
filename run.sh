#!/bin/bash
unrar x Stockholm.rar
unrar x data.rar

mkdir -p HMM  # garante que a pasta existe

for hmm_file in Stockholm/*; do
    hmm_name=$(basename "$hmm_file" .hmm)  # remove extensão .hmm
    hmm_name=$(basename "$hmm_name" .sto)  # remove extensão .sto se houver
    nhmmer \
        --tblout "HMM/${hmm_name}.tblout" \
        --cpu $(nproc) \
        "$hmm_file" \
        genome
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
        ' genome >> "PredictedHairpin/${nome_arquivo_info%.*}" &
    done < "$arquivo_info"
done

# Aguarda todas as tarefas em background terminarem
wait

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

#remover profile name
for arquivo in PredictedHairpin/*profile*; do
  novo_nome=${arquivo/profile/};
  mv "$arquivo" "$novo_nome";
done

mkdir ID
for file in PredictedHairpin/*; do
    filename=$(basename "$file")
    blastn -query "$file" -db "data/families/$filename" -out "ID/$filename" -evalue 0.01 -outfmt "6 qseqid" -word_size 15 &
done
wait

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

find PredictedCurated -size  0 -print -delete

#PredictedCurated
mkdir NonIdentical
mkdir Identical
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

# run skipredundant
executar_skipredundant() {
  local file=$1
  local threshold=$2
  local outseq_dir=$3
  local redundantoutseq_dir=$4

  filename=$(basename "$file")
  skipredundant -feature toggle -sequences "$file" [-datafile matrixf] -mode 1 -threshold "$threshold" -minthreshold 30 -maxthreshold 90 -gapopen 10 -gapextend 5 -outseq "$outseq_dir/$filename" -redundantoutseq "$redundantoutseq_dir/$filename"
}

#  comandos
thresholds=(100 70 75 80 85 90 95)
outseq_dirs=("NonIdentical" "NonRedundant70" "NonRedundant75" "NonRedundant80" "NonRedundant85" "NonRedundant90" "NonRedundant95")
redundantoutseq_dirs=("Identical" "Redundant70" "Redundant75" "Redundant80" "Redundant85" "Redundant90" "Redundant95")

# parallelo
for file in PredictedHairpin/*; do
  for ((i=0; i<${#thresholds[@]}; i++)); do
    executar_skipredundant "$file" "${thresholds[$i]}" "${outseq_dirs[$i]}" "${redundantoutseq_dirs[$i]}" &
  done
  wait
done

#Delete gff
rm -f *.gff


find NonIdentical -size  0 -print -delete
find Identical -size  0 -print -delete
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

# run skipredundant
executar_skipredundant() {
  local file=$1
  local threshold=$2
  local outseq_dir=$3
  local redundantoutseq_dir=$4

  filename=$(basename "$file")
  skipredundant -feature toggle -sequences "$file" [-datafile matrixf] -mode 1 -threshold "$threshold" -minthreshold 30 -maxthreshold 90 -gapopen 10 -gapextend 5 -outseq "$outseq_dir/$filename" -redundantoutseq "$redundantoutseq_dir/$filename"
}

# comandos
thresholds=(100 70 75 80 85 90 95)
outseq_dirs=("CuratedNonIdentical" "CuratedNonRedundant70" "CuratedNonRedundant75" "CuratedNonRedundant80" "CuratedNonRedundant85" "CuratedNonRedundant90" "CuratedNonRedundant95")
redundantoutseq_dirs=("CuratedIdentical" "CuratedRedundant70" "CuratedRedundant75" "CuratedRedundant80" "CuratedRedundant85" "CuratedRedundant90" "CuratedRedundant95")

# parallelo
for file in PredictedCurated/*; do
  for ((i=0; i<${#thresholds[@]}; i++)); do
    executar_skipredundant "$file" "${thresholds[$i]}" "${outseq_dirs[$i]}" "${redundantoutseq_dirs[$i]}" &
  done
  wait
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
