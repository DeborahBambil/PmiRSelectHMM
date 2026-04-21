#!/bin/bash
unrar x Stockholm.rar
unrar x data.rar

mkdir -p HMM  # garante que a pasta existe

for hmm_file in Stockholm/*; do
    hmm_name=$(basename "$hmm_file" .hmm)  # remove extensÃ£o .hmm
    hmm_name=$(basename "$hmm_name" .sto)  # remove extensÃ£o .sto se houver
    nhmmer \
        --tblout "HMM/${hmm_name}.tblout" \
        --cpu $(nproc) \
        "$hmm_file" \
        genome
done


# cria pasta com coordenadas extraÃ­das do .tblout
mkdir -p Coordenadas

for arquivo in HMM/*; do
    nome_arquivo=$(basename "$arquivo")
    # pega: seq_id, start, end
    awk '$1 != "#" {print $1, $9, $10}' "$arquivo" > "Coordenadas/${nome_arquivo%.*}.txt"
done

# remove arquivos vazios
find Coordenadas -size 0 -print -delete

mkdir -p PredictedHairpin

#tentativa1
mkdir -p PredictedHairpin

# Usamos o parallel para processar os arquivos de coordenadas
ls Coordenadas/ | parallel -j+0 "
    while read -r seq_id start end; do
        # Otimização CRÍTICA: o awk para de ler o arquivo assim que termina a sequência
        awk -v sid=\"\$seq_id\" -v s=\"\$start\" -v e=\"\$end\" '
            /^>/ { 
                if (found) exit; # Se já achamos e agora começou outro >, encerra o AWK
                if (\$0 ~ sid) {found=1; next} 
                else {found=0} 
            }
            found { seq = seq \$0 }
            END { 
                if (seq != \"\") {
                    print \">\" sid \":\" s \"-\" e; 
                    print substr(seq, s, e - s + 1) 
                }
            }
        ' genome >> PredictedHairpin/{}
    done < Coordenadas/{} "

# Agora ele vai chegar aqui!
echo "Limpando arquivos vazios..."
find PredictedHairpin -size 0 -print -delete

echo "Concluído com sucesso."


#COMANDO TODOS OS ARQUIVOS, ONDE FOR T troca por U

find PredictedHairpin/ -type f -exec sed -i '/^>/! s/T/U/g' {} +

#REMOVE>300
#Loop to iterate over the files in the 'PredictedHairpin' directory
for file in PredictedHairpin/*; do
    awk '/^>/ { if (seq != "" && length(seq) <= 301) { print header ORS seq } header = $0; seq = "" } !/^>/ { seq = seq $0 } END { if (length(seq) <= 301) { print header ORS seq } }' "$file" > arquivo_filtrado.fasta
    mv arquivo_filtrado.fasta "$file"
done

find PredictedHairpin -size  0 -print -delete

#remove .txt

cd PredictedHairpin && for f in *.txt; do mv "$f" "${f%.txt}"; done
cd ..

#BLASTNMATURETAB
mkdir ID
for file in PredictedHairpin/*; do
    filename=$(basename "$file")
    blastn -query "$file" -db "data/families/$filename" -out "ID/$filename" -evalue 0.01 -outfmt "6 qseqid" -word_size 15 &
done
wait
find ID -size  0 -print -delete

for arquivo in ID/*; do
  sort "$arquivo" | uniq > "$arquivo.tmp" && mv "$arquivo.tmp" "$arquivo"
done

#ID COM FA ">"
#Loop to iterate over the files in the "ID" directory.
for file in ID/*; do
    sed -i 's/^/>/' "$file"
done

#CompleteSequenceCurated
mkdir -p PredictedCurated

for idfile in ID/*; do
    fname=$(basename "$idfile")                # nome do arquivo
    hairpin_file="PredictedHairpin/$fname"     # fasta de entrada
    curated_file="PredictedCurated/$fname"     # saÃ­da
    
    # Converte cabeÃ§alhos em um "set"
    awk 'NR==FNR { if ($1 ~ /^>/) ids[$1]; next }
         /^>/ { keep = ($1 in ids) }
         keep' "$idfile" "$hairpin_file" > "$curated_file"
done


#Non-Redundant PredictedHairpin

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

find PredictedHairpin -size  0 -print -delete
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

#Delete gff
rm -f *.gff

#Non-Redundant PredictedCurated

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








