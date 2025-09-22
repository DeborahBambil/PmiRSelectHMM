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


# cria pasta com coordenadas extraídas do .tblout
mkdir -p Coordenadas

for arquivo in HMM/*; do
    nome_arquivo=$(basename "$arquivo")
    # pega: seq_id, start, end
    awk '$1 != "#" {print $1, $9, $10}' "$arquivo" > "Coordenadas/${nome_arquivo%.*}.txt"
done

# remove arquivos vazios
find Coordenadas -size 0 -print -delete

#remove profile
for f in Coordenadas/*profile*; do
    mv "$f" "${f/profile/}"
done

#remove profile
for f in Coordenadas/*.txt*; do
    mv "$f" "${f/.txt/}"
done


set -euo pipefail

GENOME="genome"   # <- ajuste aqui se seu FASTA se chama apenas "genome"
if [ ! -f "$GENOME" ]; then
    echo "Arquivo FASTA '$GENOME' não encontrado. Ajuste a variável GENOME." >&2
    exit 1
fi

mkdir -p PredictedHairpin

for arquivo_info in Coordenadas/*; do
    [ -f "$arquivo_info" ] || continue
    nome_arquivo_info=$(basename "$arquivo_info")
    out="PredictedHairpin/${nome_arquivo_info%.*}"
    : > "$out"   # zera/trunca arquivo de saída

    while read -r seq_id start end; do
        awk -v seq_id="$seq_id" -v start="$start" -v end="$end" '
        BEGIN {
            pattern = "^>" seq_id "([ \t]|$)"
            seq = ""
            id = 0
        }
        /^>/ {
            if (id) {
                # já lemos o contig desejado e chegamos a próximo header -> podemos parar
                exit
            }
            if ($0 ~ pattern) {
                id = 1
                seq = ""
                next
            }
        }
        id { seq = seq $0 }
        END {
            if (seq == "") {
                printf("WARNING: seq_id %s not found in %s\n", seq_id, ARGV[1]) > "/dev/stderr"
                exit
            }
            s = start + 0
            e = end + 0

            if (s <= e) {
                pos = s
                len = e - s + 1
            } else {
                pos = e
                len = s - e + 1
            }

            # valida e ajusta limites se necessário
            total_len = length(seq)
            if (pos < 1) {
                printf("WARNING: adjusted start <1 for %s:%s-%s\n", seq_id, start, end) > "/dev/stderr"
                pos = 1
            }
            if (pos + len - 1 > total_len) {
                printf("WARNING: adjusted end > seq length for %s:%s-%s (seq len=%d)\n", seq_id, start, end, total_len) > "/dev/stderr"
                len = total_len - pos + 1
            }
            if (len <= 0) {
                printf("WARNING: invalid region for %s:%s-%s (after adjustment)\n", seq_id, start, end) > "/dev/stderr"
                exit
            }

            subseq = substr(seq, pos, len)

            # se era reverse (start > end), faz reverse-complement
            if (s > e) {
                n = length(subseq)
                rc = ""
                for (i = n; i >= 1; i--) {
                    c = substr(subseq, i, 1)
                    if (c == "A") rc = rc "T"
                    else if (c == "a") rc = rc "t"
                    else if (c == "C") rc = rc "G"
                    else if (c == "c") rc = rc "g"
                    else if (c == "G") rc = rc "C"
                    else if (c == "g") rc = rc "c"
                    else if (c == "T") rc = rc "A"
                    else if (c == "t") rc = rc "a"
                    else if (c == "U" || c == "u") rc = rc "A"
                    else rc = rc c
                }
                subseq = rc
            }

            # imprime header e sequência com quebra a cada 60 chars
            printf(">%s:%s-%s\n", seq_id, start, end)
            L = length(subseq)
            for (i = 1; i <= L; i += 60) {
                print substr(subseq, i, 60)
            }
        }
        ' "$GENOME" >> "$out"
    done < "$arquivo_info"
done

#REMOVE>300
#Loop to iterate over the files in the 'PredictedHairpin' directory
for file in PredictedHairpin/*; do
    awk '/^>/ { if (seq != "" && length(seq) <= 301) { print header ORS seq } header = $0; seq = "" } !/^>/ { seq = seq $0 } END { if (length(seq) <= 301) { print header ORS seq } }' "$file" > arquivo_filtrado.fasta
    mv arquivo_filtrado.fasta "$file"
done

find PredictedHairpin -size  0 -print -delete

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
    curated_file="PredictedCurated/$fname"     # saída
    
    # Converte cabeçalhos em um "set"
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

DIR_COMPLETO="PredictedHairpin"

DIRETORIOS_INCOMPLETOS=(
    "NonIdentical" "Identical" "Redundant70" "NonRedundant70"
    "NonRedundant75" "Redundant75" "NonRedundant80" "Redundant80"
    "NonRedundant85" "Redundant85" "NonRedundant90" "Redundant90"
    "NonRedundant95" "Redundant95"
)

# Cria um arquivo temporário para armazenar os mapeamentos de cabeçalho
temp_file=$(mktemp)

echo "Criando mapa de cabeçalhos completos..."

# Percorre todos os arquivos no diretório completo e armazena o mapeamento
find "$DIR_COMPLETO" -name "*" -type f | while read -r arquivo; do
    # O 'awk' encontra o padrão ">" e extrai o cabeçalho incompleto e o completo
    awk '/^>/ {
        # Extrai a parte numérica (ex: 60049-60255) do cabeçalho
        # Usa o split para pegar a última parte se houver um ":"
        split($0, a, ":")
        incompleto = a[length(a)]
        
        # Remove o ">"
        sub(/^>/, "", incompleto)
        
        # Limpa espaços em branco e outros caracteres
        gsub(/[\r\t ]/, "", incompleto)
        
        # Mapeia a parte incompleta para o cabeçalho completo
        print incompleto, $0
    }' "$arquivo" >> "$temp_file"
done

# Percorre cada diretório incompleto
for dir_incompleto in "${DIRETORIOS_INCOMPLETOS[@]}"; do
    if [ ! -d "$dir_incompleto" ]; then
        echo "Aviso: Diretório $dir_incompleto não encontrado. Pulando..."
        continue
    fi
    
    echo "Processando diretório: $dir_incompleto"

    # Itera sobre cada arquivo incompleto
    find "$dir_incompleto" -name "*" -type f | while read -r arquivo; do
        # Usa 'awk' para fazer a substituição em cada arquivo
        awk '
            # Carrega o mapa de cabeçalhos
            BEGIN {
                while(getline < "'"$temp_file"'") {
                    map[$1] = $2
                }
            }
            
            # Se a linha começa com ">", processa
            /^>/ {
                # Extrai a parte numérica
                cabecalho_incompleto = $0
                sub(/^>/, "", cabecalho_incompleto)
                gsub(/[\r\t ]/, "", cabecalho_incompleto)

                # Se houver uma correspondência no mapa, substitui
                if (cabecalho_incompleto in map) {
                    print map[cabecalho_incompleto]
                    next
                }
            }
            # Se não for um cabeçalho, imprime a linha como está
            { print }
        ' "$arquivo" > "$arquivo.tmp" && mv "$arquivo.tmp" "$arquivo"
        
        echo "  - Arquivo atualizado: $arquivo"
    done
done

# Remove o arquivo temporário
rm "$temp_file"

echo "Processo concluído."


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

DIR_COMPLETO="PredictedCurated"

DIRETORIOS_INCOMPLETOS=(
    "CuratedNonIdentical" "CuratedIdentical" "CuratedRedundant70" "CuratedNonRedundant70"
    "CuratedNonRedundant75" "CuratedRedundant75" "CuratedNonRedundant80" "CuratedRedundant80"
    "CuratedNonRedundant85" "CuratedRedundant85" "CuratedNonRedundant90" "CuratedRedundant90"
    "CuratedNonRedundant95" "CuratedRedundant95"
)

# Cria um arquivo temporário para armazenar os mapeamentos de cabeçalho
temp_file=$(mktemp)

echo "Criando mapa de cabeçalhos completos..."

# Percorre todos os arquivos no diretório completo e armazena o mapeamento
find "$DIR_COMPLETO" -name "*" -type f | while read -r arquivo; do
    # O 'awk' encontra o padrão ">" e extrai o cabeçalho incompleto e o completo
    awk '/^>/ {
        # Extrai a parte numérica (ex: 60049-60255) do cabeçalho
        # Usa o split para pegar a última parte se houver um ":"
        split($0, a, ":")
        incompleto = a[length(a)]
        
        # Remove o ">"
        sub(/^>/, "", incompleto)
        
        # Limpa espaços em branco e outros caracteres
        gsub(/[\r\t ]/, "", incompleto)
        
        # Mapeia a parte incompleta para o cabeçalho completo
        print incompleto, $0
    }' "$arquivo" >> "$temp_file"
done

# Percorre cada diretório incompleto
for dir_incompleto in "${DIRETORIOS_INCOMPLETOS[@]}"; do
    if [ ! -d "$dir_incompleto" ]; then
        echo "Aviso: Diretório $dir_incompleto não encontrado. Pulando..."
        continue
    fi
    
    echo "Processando diretório: $dir_incompleto"

    # Itera sobre cada arquivo incompleto
    find "$dir_incompleto" -name "*" -type f | while read -r arquivo; do
        # Usa 'awk' para fazer a substituição em cada arquivo
        awk '
            # Carrega o mapa de cabeçalhos
            BEGIN {
                while(getline < "'"$temp_file"'") {
                    map[$1] = $2
                }
            }
            
            # Se a linha começa com ">", processa
            /^>/ {
                # Extrai a parte numérica
                cabecalho_incompleto = $0
                sub(/^>/, "", cabecalho_incompleto)
                gsub(/[\r\t ]/, "", cabecalho_incompleto)

                # Se houver uma correspondência no mapa, substitui
                if (cabecalho_incompleto in map) {
                    print map[cabecalho_incompleto]
                    next
                }
            }
            # Se não for um cabeçalho, imprime a linha como está
            { print }
        ' "$arquivo" > "$arquivo.tmp" && mv "$arquivo.tmp" "$arquivo"
        
        echo "  - Arquivo atualizado: $arquivo"
    done
done

# Remove o arquivo temporário
rm "$temp_file"

echo "Processo concluído."













