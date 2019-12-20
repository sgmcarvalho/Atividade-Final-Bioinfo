##**Execução do script rnaseq-denovo.sh**
##**1. Para a montagem denovo, sem genoma referência, esse script utiliza a ferramenta Trinity. O Trinity combina três etapas de montagem (Inchworm, Chrysalis, Butterflyé) para processar os dados, e particiona esses dados em diferentes Grafos de Bruijn para verificar independentemente a complexidade de cada locus, com as isoformas possíveis. Od Grafos de Bruijn correspondem a sequências cíclicas S de um alfabeto, de onde são derivadas as sequências de tamanho K-mers, consecutivas e que aparecem somente uma vez. E K-mers corresponde a um parâmetro crítico onde deve ser suficientemente grande para não pegas falsas sobreposições, devendo ser também, maior ou igual à especificidade.**

#São utilizados como input os arquivos gerados pelo prinseq, e estes serão renomeados (onde houver 1 e 2, será renomeado para R1 e R2), tanto para singletons como para sequências pareadas.
#Sequâncias com os sufixos 1.fastq, 2.fastq, 1singletons.fastq, 2singletons.fastq serão utilizadas como input.
#O comando para o script rnaseq-denovo.sh é 

```bash=
./rnaseq-denovo.sh ./output/processed/prinseq ./output
```
#E o script segue:
```bash=
# input - diretório contendo os arquivos de entrada no formato .fastq
input=$1
# validação do parâmetro "input"
if [ ! ${input} ]
then   
        echo "Missing input directory"
        exit
else   
        if [ ! -d ${input} ]
        then   
                echo "Wrong input directory ${input}"
                exit
        fi
fi
# output - diretório para armazenar o resultado do processo de montagem
output=$2
# validação do parâmetro "output"
if [ ! ${output} ]
then   
        echo "Missing output directory"
        exit
else   
        if [ ! -d ${output} ]
        then   
                echo "Wrong output directory ${output}"
                exit
        fi
fi
num_threads="8"
mem_gb="10G"

# Arquivos e diretórios de saída (output) 
basedir_out="${output}/"
renamed_out="${basedir_out}/renamed"
trinity_out="${basedir_out}/trinity_assembled"

mkdir -p ${renamed_out}

left=()
left_singleton=()
right=()
right_singleton=()

echo "Renaming step ..."
mkdir -p ${trinity_out}

for fastq in `ls ${input}/*.fastq`; do
	# obtendo nome do arquivo 
	fastqbn=`basename ${fastq}`;

	if [[ ! $fastqbn =~ \.bad_ ]]; then
		renamed_fastq="${renamed_out}/${fastqbn}"
		if [ ! -e ${renamed_fastq} ]; then
			echo -e "\tRenaming ${fastqbn} ..."
			if [[ ${fastqbn} =~ _1[\._] ]]; then
				awk '{ if (NR%4==1) { if ($1!~/\/1$/) { print $1"/1" } else { print $0 } } else if (NR%4==3) { print "+" } else { print $0 } }' ${fastq} > ${renamed_fastq}
			elif [[ ${fastqbn} =~ _2[\._]  ]]; then
				awk '{ if (NR%4==1) { if ($1!~/\/2$/) { print $1"/2" } else { print $0 } } else if (NR%4==3) { print "+" } else { print $0 } }' ${fastq} > ${renamed_fastq}
			else 
				echo "Warning: ${fastqbn} discarded!"
			fi
		fi
		if [[ ${fastqbn} =~ _1[\._] ]]; then
			if [[ ${fastqbn} =~ singletons ]]; then
				left_singleton=($(printf "%s\n" ${left_singleton[@]} ${renamed_fastq} | sort -u ))

			else
				left=($(printf "%s\n" ${left[@]} ${renamed_fastq}  | sort -u ))

			fi
		elif [[ ${fastqbn} =~ _2[\._] ]]; then
			if [[ ${fastqbn} =~ singleton ]]; then
				right_singleton=($(printf "%s\n" ${right_singleton[@]} ${renamed_fastq}  | sort -u ))

			else
				right=($(printf "%s\n" ${right[@]} ${renamed_fastq}  | sort -u ))
			fi

		else
			echo "Warning: ${fastqbn} discarded!"
		fi

	fi
done

if [ ! -d ${trinity_out}/Trinity.timing ]; then
	echo -e "Assembling step (Trinity) ..."

	rm -fr ${trinity_out}
	mkdir -p ${trinity_out}

	Trinity --KMER_SIZE 27 \
		--output ${trinity_out} \
		--seqType fq \
		--max_memory ${mem_gb} \
		--CPU ${num_threads} \
		--min_per_id_same_path 95 \
		--max_diffs_same_path  5 \
		--path_reinforcement_distance 5 \
		--group_pairs_distance 500 \
		--min_glue 5 \
		--min_contig_length 600 \
		--min_kmer_cov 3 \
		--left $(IFS=, ; echo "${left[*]},${left_singleton[*]}") \
		--right $(IFS=, ; echo "${right[*]},${right_singleton[*]}") \
		 > ${trinity_out}/Trinity.log.out.txt \
		2> ${trinity_out}/Trinity.log.err.txt
fi
```

#**2. Também existe a possibilidade do Trinity executar a montagem utilizando como guia os alinhamentos obtidos com STAR. P script é o rnaseq-ref-trinity.sh:** 
#O comando para o script rnaseq-ref-trinity.sh é: 

```bash=
./rnaseq-ref-trinity.sh ./output/star_out_final ./output
```
#E o script segue:
```bash=
# input - diretório contendo os arquivos de entrada no formato .fastq
input=$1

# as linhas que iniciam com cerquilha são comentários
# validação do parâmetro "input"

if [ ! ${input} ]
then   
        echo "Missing input directory"
        exit
else   
        if [ ! -d ${input} ]
        then   
                echo "Wrong input directory ${input}"
                exit
        fi
fi

# output - diretório para armazenar o resultado do processo de montagem
output=$2

# validação do parâmetro "output"
if [ ! ${output} ]
then   
        echo "Missing output directory"
        exit
else   
        if [ ! -d ${output} ]
        then   
                echo "Wrong output directory ${output}"
                exit
        fi
fi

num_threads="8"
mem_gb="16G"

# Arquivos e diretórios de entrada (input)
# Arquivos e diretórios de saída (output) 

basedir_out="${output}"
aligned_out="${basedir_out}/trinity_GG_input"

mkdir -p ${aligned_out}

trinity_out="${basedir_out}/trinity_GG_assembled"

# Criando diretórios para as saídas dos programas que serão utilizados a seguir
mkdir -p ${trinity_out}

if [ ! -e "${aligned_out}/All.sorted.bam" ]; then
	echo -e "Collecting alignments ..."

	bamfiles=()

	bamfiles=( $( find ${input} -name 'Aligned.out.sorted.bam' ) )

	samtools merge -f ${aligned_out}/All.sorted.bam ${bamfiles[*]}
	samtools sort --threads ${num_threads} ${aligned_out}/All.sorted.bam > ${aligned_out}/All.csorted.bam

	rm -f ${aligned_out}/All.sorted.bam
fi

if [ ! -d ${trinity_out}/Trinity.timing ]; then
	echo -e "Assembling step (Trinity) ..."
	rm -fr ${trinity_out}
	mkdir -p ${trinity_out}

	Trinity --KMER_SIZE 27 \
		--output ${trinity_out} \
		--seqType fq \
		--max_memory ${mem_gb} \
		--CPU ${num_threads} \
		--min_per_id_same_path 95 \
		--max_diffs_same_path  5 \
		--path_reinforcement_distance 5 \
		--group_pairs_distance 500 \
		--min_glue 5 \
		--min_contig_length 600 \
		--min_kmer_cov 3 \
		--genome_guided_bam ${aligned_out}/All.csorted.bam \
		--genome_guided_max_intron 10000 \
		 > ${trinity_out}/Trinity.log.out.txt \
		2> ${trinity_out}/Trinity.log.err.txt
        #--genome_guided_bam ${aligned_out}/All.sorted.bam \
		#--genome_guided_bam ./XM_024918093.1.bam \
fi
```

#Onde: --KMER_SIZE 27 especifica o tamanho da sequência K-mers, no caso 27pb; --output especifica os arquivos de saída gerados; --seqType fq especifica o tipo de de arquivo de reads, no caso arquivo FASTQC; --max_memory sugere o máximo de memória a ser usada pelo programa; --CPU para o número de CPUs a serem usados nas análises; --min_per_id_same_path 95 indica o mínimo de percentual de identidade para dois caminhos serem fundidos em um só, no caso 95%; --max_diffs_same_path  5 para o máximo de diferença entre dois caminhos, no caso 5%; --path_reinforcement_distance 5 para reforçar o máximo de diferença; --group_pairs_distance 500 para o máximo de comprimento esperado entre os pares de fragmentos, no caso 500pb; --min_glue 5 para o número mínimo de reads necessários para unir dois contigs, no caso 5; --min_contig_length 600 comprimento mínimo do contig montado para ser relatado, no caso 600pb; --min_kmer_cov 3 para contagem mínima para o K-mers ser montado, no caso 5; --left nomeia os arquivos de leituras à esquerda ; --right nomeia os arquivos de leituras à direita.

##**3. Para avaliação das montagens, comparação e montagem com o Trinity Genome Guided, segue:**

```bash=
makeblastdb -title "RefTrans" \
			-dbtype nucl \
			-out RefTrans \
			-in transcriptoma.fa \
			 > makebleastdb.log.out.txt \
			2> makeblastdb.log.err.txt

blastn 	-max_hsps 1 \
	-max_target_seqs 1 \
	-num_threads 8 \
	-query ./output/trinity_assembled/Trinity.fasta \
	-task blastn \
	-db ./RefTrans \
	-out ./Trinity_x_RefTrans.txt \
	-evalue 1e-5 \
	-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
	 > ./Trinity_x_RefTrans.log.out.txt \
	2> ./Trinity_x_RefTrans.log.err.txt
```

#Onde: no makeblastdb, -title entitula a base de dados BLAST como "RefTrans"; -dbtype nucl é um argumento requerido para o tipo de molécula alvo; -out RefTrans especifica o nome da base de dados a ser criada, no caso RefTrans; -in especifica o input, no caso o arquivo transcriptoma.fa.
      no blastn, -max_hsps 1 define o número máximo de HSPs por sequência de subject para salvar para cada query; -max_target_seqs 1 número máximo de sequências alinhadas para manter;	-num_threads 8 indica o número de CPUs para usar nas análises do blast; -query indica o diretório e nome do arquivo de input, no caso ./output/trinity_assembled/Trinity.fasta;	-task para executar o blastn; -db para especificar o nome da base de dados, no cacso ./RefTrans ; -out indica o nome do arquivo de saída, no caso ./Trinity_x_RefTrans.txt;	-evalue 1e-5 especifica o valor de threshold esperado para o E-value, no caso 5; -outfmt 6 para alinhamento visualizado de forma tabular.

##**4. Para visualização das sequências query (qseqid) e subject (sseqid) HSPs de todos HITs obtidos, executa-se:**

```bash=
cut -f 1,2 ./Trinity_x_RefTrans.txt
cut -f 2 ./Trinity_x_RefTrans.txt | sort | uniq -c

cut -f 1,2 ./Trinity-GG_x_RefTrans.txt
cut -f 2 ./Trinity-GG_x_RefTrans.txt | sort | uniq -c
```

#E dessa forma selecionam-se os nomes das query e subject, e as ocorrências de sequências da base de dados (subject).

##**5. A etapa seguinte realiza a estimação da abundância e encontra as isoformas ou genes com expressão diferencial. Nesse caso, utiliza-se o ambiente ${TRINITY_HOME} e verificam-se os scripts auxiliares align_and_estimate_abundance.pl e abundance_estimates_to_matrix.pl:**

```bash=
${TRINITY_HOME}/util/align_and_estimate_abundance.pl

${TRINITY_HOME}/util/abundance_estimates_to_matrix.pl
```

#Também foram executados os scripts run-DESeq2.R e rnaseq-trinity-abundance.sh, com o comando:
```bash=
rnaseq-trinity-abundance.sh	./output/renamed/ \
				./output/trinity_assembled/ \
				./output
```

#Seguindo o script:
```bash=
# input - diretório contendo os arquivos de entrada no formato .fastq
input=$1

# validação do parâmetro "input"
if [ ! ${input} ]
then   
        echo "Missing input (renamed for Trinity) directory"
        exit
else   
        if [ ! -d ${input} ]
        then   
                echo "Wrong input (renamed for Trinity) directory ${input}"
                exit
        fi
fi

# trinity_output - diretório onde foi armazenado o resultado do processo de montagem
trinity_output=$2

# validação do parâmetro "trinity_output"
if [ ! ${trinity_output} ]
then   
        echo "Missing Trinity output directory"
        exit
else   
        if [ ! -d ${trinity_output} ]
        then   
                echo "Wrong Trinity output directory ${trinity_output}"
                exit
        fi
fi

# output - diretório para armazenar o resultado da avaliação de abundância
output=$3

# validação do parâmetro "output"
if [ ! ${output} ]
then   
        echo "Missing output directory"
        exit
else   
        if [ ! -d ${output} ]
        then   
                echo "Wrong output directory ${output}"
                exit
        fi
fi

num_threads="8"
# Arquivos e diretórios de saída (output) 

abundance_out="${output}/abundance"
mkdir -p ${abundance_out}

left=()
right=()

echo "Collecting reads step ..."

left=($(find ${input} -type f -name '*_1.fastq'))

rm -f ${abundance_out}/samples.txt
rm -f ${abundance_out}/quant_files.txt
rm -f ${abundance_out}/groups.txt

echo -e "id\tname\tgroup" > ${abundance_out}/groups.txt

for l in ${left[@]}; do

	repname=`basename ${l} | sed 's/\..*$//'`
	condname=`echo ${repname} | sed 's/[0-9]\+//'`
	r=`echo ${l} | sed 's/_1.fastq/_2.fastq/'`
	right=(${right[@]} ${r})

	echo -e "${condname}\t${abundance_out}/${repname}\t${l}\t${r}" >> ${abundance_out}/samples.txt
	echo -e "${abundance_out}/${repname}/quant.sf" >> ${abundance_out}/quant_files.txt

	echo -e "${repname}\t${repname}\t${condname}" >> ${abundance_out}/groups.txt
done

#echo ${left[*]}
#echo ${right[*]}

trinity_fasta=`find ${trinity_output} -type f -name 'Trinity*.fasta'`
trinity_trans_map=`find ${trinity_output} -type f -name Trinity*.gene_trans_map`

echo "Estimating abundances ..."

${TRINITY_HOME}/util/align_and_estimate_abundance.pl 	--transcripts	${trinity_fasta} \
							--est_method	salmon \
							--salmon_add_opts "--validateMappings" \
							--samples_file	${abundance_out}/samples.txt \
							--gene_trans_map ${trinity_trans_map} \
							--prep_reference \
							--thread_count ${num_threads} \
							--seqType fq \
							--output_dir ${abundance_out} \
							 > ${abundance_out}/align_and_estimate_abundance.log.out.txt \
							2> ${abundance_out}/align_and_estimate_abundance.log.err.txt

echo "Constructing abundance matrix ..."

${TRINITY_HOME}/util/abundance_estimates_to_matrix.pl	--est_method salmon \
							--gene_trans_map ${trinity_trans_map} \
							--name_sample_by_basedir \
							--cross_sample_norm none \
							--quant_files ${abundance_out}/quant_files.txt \
							--out_prefix ${abundance_out}/abundance \
							 > ${abundance_out}/abundance_estimates_to_matrix.log.out.txt \
							2> ${abundance_out}/abundance_estimates_to_matrix.log.err.txt 

echo "Calculating Differentially Expressed Genes ..."
mkdir -p ${abundance_out}/DEG

run-DESeq2.R 	--in="${abundance_out}/abundance.gene.counts.matrix"  \
		--groups="${abundance_out}/groups.txt" \
		--out="${abundance_out}/DEG" \
		 > ${abundance_out}/DEG/run-DESeq2.log.out.txt \
		2> ${abundance_out}/DEG/run-DESeq2.log.err.txt
```

##**6. E Para visualização dos resultados seguiu-se:**
```bash=
cat ./output/abundance/DEG/SAMPLEA-SAMPLEB.txt | column -s$'\t' -t | less -S
```

#Ou pelo script do R-Studio:
```bash=
# Carregando a biblioteca xlsx para lidar com arquivos do Excel
library("xlsx")
# Carregando a biblioteca da função para gerar o HeatMap - heatmap.2()
library("gplots")
# Carregando a biblioteca com a função para criar a paleta de cores
library("RColorBrewer")

# Equivalente ao comando "pwd" no linux
getwd()
# Equivalente ao comando "cd" no linux
setwd("/state/partition1/scarvalho")

getwd()

# Carrega arquivo texto separado por TAB ("\t") com cabeçalho e sem transformar 
# strings em fatores
deg.df <- read.delim(file="./output/abundance/DEG/SAMPLEA-SAMPLEB.txt", 
                     sep="\t", 
                     header=TRUE, 
                     stringsAsFactors = FALSE)

# exibir somente as 2 primeiras linhas para checagem
head(deg.df,2)

# Dimensão do data.frame, onde o primeiro valor refere-se à quantidade de linhas e
# o segundo valor à quantidade de colunas
dim(deg.df)

# selecionar um subconjunto de linhas do data.frame deg.df em um outro data.frame
# ("subset.deg.df") contendo somente os genes que possuírem logFC >= 2 ou logFC <= -2,
# além de possuírem um p-valor corrigido ("FDR") <= 0.05
subset.deg.df <- subset(deg.df, (((logFC >= 2) | (logFC <= -2)  ) & (FDR <=0.05))  )

dim(subset.deg.df)

head(subset.deg.df,2)

# criando um vetor de nomes de colunas selecionadas,
# ou seja, colunas que não são "X", "logFC", "PValue" ou "FDR",
# as quais, sabemos previamente que contém os valores de expressão.
# Os nomes das colunas serão ordenados depois de secionados com a 
# função "setdiff", a qual obtém
# um conjunto de valores (vetor) coma a diferença entre dois conjuntos de
# valores: o de todas as colunas ("colnames(subset.deg.df)") menos os nomes
# das colunas indesejadas.

sel_columns <- sort(
                      setdiff( colnames(subset.deg.df), 
                        c("X", "logFC", "PValue", "FDR") 
                      )
                    )

# criando uma nova matriz "expression_data" que contém somente as colunas
# selecionadas
expression_data <- as.matrix( subset.deg.df[,sel_columns] )

# Atribuindo para os nomes das linhas da matriz o conteúdo da coluna "X",
# onde sabemos previamente que contém o identificador dos genes
rownames(expression_data) <- subset.deg.df$X

head(expression_data,2)

# Função para gravar o data.frame em um arquivo Excel, na planilha 
# de nome "DEGS_SAMPLEA-SAMPLEB"
write.xlsx(subset.deg.df,  
           file="./output/abundance/DEG/SAMPLEA-SAMPLEB.xlsx",
           sheetName="DEGS_SAMPLEA-SAMPLEB"
           )


# cria uma paleta personalizada de 299 cores do vermelho ao verde,
# passando pelo amarelo
my_palette <- colorRampPalette(c("red", "yellow", "green"))(n = 299)

# define as quebras das cores manualmente para transição de cores
col_breaks = c(seq(-1,0,length=100),        # for red
               seq(0.01,0.8,length=100),    # for yellow
               seq(0.81,1,length=100))      # for green


# criação de uma imagem de tamanho 5 x 5 polegadas
png("./output/abundance/DEG/heatmap_DEGS_SAMPLEA-SAMPLEB.png",    # cria arquivo do tipo PNG for the heat map        
    width = 5*300,        # 5 x 300 pixels de largura
    height = 5*300,       # 5 x 300 pixels de altura
    res = 300,            # 300 pixels por polegada
    pointsize = 8)        # tamanho da fonte

heatmap.2(expression_data, # a matriz com os valores de expressão
          main = "Correlation", # Título do HeatMap
          density.info="none",  # desabilita o gráfico de densidade dentro da legenda
          trace="none",         # desabilita as linhas dentro do HeatMap
          margins =c(12,12),     # definiação das margens no entorno do gráfico
          col=my_palette,       # nome do objeto contendo a paleta de cores criada anteriormente
          breaks=col_breaks,    # pontos de quebra para a transição de cores
          dendrogram="both",    # desenhar dendrograma para linhas e colunas
          distfun = function(x) as.dist(1-cor(t(x))), # distância baseada em correlação
          hclustfun = function(x) hclust(x, method="centroid") # método de ligação pelo centróide
          )

dev.off()               # fecha o arquivo da imagem PNG
```
