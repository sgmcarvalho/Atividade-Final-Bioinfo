#**Atividade final - Relatório**
#**Relatório com o script comentado rnaseq-ref.sh e montagem denovo**
#Escrever os materiais e métodos necessários para a análise de expressão gênica diferencial. 

#Para o cufflinks alterar os parâmetros de tamanho dos introns etc, pois o programa é bem específico para os arquivos gffs
#Na conferência dos arquivos gerados no prinseq, deve constar que eles apresentam conteúdo e tamanho 
#Comando ls -lh ./output/processed/prinseq/*_1.fastq 

#Descrever os resultados obtidos - nº de sequências submetidas, nº de sequências processadas, nº de sequências alinhadas, genes diferencialmente expressos, adicionar informações aos genes.

## **Análises preliminares - Geração de arquivos necessários para o rnaseq-ref.sh**
#**1. Anteriormente à execução do script rnaseq-ref.sh, foram gerados os seguintes arquivos:**

#1. abundance_A.txt; 2. abundance_B.txt; 3. transcriptoma.fa; 4. ACCS.txt

#Esses arquivos representam a simulação dos genes 'sequenciados e identificados', as suas respectivas abundâncias, em duas condições biológicas (A e B), e um arquivo contendo apenas os números de acesso dos genes escolhidos para ser utilizado para gerar o arquivo transcriptoma.fa. 

#Para a atual atividade, foram selecionados 5 genes de *Arabdopsis thaliana*: LEAFY (LFY), APETALA1 (AP1), PISTILLATA (PI), AGAMOUS (AG), SEEDSTICK (STK) 
#Os genes escolhidos pertencem principalmente à família gênica MADS-box, que codificam fatores de transcrição, intimamente ligados à identidade e desenvolvimento de órgãos florais em Angiospermas.

#A execução dessa etapa foi feita utilizando o pacote de ferramentas E-Direct e a linha de comando:

```bash=
rm -f transcriptoma.fa
for acc in $(cat ACCS.txt); do
        echo "Pegando FASTA para ${acc} ..."

        esearch -db nucleotide -query ${acc} | efetch \
        -format fasta >> transcriptoma.fa
done
```

#Dessa forma, foi gerado o arquivo em formato fasta transcriptoma.fa, contendo as sequências dos 5 genes descritos anteriormente.
#Os arquivos de abundância foram gerados manualmente, utilizando os editores de texto "nano" ou "vim".

#Para gerar o arquivo ACCS.txt foram utilizados os seguintes comandos e em seguida gerado o arquivo transcriptoma.fa de forma mais rápida: 

```bash=
rm -f transcriptoma.fa
for acc in XM_024918093.1 XM_024914151.1 XM_024918802.1 XM_024922757.1 XM_024914431.1 XM_024914432.1 ; do
        echo "Pegando FASTA para ${acc} ..."
        esearch -db nucleotide -query ${acc} | efetch \
        -format fasta >> transcriptoma.fa
done

rm -f transcriptoma.fa
for acc in $(cat ACCS.txt); do
        echo "Pegando FASTA para ${acc} ..."
        esearch -db nucleotide -query ${acc} | efetch \
        -format fasta >> transcriptoma.fa
done
```

```bash=
rm -f transcriptoma.fa

IFS=$'\n'
for accline in $(cat ./ACCS.txt); do
        acc=`echo ${accline} | cut -f 1`
        seqref=`echo ${accline} | cut -f 2`
        chr_start=`echo ${accline} | cut -f 3`
        chr_stop=`echo ${accline} | cut -f 4`
        strand=`echo ${accline} | cut -f 5`

        echo "Pegando FASTA para ${acc}  [${seqref}:${chr_start}-${chr_stop}(${strand})] ..."
        efetch -db nucleotide -id ${seqref}  -format fasta \
        -chr_start ${chr_start} \
        -chr_stop ${chr_stop} \
        -strand ${strand} | \
        sed "s/^>.*/>${acc}/" \
        >> transcriptoma.fa
done
```

#**2. A próxima fase da primeira etapa foi gerar arquivos .fastq que representam as réplicas biológicas com os fragmentos teoricamente gerados no sequenciamento, considerando as respectivas abundâncias dos genes:**

#Utilizou-se o comando generate_fragments.py

```bash
generate_fragments.py -r transcriptoma.fa \
   -a ./abundance_A.txt \
   -o ./tmp.frags.r1 \
   -t 25000 \
   -i 300 \
   -s 30
```

#Onde: -r indica o arquivo fasta itulizado para gerar os fragmentos; -a indica o arquivo de abundância a ser utilizado; -o indica o arquivo para saída dos fragmentos gerados; -t indica o número total de reads gerados, no caso, 25000; -i indica o tamanho médio dos fragmentos; -s indica o desvio padrão para o comprimento de inserção, no caso 30pb.

#**3. O processo continua com o renomeamento do arquivo, simulação *paired* de leituras, desintercalamento das leituras em arquivos R1 e R2, remoção e aruivos desnecessários e verificação das sequências geradas.** 

```bash=
cat ./tmp.frags.r1.1.fasta | renameSeqs.pl \
   -if FASTA \
   -of FASTA \
   -p SAMPLEA1 \
   -w 1000 | \
   sed 's/^>\(\S\+\).*/>\1/' \
   > ./frags_A1.fa
   
cat ./frags_A1.fa | simNGS -a \
AGATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
-p paired \
/usr/local/bioinfo/simNGS/data/s_4_0099.runfile \
-n 151 > ./SAMPLEA1.fastq 2> SAMPLEA1.err.txt

mkdir -p ./raw
deinterleave_pairs SAMPLEA1.fastq \
   -o ./raw/SAMPLEA1_R1.fastq \
      ./raw/SAMPLEA1_R2.fastq

rm -f ./tmp.frags.r1.1.fasta ./frags_A1.fa ./SAMPLEA1.fastq ./SAMPLEA1.err.txt
```

#Para cada réplica o mesmo processo foi executado, e para facilitar criei um script por réplica: gerandofastqA1.sh, gerandofastqA2.sh, gerandofastqB1.sh, gerandofastqB2.sh. Assim como pode ser utilizado o scripts sim.sh, disponibilizado pelo professor.
#Onde: cat é utilizado para criar arquivos ou imprimir conteúdo na tela; mkdir -p cria diretórios (no caso, foi criado o diretório ./raw para receber os arquivos R.fastq); simNGS é um software para simmulação de sequenciamento next-generation (no caso, Illumina); -a para o adaptador (seguido da sequência de nucleotídeos); -p para considerar a corrida com paired-end; -n para indicar o número de ciclos para a corrida, considerando o tamanho das leituras de 151pb; deinterleave_pairs seguido do formato de entrada para desintercalar as sequências; deinterleave_pairs -o para indicar o arquivo de saída (output; no caso foram 2 arquivos: R1.fastq e R2.fastq).

#**4. A etapa seguinte foi baixar os arquivos com as sequêcnias e anotações referência, em formato .fna e .gff. Para isso é necessário acessar pelo NCBI o genoma da espécie em análise e consultar pelo Taxonomy ID (no caso, txid3702[orgn]). E em seguida os arquivos tipo .gz foram descompactados.**

```bash=
mkdir ./refs
cd ./refs
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.gff.gz

gunzip *.gz
```

#Ainda nessa etapa, foi baixado o arquivo .sra da espécie, utilizando o comando prefetch do NCBI Toolkit

```bash=
prefetch SRR10662152

fastq-dump \
   --split-3 \
   --origfmt SRR8707026
```

#Onde: --split-3 é para dados paired; --origfmt para manter o nome de origem
#A limpeza do arquivo fasta foi realizada pelo script cleanfasta.sh, dentro do diretório criado ./refs.

```bash=
#!/bin/bash
infile=$1
if [ ! ${infile} ]; then
        echo "Missing input fasta file"
        exit
fi
if [ ! -e ${infile} ]; then
        echo "Not found input fasta file (${infile})"
        exit
fi
sed 's/^>\(\S\+\).*/>\1/' ${infile}
./cleanfasta.sh GCF_000001735.4_TAIR10.1_genomic.fna > genome.fa
```
#Onde: $1 indica a linha do último arquivo de input ou linha; -e é usado para adicionar um script a ser executado, no caso o sed, que filtra conteúdos de interesse em um arquivo; s faz a combinação entre o especificado e o conteúdo a ser pesquisado; ^ delimita que tudo que aparecer após dele deve aparecer no arquivo/conteúdo; /.../ delimitam o que vai ser selecionado com as características especificadas entre as barras; .* delimita tudo que tiver as características dentro de () será de interesse; \+ indica que haverá um ou mais conteúdos delilitados encontrados

#Para o ajuste do arquivo GFF foi utilizado o script fixNCBIgff.sh para selecionar e delimitar o conteúdo desse arquivo em conteúdo descritivo dos genes.  

```bash=
fixNCBIgff.sh GCF_000001735.4_TAIR10.1_genomic.gff genome.gff
```

#Esse arquivo pode ser convertido em outros formatos, como GTF, FASTA do transcriptoma, FASTA do proteoma

```bash=
gffread genome.gff -g genome.fa -T -o genome.gtf
gffread genome.gff -g genome.fa -w transcriptome.fa
gffread genome.gff -g genome.fa -y proteome.fa
```
#Onde: -T para formato gtf do output; -w para escrever o arquivo FASTA com exons porcessados (spliced); -y para escrever a tradução do CDS.

##**Execução do script rnaseq-ref.sh**

#O pepiline utilizado nessa etapa foi modificado em algumas etapas, que serão descritas na sequência da explicação dos comandos. Dessa forma, o nome utilizado neste tópico será rnaseq-ref.modif.sh.

#**1. A primeira etapa do script consiste em conferir a existência dos diretórios e arquivos de entrada (input) e se apresentam conteúdo. São eles:**

#1. ./raw - com arquivos SAMPLEA1_R1.fastq, SAMPLEA2_R2.fastq, SAMPLEB2_R1.fastq, SRR10441850_2.fastq, SAMPLEA1_R2.fastq, SAMPLEB1_R1.fastq, SAMPLEB2_R2.fastq, SAMPLEA2_R1.fastq, SAMPLEB1_R2.fastq. 
#2. ./output - diretório para receber todas as análises.
#3. ./refs/genome.gtf - arquivo de extensão GTF dentro do diretório refs, com as sequêcnias e anotações referência.
#4. ./refs/genome.fa - arquivo de extensão FASTA dentro do diretório refs, com as sequêcnias e anotações referência.

#Essa etapa pode ser conferida previamente utilizando as seguintes linhas de comando:

```bash=
ls  -dlh ./raw ./output ./refs/genome.gtf ./refs/genome.fa
ls -dlh ./raw/*
wc -l raw/*.fastq
find ./raw -name '*_R[12].fastq' -exec bash -c 'echo -e "$0\t$(echo "$(cat $0 | wc -l)/4" | bc)"' {}  \;
find ./raw -name '*_R[12].fastq' -exec bash -c 'echo -e "$0\t$(cat $0 | paste - - - - | cut -f 1 | wc -l)"' {}  \;
```
#Os comandos no pepiline do rnaseq-ref.modif.sh para essa etapa 1 foram:

```bash=
num_threads=4
indir=$1
# SE ${indir} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO 1 NA LINHA DE COMANDO

if [ ! ${indir} ]; then
	echo "Missing input directory."
	exit
fi
# SE ${indir} NÃO É DIRETÓRIO
if [ ! -d ${indir} ]; then
	echo "Wrong input directory (${indir})."
	exit
fi

outdir=$2
# SE ${outdir} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO 2 NA LINHA DE COMANDO
if [ ! ${outdir} ]; then
	echo "Missing output directory."
	exit
fi

# SE ${outdir} NÃO É DIRETÓRIO
if [ ! -d ${outdir} ]; then
	echo "Wrong output directory (${outdir})."
	exit
fi
refgtf=$3

# SE ${refgtf} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO 3 NA LINHA DE COMANDO
if [ ! ${refgtf} ]; then
	echo "Missing GTF file."
	exit
fi

if [ ! -e "${refgtf}" ]; then
	echo "Not found GTF file (${refgtf})."
	exit
fi
refseq=$4

# SE ${refseq} NÃO EXISTE, SE NÃO FOI PASSADO ARGUMENTO 4 NA LINHA DE COMANDO
if [ ! ${refseq} ]; then
	echo "Missing GENOME fasta file."
	exit
fi

if [ ! -e "${refseq}" ]; then
	echo "Not found GENOME fasta file (${refseq})."
	exit
fi
```

#**2. A etapa seguinte é a execução do pepiline preprocess3.sh:**

```bash=
./preprocess3.sh "${indir}" "${outdir}"

mkdir -p ${outdir}/star_index
mkdir -p ${outdir}/star_out_pe
mkdir -p ${outdir}/star_out_se
mkdir -p ${outdir}/star_out_final
mkdir -p ${outdir}/cufflinks
mkdir -p ${outdir}/cuffmerge
mkdir -p ${outdir}/stringtie
mkdir -p ${outdir}/stringmerge
mkdir -p ${outdir}/cuffcompare
mkdir -p ${outdir}/cuffquant

for r1 in `ls ${outdir}/processed/prinseq/*.atropos_final.prinseq_1.fastq`; do
	r1_singletons=`echo ${r1} | sed 's/prinseq_1.fastq/prinseq_1_singletons.fastq/'`

	if [ ! -e "${r1_singletons}" ]; then
		touch ${r1_singletons}
	fi

	r2=`echo ${r1} | sed 's/prinseq_1.fastq/prinseq_2.fastq/'`

	if [ ! -e "${r2}" ]; then
		echo "Read 2 (${r2}) paired with Read 1 ($r1) not found."
		exit
	fi

	r2_singletons=`echo ${r2} | sed 's/prinseq_2.fastq/prinseq_2_singletons.fastq/'`

	if [ ! -e "${r2_singletons}" ]; then
		touch ${r2_singletons}
	fi

	name=`basename ${r1} | sed 's/.atropos_final.prinseq_1.fastq//'`
```
#Inicialmente foram criados os diretórios que receberão todos os resultados das análises. O comando utilizado foi *mkdir -p*, já descritos anteriormente.
#Neste pepiline, foi realizado o pré-processamento utilizando FastQC (pré), ATROPOS, PrinSeq, FastQC (pós).

#Inicialmente foram criados os diretórios para receber todos os resultados das análises: ${outdir}/processed/fastqc/pre, /processed/prinseq, /processed/fastqc/pos. O comando utilizado foi *mkdir -p*, já descritos anteriormente. 

#O FastQC (pré) verifica a qualidade dos dados do sequenciamento. Para isso, analisa os dados dos arquivos de formato BAM, SAM ou FastQ; no caso da atividade atual, são arquivos de input os de nomenclatura _R1.fastq e _R2.fastq, do diretório ./raw.

#Alguns comandos importantes: -t 2 para especificar que seriam processados simultaneamente dois arquivos (files); -o para especificar o local e o nome do arquivo de saída, no caso, em ${outdir}/processed/fastqc/pre/, sendo ${outdir} para o diretório ./output. 
:::

#O ATROPOS detecta e realiza a trimagem dos adaptadores ainda presentes nas sequências provenientes do sequenciamento, no caso oa arquivos R1 e R2, podendo utilizar formatos BAM, SAM ou FastQ.

#Comandos: --aligner insert para identificar adaptadores de extremidade 3’ de sequências paired-end, em alinhamentos baseados em inserção, de forma mais acurada; -e 0.1 corresponde à taxa de erro global, usado para especificar a acurácia e ajustar a correspondência de inserção de adaptador; -n define o número de sequências (cópias) do adaptador a serem removidas de cada read; -m para especificar o tamanho mínimo de reads trimados a serem descartados; --op-order controla a ordem de aplicação da primeira operação de trimagem; --match-read-wildcards define a correspondência de sequências com bases N na read; -O especifica o tamanho mínimo de uma sobreposição entre o adaptador e a read; -q para trimar extremidades de baixa qualidade de reads sem a remoção de adaptadores; -T define o número de núcleos para processar a trimagem; --correct-mismatches conservative estabelece a correção de mismatches durante o alinhamento de acordo com a qualidade da base; define a base como inalterada quando a qualidade é igual; --pair-filter any define que qualquer sequência de uma read paired-end deve receber o critério de filtragem para que seja filtrado; -a para especificar o adaptador ligado à extremidade 3’ para ser removido da sequência de leitura 1 do par de sequências, se utilizado o ‘$’ em seguida, o adaptador é apenas encontrado;-A para especificar o adaptador ligado à extremidade 3’ para ser removido da sequência de leitura 2 do par de sequências; -o especifica arquivos de saída de sequências trimadas, em formatos FASTQ ou FASTA; -p especifica arquivos de saída para sequências paired-end de read 2; -pe1 seguido do endereço do arquivo .fastq especifica o arquivo de read 1; se utilizado ‘-1’, arquivos no formato SAM, BAM também são suportados como arquivos de entrada; -pe2 seguido do endereço do arquivo .fastq, especifica o arquivo de read 2; --untrimmed-output obtém a segunda sequência de um read tipo paired-end trimado; --untrimmed-paired-output obtém a segunda sequência de um read, quando o adaptador já foi encontrado nesse primeiro read; cat para concatenar os resultados; rm -f para excluir arquivos desnecessários já concatenados em .atropos_final.fastq, no caso os de extensão R1/R2.atropos_insert.fastq, _adapter.fastq.  

#Os resultados do ATROPOS foram concatenados e alguns arquivos desnecessário foram excluídos. O local de saída foi ./output/processed/atropos/, gerando 8 arquivos para cada SAMPLE (A1-B2) com extensão _R1 e _R2 mais.

#Em seguida foi executado o Prinseq (prinseq-lite.pl), que trima, filtra, reformata e gera qualidade de sequências derivadas de sequenciamento next-generation.

#Comandos: -fastq e -fastq2 indicam os endereços dos arquivos de entrada, no caso R1 e R2.atropos_final.fastq; -out_format 3 para mudar o formato do arquivo de saída, no caso o 3 indica a opção para arquivo tipo FASTQ; -trim_qual_window 3 indica o tamanho da janela usada para calcular os valores (scores) de qualidade; -trim_qual_step 1 determina o tamanho da etapa usada para mover a janela de deslizamento, no caso 1 indica mover a janela sobre todas as contagens da qualidade sem faltar nenhum, com o tamanho menor/igual a janela (window); -trim_qual_right 30 trima a sequência pelo score de qualidade a partir da extremidade 3', no caso o valor de threshold é 30, ou seja, o valor mínimo ideal de qualidade; -trim_qual_type  mean indica que o tipo de cálculo de pontuação a ser utilizado é a média; 	-trim_qual_rule lt indica que a regra a ser usada para comparar o score	é menor que (less than); -out_good seguido de um endereço de diretório especifica onde o arquivo de saída será criado; -out_bad null para previnir que o programa não gere os arquivos output sem passar filtros no conteúdo de dados; -lc_method dust para calcular a complexidade das sequências usando o método dust; -lc_threshold 30 especifica o valor de threshold para filtrar as sequências pela complexidade, no caso o valor de threshold é 30; -min_len 20 filtra sequências com tamanho mínimo de 20pb; -trim_tail_right 5 usado para trimar a calda poliA/T com o tamanho mínimo de 5 na extremidade 3';  -trim_tail_left 5 usado para trimar a calda poliA/T com o tamanho mínimo de 5 na extremidade 5'; -ns_max_p 80 filtra sequências com mais de 80% de Ns; -noniupac filtra sequências outros caracteres além A, C, G, T ou N.

#Por fim, no script de preprocess3.sh, foi executado o FastQC (pós), que verifica a qualidade dos dados a partir dos arquivos  do prinseq, com extenção .atropos_final.prinseq_1.fastq. 

#**3. Na sequência do rnaseq-ref.modif.sh, ocorre a indexação do genoma para sequente execução do alinhador STAR.
Nessa etapa, as sequências geradas pelo sequenciamento são alinhadas com o genoma referência.** 

```bash=
if [ ! -e "${outdir}/star_index/SAindex" ]; then
		echo "Indexing genome (${refseq}) ..."		# --genomeSAindexNbases 12 (sugestão do alinhador)
		# --sjdbOverhang 149 (sugestão do manual)	

	STAR    --runThreadN        ${num_threads} \
			--runMode           genomeGenerate \
			--genomeFastaFiles  ${refseq} \
			--genomeDir         ${outdir}/star_index \
			--sjdbGTFfile       ${refgtf} \
			--genomeSAindexNbases 12 \
			--sjdbOverhang      149 \
		 > ${outdir}/star_index/STAR.index.log.out.txt \
		2> ${outdir}/star_index/STAR.index.log.err.txt

	fi
	echo "STAR alignment PE with sample ${name}: ${r1} & ${r2} ..."

	# --outSAMstrandField intronMotif 
	# --outFilterIntronMotifs RemoveNoncanonical 
	# (parâmetros recomendados pelo Manual para manter a compatibilidade com Cufflinks)

	mkdir -p ${outdir}/star_out_pe/${name}
	STAR	--runThreadN        ${num_threads} \
		    --genomeDir         ${outdir}/star_index \
		    --readFilesIn       ${r1} ${r2} \
		    --outSAMstrandField intronMotif \
		    --outFilterIntronMotifs RemoveNoncanonical \
            --sjdbGTFfile       ${refgtf} \
		    --outFilterMultimapNmax 20 \
		    --outFileNamePrefix ${outdir}/star_out_pe/${name}/ \
		    --outSAMtype        BAM Unsorted \
		 > ${outdir}/star_out_pe/${name}/STAR.alignment_pe.log.out.txt \
		2> ${outdir}/star_out_pe/${name}/STAR.alignment_pe.log.err.txt

	echo "STAR alignment SE with sample ${name}: ${r1_singletons} & ${r2_singletons} ..."

	mkdir -p ${outdir}/star_out_se/${name}

	STAR --runThreadN        ${num_threads} \
		 --genomeDir         ${outdir}/star_index \
		 --readFilesIn       ${r1_singletons},${r2_singletons} \
		 --sjdbGTFfile       ${refgtf} \
		 --outSAMtype        BAM Unsorted \
		 --outFilterMultimapNmax 20 \
		 --outSAMstrandField intronMotif \
		 --outFileNamePrefix ./$outdir/star_out_se/${name}/ \
		 > ./${outdir}/star_out_se/${name}/STAR.alignment_se.log.out.txt \
		2> ./${outdir}/star_out_se/${name}/STAR.alignment_se.log.err.txt

	echo "Merging STAR alignment PE & SE (${name}) ..."
	mkdir -p ${outdir}/star_out_final/${name}
        # Combinar resultados do alinhamento com reads paired-end e alinhamento com reads single-end (singletons)   

	 samtools merge -@ ${num_threads} -f -n  ${outdir}/star_out_final/${name}/Aligned.out.bam \              ${outdir}/star_out_pe/${name}/Aligned.out.bam \
 ${outdir}/star_out_se/${name}/Aligned.out.bam \
	 >  ${outdir}/star_out_final/${name}/samtools.merge.log.out.txt \
	2>  ${outdir}/star_out_final/${name}/samtools.merge.log.err.txt

	echo "Sorting STAR alignment final (${name}) ..."
        # Ordenando o resultado do alinhamento por coordenadas genômicas
        # - exigência para executar o cufflinks

	 samtools sort -@ ${num_threads} -o      ${outdir}/star_out_final/${name}/Aligned.out.sorted.bam \        ${outdir}/star_out_final/${name}/Aligned.out.bam \
	 >  ${outdir}/star_out_final/${name}/samtools.sort.log.out.txt \
	2>  ${outdir}/star_out_final/${name}/samtools.sort.log.err.txt

	echo "Collecting alignment statistics (${name}) ..."

	SAM_nameSorted_to_uniq_count_stats.pl  ${outdir}/star_out_final/${name}/Aligned.out.bam >  ${outdir}/star_out_final/${name}/Aligned.stats.txt
```
#Onde: --runThreadN para o número de seguimentos, no caso 2;	--genomeDir especifica o diretório onde os arquivos analisados serão criados; --readFilesIn especifica os locais dos arquivos de entrada (input), no caso _R1 e _R2; --outSAMstrandField intronMotif para filtrar introns inconsistentes, sendo uma flag para o Cufflinks, que virá a seguir da execução do STAR; --outFilterIntronMotifs RemoveNoncanonical para filtrar todos os alinhamentos usando os seus próprios motifs, recomendação para manter compatibilidade com o Cufflinks; --sjdbGTFfile especifica o caminho para o arquivo GTF, no caso, dentro do diretório ./refs; outFilterMultimapNmax 20 especifica o número máximo de loci para mapear que serão exigidos para serem gerados os output, no caso 20; --outFileNamePrefix especifica o nome do arquivo output, no caso ${outdir}/star_out_pe/${name}; --outSAMtype BAM Unsorted especifica o tipo de formato do arquivo output a ser gerado, no caso BAM sem coordenadas padrão diferenciadas; samtools merge usado para combinar os resultados os reads paired-end e o alinhamento com reads single-end (singletons); samtools sort para ordenar os resultados por coordenadas genômicas; - é uma exigência para o Cufflinks.

#**4. A execução do pepiline segue com o Cufflinks.**

#Cufflinks é utilizado principalmente para montar os transcritos por alinhamento dos reads no transcriptoma referência, para cada amostra (SAMPLE) nas suas respectivas condições biológicas. Além disso testar a expressão diferencial dos dos transcritos por suas abundâncias.

#**Nesta etapa, foi necessária a modificação de alguns parâmetros, baseando-se nas informações de Wang (2005), sobre o tamanho médio dos introns em *Arabidopsis thaliana*. O tamanho médio dos fragmentos (--frag-len-mean) foi alterado de 300 para 173, e o tamanho mínimo do intron foi modificado de 100 para 50.**

```bash=
echo "Running Cufflinks (${name}) ..."
	mkdir -p ${outdir}/cufflinks/${name}

	cufflinks --output-dir ${outdir}/cufflinks/${name} \
		  --num-threads ${num_threads} \
		  --GTF-guide ${refgtf} \
		  --frag-bias-correct ${refseq} \
		  --multi-read-correct \
		  --library-type fr-unstranded \
		  --frag-len-mean 173 \
		  --frag-len-std-dev 50 \
		  --total-hits-norm \
		  --max-frag-multihits 20 \
		  --min-isoform-fraction 0.20 \
		  --max-intron-length 10000 \
		  --min-intron-length 50 \
		  --overhang-tolerance 4 \
		  --max-bundle-frags 999999 \
		  --max-multiread-fraction 0.45 \
		  --overlap-radius 10 \
		  --3-overhang-tolerance 300 \
		  ${outdir}/star_out_final/${name}/Aligned.out.sorted.bam \
		 > ${outdir}/star_out_final/${name}/cufflinks.log.out.txt \
		2> ${outdir}/star_out_final/${name}/cufflinks.log.err.txt
```
#Comandos: --output-dir especifica o local para gerar os arquivos de saída; --num-threads especifica o número de segmentos; --GTF-guide para usar o arquivo de anotação dos transcritos referência para guiar a montagem; --frag-bias-correct para correção do viés de uso, com o arquivo FASTA como referência; --multi-read-correct especifica oo uso do "método de resgate" para os múltiplos reads, dando maior acurácia; --library-type fr-unstranded para definir o tipo de livraria de reads para input, no caso fragmentos não ociosos; --frag-len-mean 173 para definir o tamanho médio dos fragmentos não pareados, no caso 173pb; --frag-len-std-dev 50 para definir o desvio do tamanho dos fragmentos std não pareados, no caso 50pb; --total-hits-norm para contar todos os hits para a normalização; --max-frag-multihits 20 para definir o número máximo de alinhamentos  por fragmento, no caso 20; --min-isoform-fraction 0.20 define o nível mínimo de abundância  dos transcritos para suprimir tas ranscrições abaixo deste nível, no caso 0.20; --max-intron-length 10000 define o tamanho dos gaps a serem ignorados, no caso ignora-se introns maiores que 10000pb; --min-intron-length 50 define o tamanho mínimo dos introns permitidos no genoma, no caso 50pb; --overhang-tolerance 4 define o número do exon terminal para tolerar em introns, no caso 4; --max-bundle-frags 999999 fragmentos máximos permitidos em um pacote antes de ignorar; --max-multiread-fraction 0.45 define a fração máxima de leituras múltiplas permitidas por transcrição, no caso 0.45; --overlap-radius 10 para o tamanho máximo dos gaps para preencher entre os fragmentos, no caso 10pb; --3-overhang-tolerance 300 define o tamanho do oveerhang ao final da extremidade 3' quando fundido à referência.

#**5. Em seguida é executado o StringTie.**
#StringTie é um montador que gera arquivos que podem ser processados pelo Cuffdiff ou outros programas.

```bash=
echo "Running StringTie (${name}) ..."
	mkdir -p ${outdir}/stringtie/${name}

	stringtie ${outdir}/star_out_final/${name}/Aligned.out.sorted.bam \
		-G ${refgtf} \
		-o ${outdir}/stringtie/${name}/transcripts.gtf \
		-p ${num_threads} \
		-f 0.20 \
		-a 10 \
		-j 3 \
		-c 2 \
		-g 10 \
		-M 0.45 \
    	-A ${outdir}/stringtie/${name}/gene_abundance.txt
done
```
#Onde: -G especifica a referência de anotação a ser usada para guiar a montagem, no caso o arquivo GTF; -o especifica o local e o nome do arquivo de saída, no caso ./output/stringtie/${name}/transcripts.gtf; -p para o número de sgmentos a serem usados; -f 0.20 define que a fração mínima da isoforma é 0.20; -a 10 define o comprimento mínimo de ancoragem para junções, no caso 10; -j 3 define o mínimo de cobertura das junções, no caso 3; -c 2 mínimo de cobertura de reads por pb para considerar a montagem; -g 10 para o gap entre os reads mapeados marcando um novo novo agrupamento, no caso 10; -M 0.45 para a fração de agrupamento autorizado a ser coberto por reads multi-hits; -A para o arquivo output de abundância e o local de saída.

#**6. O script segue para a execução do Cuffmerge e do Stringtie merge.**
#Cuffmerge faz a fusão entre os arquivos GTFs gerados pelo Cufflinks, e cria um arquivo correspondente a um catálogo de transcritos. É possível também, utilizar um arquivo GTF referência, onde o Cuffmerge irá anexar as informações dos genes no catálogo.

#O Stringtie merge permite fundir transcritos montados de múltiplos arquivos de entrada (input) em um conjunto de isoformas não redundantes.

```bash=
echo "Running cuffmerge ..."

find ${outdir}/cufflinks/ -name 'transcripts.gtf' > ${outdir}/cuffmerge/assembly_GTF_list.txt

cuffmerge -o ${outdir}/cuffmerge/ \
	--ref-gtf ${refgtf} \
	--ref-sequence ${refseq} \
	--min-isoform-fraction 0.20 \
	--num-threads ${num_threads} \
	   ${outdir}/cuffmerge/assembly_GTF_list.txt \
	 > ${outdir}/cuffmerge/cuffmerge.log.out.txt \
	2> ${outdir}/cuffmerge/cuffmerge.log.err.txt

echo "Running stringtie merge ..."
find ${outdir}/stringtie/ -name 'transcripts.gtf' > ${outdir}/stringmerge/assembly_GTF_list.txt

stringtie --merge \
	-G ${refgtf} \
	-o ${outdir}/stringmerge/merged.gtf \
	-c 1 \
	-T 1 \
	-f 0.20 \
	-g 10 \
	-i \
	${outdir}/stringmerge/assembly_GTF_list.txt
```
#Comandos: no Cuffmerge -o para especificar o diretório de saída; --ref-gtf indica uso de arquivo GTF referência; --ref-sequence especifica o arquivo com as sequências de DNA para referência; --min-isoform-fraction 0.20 define a abundância limite para discarte das isorformas, no caso isoformas com menos de 0.20 de abundância foram discartadas; --num-threads usado para fazer a fusão das montagens.
          no Stringtie merge -G especifica o arquivo de anotações referência para incluir na fusão, no caso GTF; -o especifica o local e o nome do arquivo de saída, no caso ./output/stringmerge/merged.gtf; -c 1 indica a cobertura mínima de transcritos para incluir na fusão, no caso 1; -T 1 indica o mínimo de TPM para incluir na fusão; -f 0.20 indica a fração mínima de isoformas, no caso 0.20; -g 10 indica o gap entre transcritos a ser considerado para fundir juntos, no caso 10; -i permite manter transcrições fundidas com introns retidos.

#**7. As execuções final envolvem o Cuffcompare, Cuffquant, Cuffnorm e cuffdiff**   

#O Cuffcompare utiliza os arquivos do alinhamento com STAR (GTFs) para quantificar todas as réplicas.  
#O Cuffquant quantifica a expressão gênica ou de transcritos do RNA-seq de cada réplica, para ser analizada posteriormente pelo Cuffdiff ou Cuffnorm. Para suas análises utiliza os arquivos gerados pelo STAR.
#O Cuffnorm normaliza os níveis de expressão dos genes/transcritos de RNA-Seq, deixando o conjunto de dados na mesma escala. Os arquivos gerados, matrizes de contagem, são em formato de texto.
#O Cuffdiff executado na sequência realiza a análise de cálculo da expressão diferencial.

```bash=
cuffcompare	-r ${refgtf} \
		    -s ${refseq} \
		    -o ${outdir}/cuffcompare/stringmerge \
		   ${outdir}/stringmerge/merged.gtf \
		 > ${outdir}/stringmerge/cuffcompare.log.out.txt \
		2> ${outdir}/stringmerge/cuffcompare.log.err.txt

biogroup_label=()
for bamfile in `ls ${outdir}/star_out_final/*/Aligned.out.sorted.bam`; do
	name=`basename $(dirname ${bamfile})`

echo "Running cuffquant using sample ${name} with ${outdir}/stringmerge/merged.gtf as reference ..."
	mkdir -p ${outdir}/cuffquant/${name}

	cuffquant 	--output-dir ${outdir}/cuffquant/${name} \
			    --frag-bias-correct ${refseq} \
			    --multi-read-correct \
			    --num-threads ${num_threads} \
			    --library-type fr-unstranded \
			    --frag-len-mean 300 \
			    --frag-len-std-dev 50 \
			    --max-bundle-frags 9999999 \
			    --max-frag-multihits 20 \
		   ${outdir}/stringmerge/merged.gtf \
		   ${bamfile} \
		 > ${outdir}/cuffquant/${name}/cuffquant.log.out.txt \
		2> ${outdir}/cuffquant/${name}/cuffquant.log.err.txt

	groupname=`echo ${name} | sed 's/[0-9]\+$//'`
	biogroup_label=($(printf "%s\n" ${biogroup_label[@]} ${groupname} | sort -u ))

done

biogroup_files=()

echo "Running Differential Expression Analysis ..."
for label in ${biogroup_label[@]}; do
	echo -e "\tCollecting .cxb files for ${label} ..."
	group=()

	for cxbfile in `ls ${outdir}/cuffquant/${label}*/abundances.cxb`; do
		echo -e "\t\tFound ${cxbfile}"
		group=(${group[@]} "${cxbfile}")
biogroup_files=()

echo "Running Differential Expression Analysis ..."
for label in ${biogroup_label[@]}; do
	echo -e "\tCollecting .cxb files for ${label} ..."
	group=()

	for cxbfile in `ls ${outdir}/cuffquant/${label}*/abundances.cxb`; do
		echo -e "\t\tFound ${cxbfile}"
		group=(${group[@]} "${cxbfile}")
	done

	biogroup_files=(${biogroup_files[@]} $(IFS=, ; echo "${group[*]}") )
done

echo -e "\tRunning cuffnorm & cuffdiff ..."
echo -e "\t\tLabels.: " $(IFS=, ; echo "${biogroup_label[*]}")
echo -e "\t\tFiles..: " ${biogroup_files[*]}
echo -e "\t\t\tGenerating abundance matrices (cuffnorm) ..."

mkdir -p ${outdir}/cuffnorm/

cuffnorm 	--output-dir ${outdir}/cuffnorm \
 		    --labels $(IFS=, ; echo "${biogroup_label[*]}") \
 		    --num-threads ${num_threads} \
		    --library-type fr-unstranded \
 		    --library-norm-method geometric \
		    --output-format simple-table \
 		   ${outdir}/stringmerge/merged.gtf \
 		   ${biogroup_files[*]} \
 	 	 > ${outdir}/cuffnorm/cuffdiff.log.out.txt \
 		2> ${outdir}/cuffnorm/cuffdiff.log.err.txt

echo -e "\t\t\tAnalysing differential expression (cuffdiff) ..."

mkdir -p ${outdir}/cuffdiff/

cuffdiff 	--output-dir ${outdir}/cuffdiff \
 		    --labels $(IFS=, ; echo "${biogroup_label[*]}") \
 		    --frag-bias-correct ${refseq} \
 		    --multi-read-correct \
 		    --num-threads ${num_threads} \
 		    --library-type fr-unstranded \
 		    --frag-len-mean 173 \
 		    --frag-len-std-dev 50 \
 		    --max-bundle-frags 9999999 \
 		    --max-frag-multihits 20 \
 		    --total-hits-norm \
	done

	biogroup_files=(${biogroup_files[@]} $(IFS=, ; echo "${group[*]}") )
done

echo -e "\tRunning cuffnorm & cuffdiff ..."
echo -e "\t\tLabels.: " $(IFS=, ; echo "${biogroup_label[*]}")
echo -e "\t\tFiles..: " ${biogroup_files[*]}
echo -e "\t\t\tGenerating abundance matrices (cuffnorm) ..."

mkdir -p ${outdir}/cuffnorm/

cuffnorm 	--output-dir ${outdir}/cuffnorm \
 		    --labels $(IFS=, ; echo "${biogroup_label[*]}") \
 		    --num-threads ${num_threads} \
		    --library-type fr-unstranded \
 		    --library-norm-method geometric \
		    --output-format simple-table \
 		   ${outdir}/stringmerge/merged.gtf \
 		   ${biogroup_files[*]} \
 	 	 > ${outdir}/cuffnorm/cuffdiff.log.out.txt \
 		2> ${outdir}/cuffnorm/cuffdiff.log.err.txt

echo -e "\t\t\tAnalysing differential expression (cuffdiff) ..."

mkdir -p ${outdir}/cuffdiff/

cuffdiff 	--output-dir ${outdir}/cuffdiff \
 		    --labels $(IFS=, ; echo "${biogroup_label[*]}") \
 		    --frag-bias-correct ${refseq} \
 		    --multi-read-correct \
 		    --num-threads ${num_threads} \
 		    --library-type fr-unstranded \
 		    --frag-len-mean 173 \
 		    --frag-len-std-dev 50 \
 		    --max-bundle-frags 9999999 \
 		    --max-frag-multihits 20 \
 		    --total-hits-norm \
 		    --min-reps-for-js-test 2 \
 		    --library-norm-method geometric \
 		    --dispersion-method per-condition \
 		    --min-alignment-count 10 \

 		  ${outdir}/stringmerge/merged.gtf \
 		  ${biogroup_files[*]} \
 	 	> ${outdir}/cuffdiff/cuffdiff.log.out.txt \
 	   2> ${outdir}/cuffdiff/cuffdiff.log.err.txt
```

#Onde: no Cuffcompare, -r especifica o arquivo de input GTF; -s especifica o arquivo com as sequências referência, podendo ser multi-fasta; -o especifica o diretório de saída, no caso .output/cuffcompare/stringmerge. 
      no Cuffquant, --output-dir especifica o diretório de saída dos arquivos gerados; --frag-bias-correct inidica correção do viés de uso com arquivo FASTA referência; --multi-read-correct utiliza o método de resgate para reads múltiplos; --num-threads especifica o número de segmentos a serem usados na quantificação; --library-type fr-unstranded para definir o tipo de livraria de reads para input, no caso fragmentos não ociosos; --frag-len-mean 300 indica o tamanho médio do fragmento, considerando apenas reads não pareados, no caso 300pb; --frag-len-std-dev 50 para definir o desvio do tamanho dos fragmentos std não pareados, no caso 50pb; --max-bundle-frags 9999999 ragmentos máximos permitidos em um pacote antes de ignorar; --max-frag-multihits 20 para o número máximo de alinhamentos por fragmento, no caso 20.
      no Cuffnorm, --output-dir especifica o local para gerar os arquivos de saída; --labels para lista  de etiquetas; --num-threads especifica o número de segmentos; --library-type fr-unstranded para definir o tipo de livraria de reads para input, no caso fragmentos não ociosos; --library-norm-method geometric especifica o método geométrico utilizado para normalizar as livrarias de tamanhos; --output-format simple-table especifica a tabela simples como tipo de arquivo de saída.

#A visualização dos resultados do script rnaseq-ref.modif.sh pode ser feita pelos comandos:
#a. Contagem das sequências durante as etapa de processamento
```bash=
#para os dados do diretório raw
find ./raw -name '*_R[12].fastq' -exec bash -c 'echo -e "$0\t$(echo "$(cat $0 | wc -l)/4" | bc)"' {}  \;

#para os dados do atropos
find ./atropos -name '*_R[12].atropos_final.fastq' -exec bash -c 'echo -e "$0\t$(echo "$(cat $0 | wc -l)/4" | bc)"' {}  \;

#para os dados do prinseq
 find ./ -name '*.atropos_final.prinseq_[12]*.fastq' -exec bash -c 'echo -e "$0\t$(echo "$(cat $0 | wc -l)/4" | bc)"' {}  \;
```

#b. Alinhamento das sequências em relação as referências
```bash=
#Contagem de reads alinhanhas na referência
 samtools view Aligned.out.bam | cut -f1,3| sed '/*/d' | wc -l

#genes diferencialmente expressos:
cd cuffdiff
grep "yes" gene_exp.diff
```
