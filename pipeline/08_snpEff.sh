#!/usr/bin/bash
#SBATCH --mem=128G -p intel,batch --nodes 1 --ntasks 8 --out logs/snpEff.log
module unload miniconda2
module load miniconda3
module load snpEff
module load tabix
module load yq

# THIS IS AN EXAMPLE OF HOW TO MAKE SNPEFF - it is for A.fumigatus
GFFGENOME=$SNPEFFGENOME.gff
CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
  CPU=1
fi
MEM=64g
CHROMPREF=JOWA01
# this module defines SNPEFFJAR and SNPEFFDIR
if [ -f config.txt ]; then
  source config.txt
fi
mapDomains=$(realpath scripts/map_snpEff2domains.py)
snpEffTab=$(realpath scripts/snpEff_2_tab.py)

DOMAINS=genome/$(basename ${GFFGENOME} .gff)_InterproDomains.txt
DOMAINS=$(realpath $DOMAINS)
GFFGENOMEFILE=$GENOMEFOLDER/$GFFGENOME
FASTAGENOMEFILE=$GENOMEFOLDER/$GENOMEFASTA
REFGENOME=$(realpath $REFGENOME)
echo "refgenome is $REFGENOME"
if [ -z $SNPEFFJAR ]; then
  echo "need to defined \$SNPEFFJAR in module or config.txt"
  exit
fi

if [ -z $SNPEFFDIR ]; then
  echo "need to defined \$SNPEFFDIR in module or config.txt"
  exit
fi
SNPEFFDIR=$(realpath $SNPEFFDIR)
# could make this a config

if [ -z $FINALVCF ]; then
  echo "need a FINALVCF variable in config.txt"
  exit
fi
FINALVCF=$(realpath $FINALVCF)
if [[ -z $POPYAML || ! -s $POPYAML ]]; then
  echo "Need populations.yaml file and variable POPYAML set in the config.txt"
  exit
fi
SNPEFFOUT=$(realpath $SNPEFFOUT)
mkdir -p $SNPEFFOUT
## NOTE YOU WILL NEED TO FIX THIS FOR YOUR CUSTOM GENOME
if [ ! -e $SNPEFFOUT/$snpEffConfig ]; then
  rsync -a $SNPEFFDIR/snpEff.config $SNPEFFOUT/$snpEffConfig
  echo "# Sced.fungidb " >> $SNPEFFOUT/$snpEffConfig
  # CHANGE Aspergillus fumigatus Af293 FungiDB to your genome name and source - though this is really not important - $SNPEFFGENOME.genome is really what is used
  echo "$SNPEFFGENOME.genome : SapiospermumIHEM14462 FungiDB" >> $SNPEFFOUT/$snpEffConfig
  chroms=$(grep ">" $REFGENOME | awk '{print $1}' | perl -p -e "s/>${CHROMPREF}0+/Chr/" | sort | perl -p -e 's/\n/, /' |  perl -p -e 's/,\s+$/\n/')
  echo -e "\t$SNPEFFGENOME.chromosomes: $chroms" >> $SNPEFFOUT/$snpEffConfig
  # THIS WOULD NEED SPEIFIC FIX BY USER - IN A.fumigatus the MT contig is called mito_A_fumigatus_Af293
  #echo -e "\t$SNPEFFGENOME.mito_SapiospermumIHEM14462.codonTable : Mold_Mitochondrial" >> $SNPEFFOUT/$snpEffConfig
  mkdir -p $SNPEFFOUT/data/$SNPEFFGENOME
  perl -p -e "s/${CHROMPREF}0+/Chr/" $GFFGENOMEFILE | gzip -c > $SNPEFFOUT/data/$SNPEFFGENOME/genes.gff.gz
  #rsync -aL $REFGENOME $SNPEFFOUT/data/$SNPEFFGENOME/sequences.fa
  perl -p -e "s/>${CHROMPREF}0+/>Chr/" $REFGENOME > $SNPEFFOUT/data/$SNPEFFGENOME/sequences.fa

  java -Xmx$MEM -jar $SNPEFFJAR build -datadir $SNPEFFOUT/data -c $SNPEFFOUT/$snpEffConfig -gff3 -v $SNPEFFGENOME
fi
# get full path to YAML file
POPYAML=$(realpath $POPYAML)
echo "refgenome is $REFGENOME YAML is $POPYAML"
pushd $SNPEFFOUT

makeMatrix() {
  POPNAME=$1
  echo "POPNAME is $POPNAME"
  mkdir -p $POPNAME
  module load bcftools/1.12
  COMBVCF="$FINALVCF/$PREFIX.$POPNAME.SNP.combined_selected.vcf.gz $FINALVCF/$PREFIX.$POPNAME.INDEL.combined_selected.vcf.gz"
  echo "COMBVCF is '$COMBVCF'"
  for n in $COMBVCF
  do
    echo $n
    st=$(echo $n | perl -p -e 's/\.gz//')
    if [ ! -f $n ]; then
      bgzip $st
    fi
    if [ ! -f $n.tbi ]; then
      tabix $n
    fi
  done

  pushd $POPNAME
  INVCF=$PREFIX.$POPNAME.allvariants_combined_selected.vcf
  OUTVCF=$PREFIX.$POPNAME.snpEff.vcf
  SEGOUTVCF=$PREFIX.$POPNAME.snpEff.segregating.vcf
  OUTTAB=$PREFIX.$POPNAME.snpEff.tab
  OUTMATRIX=$PREFIX.$POPNAME.snpEff.matrix.tsv
  SEGOUTMATRIX=$PREFIX.$POPNAME.snpEff.segregating.matrix.tsv
  DOMAINVAR=$PREFIX.$POPNAME.snpEff.domain_variant.tsv
  SEGDOMAINVAR=$PREFIX.$POPNAME.snpEff.segregating.domain_variant.tsv
  if [[ ! -s $INVCF || $FINALVCF/$PREFIX.$POPNAME.SNP.combined_selected.vcf.gz -nt $INVCF ]]; then
    bcftools concat -a -d both -o $INVCF -O v $COMBVCF
    perl -i -p -e "s/${CHROMPREF}0+/Chr/" $INVCF
  fi
  if [[ ! -f $OUTVCF || $INVCF -nt $OUTVCF ]]; then 
  #echo " java -Xmx$MEM -jar $SNPEFFJAR eff -c $SNPEFFOUT/$snpEffConfig -dataDir $SNPEFFOUT/data -v $SNPEFFGENOME $INVCF > $OUTVCF"
  java -Xmx$MEM -jar $SNPEFFJAR eff -c $SNPEFFOUT/$snpEffConfig -dataDir $SNPEFFOUT/data -v $SNPEFFGENOME $INVCF > $OUTVCF
  fi

  perl -i -p -e "if (/Chr(\d+)/) { my \$n=sprintf('%s%06d',$CHROMPREF,\$1); s/Chr(\d+)/\$n/}" $OUTVCF
  bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT{0}[\t%TGT]\t%INFO/ANN\n' $OUTVCF > $OUTTAB
  bcftools view -e 'QUAL < 1000 || AF=1' -c 1 -o $SEGOUTVCF $OUTVCF

  # this requires python3 and vcf script
  # this assumes the interpro domains were downloaded from FungiDB and their format - you will need to generalize this
  $mapDomains --vcf $OUTVCF --domains $DOMAINS --output $DOMAINVAR
  $mapDomains --vcf $SEGOUTVCF --domains $DOMAINS --output $SEGDOMAINVAR

  # this requires Python and the vcf library to be installed
  $snpEffTab $OUTVCF $REFGENOME > $OUTMATRIX
  $snpEffTab $SEGOUTVCF $REFGENOME > $SEGOUTMATRIX
  popd
}
source $(which env_parallel.bash)
export -f makeMatrix
env_parallel -j 1 --env _ makeMatrix ::: $(yq eval '.Populations | keys' $POPYAML | perl -p -e 's/^\s*\-\s*//' )
