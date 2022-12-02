progdir=~/Lab/Tools/hapgen2
refdir=~/Lab/Reference/hapmap3_r2_b36
hapdir=~/Lab/GenotypeDataSimulation/Data/Hapgen/VerySmallBatch

# Can't get hapgen to work witout the -dl flag
# so create dummy variable, that can be passed on

for chr in $(seq 1 22); do
	dummyDL=`sed -n '2'p $refdir/hapmap3_r2_b36_chr${chr}.legend | cut -d ' ' -f 2`
	$progdir/hapgen2 -m $refdir/genetic_map_chr${chr}_combined_b36.txt \
		-l $refdir/hapmap3_r2_b36_chr${chr}.legend \
		-h $refdir/hapmap3_r2_b36_chr${chr}.haps \
		-o $hapdir/genotypes_chr${chr}_hapgen \
		-n 200 0 \
		-dl $dummyDL 0 0 0 \
		-no_haps_output
done

#chr5=`wc -l $outdir/genotypes_chr5_hapgen.controls.gen |cut -d " " -f 1` - doent work?
#chr5=`wc -l $hapdir/genotypes_chr5_hapgen.controls.gen 
#chr7=`wc -l $hapdir/genotypes_chr7_hapgen.controls.gen  
#chr11=`wc -l $hapdir/genotypes_chr11_hapgen.controls.gen  
#echo "Chr5 SNPs: $chr5 -- Chr7 SNPs: $chr7 -- Chr11 SNPs: $chr11" >>  $hapdir/VARIANT_NUMBER.txt


cd $hapdir

# convert to plink format and prune SNPs
for chr in `seq 1 22`; do
	plink --data genotypes_chr${chr}_hapgen.controls \
		--oxford-single-chr $chr \
		--make-bed \
		--out genotypes_chr${chr}_hapgen.controls


	plink --bfile genotypes_chr${chr}_hapgen.controls \
		--indep-pairwise 50kb 1 0.5 \
		--out genotypes_chr${chr}_hapgen.controls

	plink --bfile genotypes_chr${chr}_hapgen.controls \
		--extract genotypes_chr${chr}_hapgen.controls.prune.in \
		--make-bed \
		--out genotypes_chr${chr}_hapgen.controls.pruned
	echo "genotypes_chr${chr}_hapgen.controls.pruned" >> file_list
done

# Merge chromsome-wide files into a single, genome-wide file
plink --merge-list file_list --make-bed --allow-no-sex --out genotypes_genome_hapgen.controls


#for chr in `seq 1 22`; do
#	plink --bfile genotypes_chr${chr}_hapgen.controls.pruned \
#		--exclude genotypes_genome_hapgen.controls-merge.missnp \
#		--make-bed \
#		--out genotypes_chr${chr}_hapgen.controls.nodup
#	echo "genotypes_chr${chr}_hapgen.controls.nodup" >> file_list_nodup
#done

# compute kinship
plink --bfile genotypes_genome_hapgen.controls\
	--make-rel square \
	--out genotypes_genome_hapgen.controls.grm
