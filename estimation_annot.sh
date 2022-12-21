SGE_TASK_ID=1
source /u/local/Modules/default/init/modules.sh
module load anaconda3

#genName="filter4"
#genPath=/u/project/sgss/UKBB/data/cal/
genName="small"
genPath=/u/home/a/aanand2/rGxG/GENIE_gxg/example/
# pcPath='/u/scratch/b/boyang19/Angela/data/traits/'
code=/u/home/a/aanand2/PGS/build_dev/GENIE_GxG
annotPath=/u/scratch/a/aanand2/PGS_test/annot_temp/
metacode=/u/home/a/aanand2/PGS/pgsFcreate.py

gxgbin=0
# Type="unweighted"
Type="Ridge"
# c_index=0
# for c_index in "${!carrier[@]}"; do
phenPath="/u/scratch/a/aanand2/no_epistasis_pheno/sim_geno/cau_ratio_0.01_h2_0.25_pheno/"
# cov="/u/scratch/b/boyang19/Angela/data/traits/${covarites[$c_index]}"
sigfile="/u/scratch/a/aanand2/PGS_test/sig_cand_mkup/"
out=/u/scratch/a/aanand2/PGS_test/
gen=${genPath}/${genName}
# annot=/u/scratch/b/boyang19/Angela/data/annot/rm_gene/${genes[$c_index]}.gene.annot

if [ ! -d "$out" ]; then
  mkdir ${out}
fi
#var=$(awk "NR==${SGE_TASK_ID}" ${sigfile}/1_sig_summary_all.txt)
var=$(awk "NR==${SGE_TASK_ID}" ${sigfile}/1_sig_test.txt)
IFS=' ' read phen Index Chr <<< $var

if [ ! -f "$out/${outfile}" ]; then
  outfile=${Index}.${Chr}.${phen}.${gxgbin}.${Type}.out.txt
  target=${annotPath}${phen}-${Index}-${Chr}.${Type}.target
  annot=${annotPath}${phen}-${Index}-${Chr}.${Type}.annot
  python ${metacode} ${annotPath} ${Type}
  $code  -g $gen -p ${phenPath}1.pheno -c ${annotPath}${phen}-${Index}-${Chr}.covar -gxgbin ${gxgbin} -snp -1 -k 10  -jn 100 -tp $target -o $out/${outfile}  -annot $annot
  rm ${annotPath}${phen}-${Index}-${Chr}.${Type}.covar
  rm ${annotPath}${phen}-${Index}-${Chr}.${Type}.target
  rm ${annotPath}${phen}-${Index}-${Chr}.${Type}.annot
fi
