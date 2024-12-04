module load HGI/common/VerifyBamID/2.0.1

mkdir -p /lustre/scratch123/hgi/projects/birth_cohort_wes/qc/BiB/verifybamid/

find /lustre/scratch126/humgen/projects/birth_cohort_wes/BiB/crams/ -name "*.cram" > cramfiles.txt

parallel --dryrun -j 4 -a cramfiles.txt \
	VerifyBamID \
		--SVDPrefix /software/hgi/installs/verifyBamID2/VerifyBamID-2.0.1/resource/1000g.phase3.100k.b38.vcf.gz.dat \
		--Reference /lustre/scratch125/humgen/resources/ref/Homo_sapiens/HS38DH/hs38DH.fa  \
		--DisableSanityCheck \
		--BamFile {} \
		--Output /lustre/scratch123/hgi/projects/birth_cohort_wes/qc/BiB/verifybamid/{/.} \
	> verifybamid.cmds

wr add -f verifybamid.cmds -i "ip13.verifybamid.bib" --req_grp "verifybamid" --memory "750M" --cpus 1 --limit_grps "verifybamid:500"

rm cramfiles.txt verifybamid.cmds
