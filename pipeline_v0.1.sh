metafunkdir="/home/projects/pr_46704/people/antalb/metafunk_v0.1"
xqsub -V -A pr_46704 -W group_list=pr_46704 -d `pwd` -e MetaFunk.err -o MetaFunk.out -l nodes=1:ppn=16,mem=32gb,walltime=0:06:00:00 -N MetaFunk -de sh /home/projects/pr_46704/people/antalb/metafunk_v0.1/metafunk.sh ${metafunkdir}

cat MetaFunk.err
cat test_20180322/test1/run.log
li test_20180322/test/
