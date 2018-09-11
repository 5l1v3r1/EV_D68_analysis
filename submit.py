import glob,os

#samples = glob.glob('MIC4106/*')
#samples = glob.glob('MIC3108/*')
samples = glob.glob('scilife_data/*/*')
gen_outdir = "../../../mapped_data"
cwd = os.getcwd()

for sample in samples:
	sname = sample.split('/')[-1]
	os.chdir(sample)
	call = 'nohup snakemake --snakefile ' + cwd + '/SVVC/Snakefile ' \
			+ gen_outdir + '/' + sname + '/minor.fasta'\
			+' --jobs 32 --cluster "sbatch -t 05:59:00 --mem=8G" 1>log'+sname + ' &'

	print(call)
	os.system(call)
	os.chdir(cwd)
