import os
from snakemake.shell import shell

#get var caller name
var=[]
if snakemake.params.var_cal=={"deepvariant"}:
        var.append('dv')
elif snakemake.params.var_cal=={"haplotypecaller"}:
        var.append('hc')

# Get sample list from config file
samps=[]
configfile = "config.yaml"
fh_c=open(configfile,'r')
for line in fh_c:
	if line.startswith(' ') and ':' in line:
		line=line.strip()
		s=line.split(':')[0]
		samps.append(s)

# cd into benchmarks dir
os.chdir('./benchmarks')

# For each sample, create an output file
for samp in samps:
	vc=var[0]
	outfile=samp+'_'+vc+'_total_time.txt'
	times=[]
	fh_out=open(outfile,'w')
	fh_out.write('sample'+'\t'+'seconds'+'\n')
	
	out_an = samp+'_time_per_analysis.txt'
	fh_outan=open(out_an,'w')
	fh_outan.write('analysis'+'\t'+'seconds'+'\n')

	
	# Then grab files that start with that sample, and add them to a list of times for each analysis
	for file in os.listdir():
		if file.startswith(samp):
			fh=open(file,'r')
			for line in fh:
				line=line.strip()
				if not line.startswith('s'):
					secs=line.split('\t')[0]
					times.append(float(secs))
					#write the analysis name and the time to the output file
					fh_outan.write(file.split('.')[1]+'\t'+str(secs)+'\n')
	
	#add up the total and write to out file
	total = float(sum(times))
	fh_out.write(samp+'\t'+str(total)+'\n')

