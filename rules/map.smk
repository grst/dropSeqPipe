"""Align the data with STAR."""

ruleorder: plot_knee_plot_whitelist > plot_knee_plot


#Which rules will be run on the host computer and not sent to nodes
localrules: multiqc_star, plot_yield, plot_knee_plot, plot_knee_plot_whitelist, extend_barcode


rule STAR_align:
	input:
		fq1="data/{sample}/trimmmed_repaired_R2.fastq.gz",
		index=lambda wildcards: star_index_prefix + '_' + str(samples.loc[wildcards.sample,'read_length']) + '/SA'
	output:
		temp('data/{sample}/Aligned.out.bam')
	log:
		'data/{sample}/Log.final.out'
	params:
		extra="""--outReadsUnmapped Fastx\
				--outSAMattributes NH HI NM AS MD\
			 	--outFilterMismatchNmax {}\
			 	--outFilterMismatchNoverLmax {}\
			 	--outFilterMismatchNoverReadLmax {}\
			 	--outFilterMatchNmin {}\
			 	--outFilterScoreMinOverLread {}\
			 	--outFilterMatchNminOverLread {}""".format(
				config['MAPPING']['STAR']['outFilterMismatchNmax'],
				config['MAPPING']['STAR']['outFilterMismatchNoverLmax'],
				config['MAPPING']['STAR']['outFilterMismatchNoverReadLmax'],
				config['MAPPING']['STAR']['outFilterMatchNmin'],
				config['MAPPING']['STAR']['outFilterMatchNminOverLread'],
				config['MAPPING']['STAR']['outFilterScoreMinOverLread'],),
		index=lambda wildcards: star_index_prefix + '_' + str(samples.loc[wildcards.sample,'read_length']) + '/'
	threads: 8
	wrapper:
		"0.22.0/bio/star/align"

rule multiqc_star:
	input:
		expand('data/{sample}/Log.final.out', sample=samples.index)
	output:
		html='reports/star.html'
	params: '-m star'
	wrapper:
		'0.21.0/bio/multiqc'


rule MergeBamAlignment:
	input:
		mapped='data/{sample}/Aligned.out.bam',
		R1_ref = "data/{sample}/trimmmed_repaired_R1.fastq.gz"
	output:
		temp('data/{sample}/Aligned.merged.bam')
	params:
		BC_start=config['FILTER']['cell-barcode']['start']-1,
		BC_end=config['FILTER']['cell-barcode']['end'],
		UMI_start=config['FILTER']['UMI-barcode']['start']-1,
		UMI_end=config['FILTER']['UMI-barcode']['end'],
		discard_secondary_alignements=True
	conda: '../envs/merge_bam.yaml'
	script:
		'../scripts/merge_bam.py'

rule extend_barcode:
	input:
		whitelist='barcodes.csv'
	output:
		barcode_ref='summary/barcode_ref.pkl',
		barcode_ext_ref='summary/barcode_ext_ref.pkl',
		barcode_mapping='summary/barcode_mapping.pkl'
	script:
		'../scripts/generate_extended_ref.py'

rule repair_barcodes:
	input:
		bam='data/{sample}/Aligned.merged.bam',
		barcode_ref='summary/barcode_ref.pkl',
		barcode_ext_ref='summary/barcode_ext_ref.pkl',
		barcode_mapping='summary/barcode_mapping.pkl'
	conda: '../envs/merge_bam.yaml'
	output:
		bam='data/{sample}/Aligned.repaired.bam',
		barcode_mapping_counts='data/{sample}/barcode_mapping_counts.pkl'
	script:
		'../scripts/repair_barcodes.py'

rule TagReadWithGeneExon:
	input:
		data='data/{sample}/Aligned.repaired.bam',
		refFlat='{}.refFlat'.format(annotation_prefix)
	params:
		memory=config['LOCAL']['memory'],
		temp_directory=config['LOCAL']['temp-directory']
	output:
		temp('data/{sample}/gene_exon_tagged.bam')
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && TagReadWithGeneExon -m {params.memory}\
		INPUT={input.data}\
		OUTPUT={output}\
		ANNOTATIONS_FILE={input.refFlat}\
		TAG=GE\
		CREATE_INDEX=true
		"""

rule bead_errors_metrics:
	input:
		'data/{sample}/gene_exon_tagged.bam'
	output:
		'data/{sample}/final.bam'
	log:
		out_stats='logs/dropseq_tools/{sample}_synthesis_stats.txt',
		summary='logs/dropseq_tools/{sample}_synthesis_stats_summary.txt',
	params:
		barcodes=lambda wildcards: int(samples.loc[wildcards.sample,'expected_cells']) * 2,
		memory =config['LOCAL']['memory'],
		temp_directory=config['LOCAL']['temp-directory'],
		SmartAdapter = config['FILTER']['5-prime-smart-adapter']
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && DetectBeadSynthesisErrors -m {params.memory}\
		INPUT={input}\
		OUTPUT={output}\
		OUTPUT_STATS={log.out_stats}\
		SUMMARY={log.summary}\
		NUM_BARCODES={params.barcodes}\
		PRIMER_SEQUENCE={params.SmartAdapter}
		"""

rule bam_hist:
	input:
		'data/{sample}/final.bam'
	params:
		memory=config['LOCAL']['memory'],
		temp_directory=config['LOCAL']['temp-directory']
	output:
		'logs/dropseq_tools/{sample}_hist_out_cell.txt'
	conda: '../envs/dropseq_tools.yaml'
	shell:
		"""export _JAVA_OPTIONS=-Djava.io.tmpdir={params.temp_directory} && BAMTagHistogram -m {params.memory}\
		TAG=XC\
		I={input}\
		READ_QUALITY=10\
		O={output}
		"""


rule plot_yield:
	input:
		R1_filtered=expand('logs/cutadapt/{sample}_R1.qc.txt', sample=samples.index),
		R2_filtered=expand('logs/cutadapt/{sample}_R2.qc.txt', sample=samples.index),
		repaired=expand('logs/bbmap/{sample}_repair.txt', sample=samples.index),
		STAR_output=expand('data/{sample}/Log.final.out', sample=samples.index),
	params:
		BC_length=config['FILTER']['cell-barcode']['end'] - config['FILTER']['cell-barcode']['start']+1,
		UMI_length=config['FILTER']['UMI-barcode']['end'] - config['FILTER']['UMI-barcode']['start']+1,
		sample_names=lambda wildcards: samples.index,
		batches=lambda wildcards: samples.loc[samples.index, 'batch']
	conda: '../envs/plots.yaml'
	output:
		pdf='plots/yield.pdf'
	script:
		'../scripts/plot_yield.R'

rule plot_knee_plot:
	input:
		'logs/dropseq_tools/{sample}_hist_out_cell.txt'
	params:
		cells=lambda wildcards: samples.loc[wildcards.sample,'expected_cells'],
		edit_distance=config['EXTRACTION']['UMI-edit-distance']
	conda: '../envs/plots.yaml'
	output:
		pdf='plots/knee_plots/{sample}_knee_plot.pdf'
	script:
		'../scripts/plot_knee_plot.R'

rule plot_knee_plot_whitelist:
	input:
		data='logs/dropseq_tools/{sample}_hist_out_cell.txt',
		barcodes='barcodes.csv'
	params:
		cells=lambda wildcards: samples.loc[wildcards.sample,'expected_cells']
	conda: '../envs/plots.yaml'
	output:
		pdf='plots/knee_plots/{sample}_knee_plot.pdf'
	script:
		'../scripts/plot_knee_plot.R'
