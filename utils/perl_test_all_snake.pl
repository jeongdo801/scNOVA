#!/usr/bin/perl

`sed -n '1p' input_user/strandphaser_output.txt > input_user/strandphaser_colnames.txt`;
`mkdir nucleosome_sampleA`;
`mkdir nucleosome_sampleB`;

@chrom_all = ("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX");

foreach $chrom_ind (@chrom_all) {

	`grep "$chrom_ind\t" input_user/strandphaser_output.txt > input_user/strandphaser_output.$chrom_ind.pre.txt`;
	`cat input_user/strandphaser_colnames.txt input_user/strandphaser_output.$chrom_ind.pre.txt > input_user/strandphaser_output.$chrom_ind.txt`;
	`rm input_user/strandphaser_output.$chrom_ind.pre.txt`;



	`mkdir bam/$chrom_ind.H1`;
	`mkdir bam/$chrom_ind.H2`;

	open (FILE, "input_user/strandphaser_output.$chrom_ind.txt");
	while (<FILE>) {
	chomp;
	($sample, $cell, $chrom, $start, $end, $class, $hap1.cis.simil, $hap1.trans.simil, $hap2.cis.simil, $hap2.trans.simil) = split("\t");
    $cell=~s/.bam//;
	if ($chrom eq $chrom_ind){

		print "cell: $cell\n";
		print "chrom: $chrom\n";
		print "start: $start\n";
		print "end: $end\n";
		print "class: $class\n";
		print "---------\n";

		`samtools view -b bam/$cell.sc_pre_mono_sort_for_mark_uniq.bam.W.bam "$chrom:$start-$end" > bam/$cell.sc_pre_mono_sort_for_mark_uniq.bam.W.$chrom.$start.$end.seg.bam`;
		`samtools view -b bam/$cell.sc_pre_mono_sort_for_mark_uniq.bam.C.bam "$chrom:$start-$end" > bam/$cell.sc_pre_mono_sort_for_mark_uniq.bam.C.$chrom.$start.$end.seg.bam`;

		if ($class eq "WC") {
			`mv bam/$cell.sc_pre_mono_sort_for_mark_uniq.bam.W.$chrom.$start.$end.seg.bam bam/$chrom_ind.H1/`;
			`mv bam/$cell.sc_pre_mono_sort_for_mark_uniq.bam.C.$chrom.$start.$end.seg.bam bam/$chrom_ind.H2/`;
		} 

		if ($class eq "CW") {
			`mv bam/$cell.sc_pre_mono_sort_for_mark_uniq.bam.C.$chrom.$start.$end.seg.bam bam/$chrom_ind.H1/`;
			`mv bam/$cell.sc_pre_mono_sort_for_mark_uniq.bam.W.$chrom.$start.$end.seg.bam bam/$chrom_ind.H2/`;
		}

	}
	}	
	close (FILE);

`samtools merge bam/$chrom_ind.H1.bam bam/$chrom_ind.H1/*.seg.bam`;
`samtools index bam/$chrom_ind.H1.bam`;
`rm -r bam/$chrom_ind.H1`;

`samtools merge bam/$chrom_ind.H2.bam bam/$chrom_ind.H2/*.seg.bam`;
`samtools index bam/$chrom_ind.H2.bam`;
`rm -r bam/$chrom_ind.H2`;

close (Chr_FILE);

}

`samtools merge nucleosome_sampleA/result.H1.bam bam/*.H1.bam`;
`samtools merge nucleosome_sampleB/result.H2.bam bam/*.H2.bam`;
`samtools index nucleosome_sampleA/result.H1.bam`;
`samtools index nucleosome_sampleB/result.H2.bam`;
`cp input_user/strandphaser_output.txt input_user/strandphaser_output_copy.txt`;
