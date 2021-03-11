#!/usr/bin/perl

`mkdir bam_merge`;
`sed -n '1p' input_user/input_subclonality.txt > input_user/input_subclonality_colnames.txt`;

@clone_all = ("clone1", "clone2", "clone3", "clone4", "clone5");
foreach $clone_ind (@clone_all) {

    `grep "$clone_ind" input_user/input_subclonality.txt | wc -l > line.txt`;
    open (FILE, "line.txt");
	while (<FILE>){
	chomp;
	($x) = split("\t");
	print "$x\n";
	}
    close (FILE);

    if ($x > 0) {

		`grep "$clone_ind" input_user/input_subclonality.txt > input_user/input_subclonality.$clone_ind.pre.txt`;
		`cat input_user/input_subclonality_colnames.txt input_user/input_subclonality.$clone_ind.pre.txt > input_user/input_subclonality.$clone_ind.txt`;
		`rm input_user/input_subclonality.$clone_ind.pre.txt`; 
		`mkdir bam_merge/$clone_ind`;

    	open (FILE, "input_user/input_subclonality.$clone_ind.txt");
    	while (<FILE>) {
		chomp;
		($cell, $clone) = split("\t");

		print "cell: $cell\n";
		print "clone: $clone\n";
		print "---------\n";

		`cp bam/$cell.sc_pre_mono_sort_for_mark_uniq.bam bam_merge/$clone_ind`;
		`cp bam/$cell.sc_pre_mono_sort_for_mark_uniq.bam.bai bam_merge/$clone_ind`;

		}
		close (FILE);
		`samtools merge bam_merge/$clone_ind.merge.bam bam_merge/$clone_ind/*.bam`;
		`samtools index bam_merge/$clone_ind.merge.bam`;
		
    }

} 
