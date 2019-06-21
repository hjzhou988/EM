#!/usr/bin/perl
############################################################################
# SCRIPT NAME:  make_GL_block
# DESCRIPTION:  creates GLID file and donor input files into input for
#               NEMO EM algorithm
# DATE WRITTEN: July 2, 2011
# WRITTEN BY:   Loren Gragert
#
# REVISION HISTORY:
# REVISION DATE         REVISED BY      DESCRIPTION
# ------- ----------    --------------  -------------------------------------
# 1.1	2013-05-02	jfreeman	use getopt
# 1.1.1 2019-6-20       Huajun Zhou     Get rid of opt_g since it was never used
#                                       Make configfile (association file) DRBX_DRB1_assoc.cfg a command-line argument
#                                       Separate input directory and output directory (change s_dir to i_dir and o_dir)
#                                       Deleted loadFreqs, since it was not used.
#
##############################################################################
use strict; # always
use warnings; # always
#use FindBin;
#
#use lib "$FindBin::Bin/../lib";
#use lib "/vol/bio/users/mhalagan/9LocusAnalysis/lib";
use lib "lib/";
use EM;
use Data::Dumper;
use Getopt::Std;

my $blank_geno;
my %h_blankhaps;
my %h_dptyped;
my %h_blank_dpa;
my %h_a;
my %h_haps;
my $f_blank_cutoff = 0.000001;

my %h_blank;
my %h_dp_pull;
my %h_dp_glid;
my %h_glstring; # stores GLstring per GLID
my %h_glstring2; # stores GLstring per GLID
my $GLID_counter = 1000000;
my %glstring_DRBX; # stores GLstring per GLID
my %glstring_DRB1;
my %glstring_DRBX_new;
my %glstring_DRB1_new;
my %GLID_DRBX; # GLID make from block for each donor
my %GLID_DRB1; # GLID make from block for each donor
my %genos_id_DRBX;
my %genos_id_DRB1;
my %DRBX_DRB1_assoc;
my %ID; # stores GLIDs per ID
my %GLID_pile;
my %GLIDlist_count; # stores count for set of GLIDs

our ($opt_c,$opt_i,$opt_o, $opt_r, $opt_f);
getopt("c:i:o:r:f:");

#my $config_dir = "/home/ec2-user/efs/9Locus/deidentified/configs";

die "Usage: make_GL_block -c loci_combo -r race -i input directory -o output directory -f association file" unless ($opt_c && $opt_r && $opt_i && $opt_o && $opt_f);
print "make_GL_block -c $opt_c -r $opt_r\n";

my $i_dir     = $opt_i;
my $o_dir     = $opt_o;
my $em_folder = $opt_c;
my $race      = $opt_r;
my $configfile = $opt_f;

print STDERR "$o_dir\n";

my @loci = split /,/,$EM::h_loci{$opt_c};
print "loci = @loci\n";

my %loc_positions;
foreach my $s_locus (@loci) {
    $loc_positions{$EM::h_loci_pos{$s_locus}}++;
}
print "loc_positions = " . join (",", keys (%loc_positions)) . "\n";

print STDERR "opt_c = $opt_c\n";



if ($opt_c eq "XRQSDP") {

    #loadPull($race);
   # loadPull($race);
# no need to read GLID file
# need to read two imputation files
    my $rh_imputed_genos = read_impute($race, $opt_c);
    my $rh_GLID_block_ACB = block($race, $opt_c, $rh_imputed_genos);
# hack
    $rh_imputed_genos = read_impute($race, "DPXRQS");
    my $rh_GLID_block_XRQ = block($race, $opt_c, $rh_imputed_genos);
    new_glid($race);
    donor_2Blocks($race, $rh_GLID_block_ACB, $rh_GLID_block_XRQ);
    new_pull($race);
    nemo ($race);

} elsif ($opt_c eq "XRQS") {
   # loadPull($race);
# no need to read GLID file
# need to read two imputation files
    my $rh_imputed_genos = read_impute($race, $opt_c);
    my $rh_GLID_block_ACB = block($race, $opt_c, $rh_imputed_genos,"QS");
# hack
    $rh_imputed_genos = read_impute($race, "QSXR");
    my $rh_GLID_block_XRQ = block($race, $opt_c, $rh_imputed_genos);
    new_glid($race);
    donor_2Blocks($race, $rh_GLID_block_ACB, $rh_GLID_block_XRQ);
    new_pull($race);
    nemo ($race);

} elsif ($opt_c eq "ACBXR") {

# no need to read GLID file
# need to read two imputation files
    my $rh_imputed_genos = read_impute($race, $opt_c);
    my $rh_GLID_block_ACB = block($race, $opt_c, $rh_imputed_genos);
# hack
    $rh_imputed_genos = read_impute($race, "XRACB");
    my $rh_GLID_block_XRQ = block($race, $opt_c, $rh_imputed_genos);
    new_glid($race);
    donor_2Blocks($race, $rh_GLID_block_ACB, $rh_GLID_block_XRQ);
    new_pull($race);
    nemo ($race);

} elsif ($opt_c eq "ACBXRQ") {

# no need to read GLID file
# need to read two imputation files
    my $rh_imputed_genos = read_impute($race, $opt_c);
    my $rh_GLID_block_ACB = block($race, $opt_c, $rh_imputed_genos);
# hack
    $rh_imputed_genos = read_impute($race, "XRQACB");
    my $rh_GLID_block_XRQ = block($race, $opt_c, $rh_imputed_genos);
    new_glid($race);
    donor_2Blocks($race, $rh_GLID_block_ACB, $rh_GLID_block_XRQ);
    new_pull($race);
    nemo ($race);

} elsif ($opt_c eq "XR") {
   # print STDERR "XR Blocking!\n";
    read_glid_DRB($race);
    donor($race);
    xr($race);
    xr_block();
    new_glid_DRB($race);
    new_pull_DRB($race);
    nemo ($race);
} else { # "CB","ACB","XRQ",etc..
    #loadPull($race) if $opt_c eq "DP" || $opt_c eq "QS";
    read_glid($race);
    my $rh_imputed_genos = read_impute($race, $opt_c);
    my $rh_GLID_block = block($race, $opt_c, $rh_imputed_genos);
    new_glid($race);
    donor($race, $rh_GLID_block);
    new_pull($race);
    nemo ($race);
}

################################################################################
# read GLID file
################################################################################
sub loadPull {

    my $race = shift;
    my $donor_file = $i_dir."/$race.pull";
    print "<pull_file = $donor_file\n";

    #my($i,$j)             = $s_combo eq "DP" ? (7,8) : (5,6);
    #my($i_blank,$j_blank) = $s_combo eq "DP" ? (11111111111,777777777) : (55555555,999999999);
    open (my $fh_donor,"<",$donor_file);
    while (<$fh_donor>) {
        chomp;
        my ($id,@glids) = split /:/;

        if($glids[7] == 11111111111){
          #  $h_blank_dpa{$id}++;
        }
        if ($glids[8] == 777777777 || $glids[6] == 999999999){
           # $h_blank{$id}++;
           # next;
        }

        foreach(@glids){ $h_dptyped{$_}++;}

    }
    close $fh_donor;

    my $num_blanks = scalar keys %h_blank;
    print "num blanks: $num_blanks\n";
}


################################################################################
# read GLID file
################################################################################
sub read_glid {
    my $race = shift;
    my $glid_file = $i_dir."/$race.glid";
    print "<glid_file = $glid_file\n";
    open (my $fh_glidin,"<",$glid_file);
    while (<$fh_glidin>) {
    	chomp;
    	my ($glid,$s_glstring) = split /\;/;
        #next if !defined $h_dptyped{$glid};
    	foreach my $loc (@loci) {
    	    my $loc_length = length ($loc);
    	    if ((substr($s_glstring,0,$loc_length) eq $loc)) {
    		  $h_glstring{$glid} = $s_glstring;
    	    }
    	}
    }
    close $fh_glidin;
    print "read_glid: $glid_file\n";
    print "# glstring = " . (scalar keys %h_glstring) . "\n";

#print Dumper(\%glstring);
}
################################################################################
# special case for DRBX, maybe this can be folded into read_glid
################################################################################
sub read_glid_DRB {
    my $race = shift;
    my $glid_file = $i_dir."/$race.glid";
    print "<glid_file = $glid_file\n";
    open (GLID,$glid_file);
    while (<GLID>) {
      chomp;
	 my $cnt =0;
      my ($glid,$glstring) = split /\;/;
      # only include GLIDs for relevant loci
      if ((substr($glstring,0,4) eq "DRB3") ||
          (substr($glstring,0,4) eq "DRB4") ||
          (substr($glstring,0,4) eq "DRB5") ||
          (substr($glstring,0,4) eq "DRBX")) {
        $glstring_DRBX{$glid} = $glstring;
		$cnt++;
      }
      if ((substr($glstring,0,4) eq "DRB1")) {
        $glstring_DRB1{$glid} = $glstring;
		$cnt++;
      }
	  # my $loc = substr($glstring,0,4);
	  # if($cnt == 0){

		# #print STDERR "$loc $glid $glstring\n";
	  # }
      #
    }
    close GLID;

}

################################################################################
# print new GLID file
################################################################################
sub new_glid {
    my $race = shift;
    my $new_glid_file = $o_dir."/$em_folder/$race.glid";
    print ">new_glid_file = $new_glid_file\n";
    open (NEWGLID,">$new_glid_file") or die "No can haz open $new_glid_file: $!";
    print "# new glstring = " . (scalar keys %h_glstring) . "\n";
    print "new glid: $new_glid_file\n";
    foreach my $glid (sort keys %h_glstring) {
	print NEWGLID "$glid;$h_glstring{$glid}\n";
    }
    close NEWGLID;
}

################################################################################
# read donor file
################################################################################
sub donor {
    my ($race, $rh_GLID_block) = @_;
    my $donor_file = $i_dir."/$race.pull";
    print "<donor_file = $donor_file\n";

    open (my $fh_donorin,"<",$donor_file);
    while (<$fh_donorin>) {
    	chomp;
    	my ($id,@glids) = split /:/;

      #  next if defined $h_blank{$id};

    	my $glid_index = 0;
    	my @new_glids;
    	my $skip = 0;
    	foreach my $glid (@glids) {
    	    if (exists $loc_positions{$glid_index}) {
    		  push @new_glids,$glid;
    	    }
    	    $glid_index++;
    	}

    	my $glidlist = join ":",@new_glids;
    	if ($opt_c eq "XR" || $opt_c eq "CB" || $opt_c eq "DP" || $opt_c eq "QS" || $opt_c eq "RD") {
    		$glidlist = join ":",$new_glids[1],$new_glids[0];
            $glidlist = join ":",$new_glids[0],$new_glids[1] if $opt_c eq "XR";
            $glidlist = join ":",$new_glids[0],$new_glids[1] if $opt_c eq "DP" || $opt_c eq "RD";
    	} else {
            next if !defined $$rh_GLID_block{$id};
            if($opt_c eq "ACB"){
                $glidlist = "$glidlist:$$rh_GLID_block{$id}";
            }else{
                $glidlist = "$$rh_GLID_block{$id}:$glidlist";
            }
    	}

    	$ID{$id} = $glidlist;
    	$GLIDlist_count{$glidlist}++;
    	push @{$GLID_pile{$glidlist}}, $id
    }
    close $fh_donorin;
}

################################################################################
# read ACBXRQ donor file
################################################################################
sub donor_2Blocks {
    my ($race, $rh_GLID_block_ACB, $rh_GLID_block_XRQ) = @_;

    my $donor_file = $o_dir."/$EM::h_block{$opt_c}/$race.pull";
    print "<donor_file = $donor_file\n";

    open (DONOR,$donor_file);
    while (<DONOR>) {
    	chomp;
    	my ($id,@glids) = split /:/;

        next unless ($$rh_GLID_block_ACB{$id});
        next unless ($$rh_GLID_block_XRQ{$id});
        #next if defined $h_blank{$id};

    	my $glidlist = $$rh_GLID_block_XRQ{$id} . ":" . $$rh_GLID_block_ACB{$id};

    	$ID{$id} = $glidlist;
    	$GLIDlist_count{$glidlist}++;
    }
    close DONOR;
    print "# donor_2Blocks = " . (scalar keys %ID) . "\n";
}
################################################################################
# write new donor pull file
################################################################################
sub new_pull {
    my $race = shift;
    my $new_pull_file = $o_dir."/$em_folder/$race.pull";
    print ">new_pull_file = $new_pull_file\n";
    print "# new_pull donors = " . (scalar keys %ID) . "\n";
    open (NEWPULL,">$new_pull_file") or die "No can haz open $new_pull_file: $!";
    foreach my $id (sort keys %ID) {
	print NEWPULL "$id:$ID{$id}\n";
    }
    close NEWPULL;
}
################################################################################
# write EM input file
################################################################################
sub nemo {
    my $race = shift;
    my $nemo_file = $o_dir."/$em_folder/$race.nemo";
    print ">nemo_file = $nemo_file\n";
    open (NEMO,">$nemo_file") or die "No can haz open $nemo_file: $!";
    foreach my $glidlist (sort keys %GLIDlist_count) {
    #print "$race,$glidlist,$GLIDlist_count{$glidlist}\n";
	print NEMO "$glidlist,$GLIDlist_count{$glidlist}\n";
    }
    close NEMO;
}

################################################################################
# Generate two locus blocks
################################################################################
sub block {
    my ($race, $combo, $rh_imputed_genos, $s_block) = @_;

    print STDERR "blocking $combo - ";
    return unless $EM::h_block{$combo};
    print STDERR $EM::h_block{$combo},"\n";

    my %GLstrings_block;
    my %GLID_block;
# create new GLIDs for genotype lists

    print "Blank Block: " . $race . " " . $combo . " " . $s_block,"\n" if defined $s_block && $s_block =~ /\S/;
    my $num_blank_subs = 0;
    my $n_blank_observed = 0;
    foreach my $id (keys %{$rh_imputed_genos}) {

       # if($n_blank_observed == 1 && defined $h_blank{$id} && defined $s_block && $s_block =~ /\S/){
       #     $GLID_block{$id} = $GLstrings_block{$blank_geno};
       # }else{

        #my $s_glstring;
        #if($n_blank_observed == 0 && defined $h_blank{$id} && defined $s_block && $s_block =~ /\S/){
        #    $s_glstring       = $blank_geno;
        #    $n_blank_observed = 1;
        #}else{

        my $s_glstring = join "|",keys %{$$rh_imputed_genos{$id}};;
        #}
           # if GLstring is new, create new GLID
        if (!exists $GLstrings_block{$s_glstring}) {
            $h_glstring{$GLID_counter} = $s_glstring;
            $GLID_block{$id} = $GLID_counter;
            $GLstrings_block{$s_glstring} = $GLID_counter;
            #last if $GLID_counter > 1000027;
            $GLID_counter++;
        }
        else {
            $GLID_block{$id} = $GLstrings_block{$s_glstring};
        }


    } # end foreach imputed_genos
    $GLID_counter = 2000000;
    print "# GLID_block = " . (scalar keys %GLID_block) . "\n";
    print "# Blank subjects = " . $num_blank_subs . "\n";
    return \%GLID_block;
}

################################################################################
#  read imputation file
################################################################################
sub read_impute {
    my ($race, $combo) = @_;
    return unless $EM::h_block{$combo};

# get block genotypes from imputation results
    my %imputed_genos;
    my %h_glstring;
    my $impute_file = $o_dir."/$EM::h_block{$combo}/${race}.impute";
    print "<impute_file = $impute_file\n";
    my %freqnorm_highest;
    open (IMP,$impute_file) or die "No can haz open $impute_file: $!\n";;

    while(<IMP>) {
	chomp;

	my ($i_nmdpid, $i_rank, $s_h1,$s_h2, $f_freqnorm) = split /,/;
    next if defined $h_blank{$i_nmdpid};

    foreach($i_nmdpid, $i_rank, $s_h1,$s_h2, $f_freqnorm){
        $_ =~ s/ //g;
    }


    #print "impute1: ",$i_nmdpid,"\t",$f_freqnorm," | ",$_,"\n";
# skip less probable genotypes

    # *** No threshold *** #


# skip less probable genotypes
    if (!exists $freqnorm_highest{$i_nmdpid}) {
         $freqnorm_highest{$i_nmdpid} = 0;
    }
    if ($f_freqnorm > $freqnorm_highest{$i_nmdpid}) {
        $freqnorm_highest{$i_nmdpid} = $f_freqnorm;
    }
    if (($f_freqnorm < $EM::impute_threshold) &&
        $freqnorm_highest{$i_nmdpid} > $EM::impute_threshold) {
        next;
    }

       # sort genotypes
	if ($s_h1 gt $s_h2) {
	    my $temp = $s_h1;
	    $s_h1 = $s_h2;
	    $s_h2 = $temp;
	}

    # if($combo eq "DP" || $combo eq "QS" && defined $h_blank_dpa{$i_nmdpid}){
    #     my($s_a1,$s_b1) = split(/\~/,$s_h1);
    #     my($s_a2,$s_b2) = split(/\~/,$s_h2);

    #     foreach my $s_aa1 (keys %h_a){
    #         foreach my $s_aa2 (keys %h_a){
    #             my $h1 = join("~",$s_aa1,$s_a1);
    #             next if !defined $h_haps{$h1};
    #             my $h2 = join("~",$s_aa2,$s_a2);
    #             next if !defined $h_haps{$h2};
    #             $imputed_genos{$i_nmdpid}{join(",",sort $h1,$h2)} = 1;
    #         }
    #     }

    # }else{

        my ($s_genotype) = join ",",$s_h1,$s_h2;
        $imputed_genos{$i_nmdpid}{$s_genotype} = 1;

    #}


     #print "impute2: ",$i_nmdpid,"\t",$s_genotype,"\n";
    }



    close IMP;
    print "# imp = " . (scalar keys %imputed_genos) . "\n";
#print Dumper(\%imputed_genos);

    return \%imputed_genos;
}
################################################################################
# DRBX~DRB1
################################################################################
sub xr{

  my $race = shift;

#  load DRBX~DRB1 association rules file
    print "<assoc_file = $configfile\n";
    my $drx = $configfile;
    open (ASSOC,"$drx") || die "DRBX_DRB1_assoc.cfg";
    while (<ASSOC>) {
    chomp;
    my ($drb1,$drbx) = split /,/;
    if ($drb1 eq "DRB1") { next; }
    $DRBX_DRB1_assoc{$drb1} = $drbx;
 print STDERR "$drb1 $drbx $DRBX_DRB1_assoc{$drb1}\n";
    }
    close ASSOC;


# make all possible combo of DRBX genos with DRB1 genos
    #my $restricted_filename = $s_dir."/restricted/$race.DRBXDRB1.restricted";
   # print ">restricted_file = $restricted_filename\n";
    #open (my $fh_restricted,">",$restricted_filename);

    my %h_geno_glid;

    my %h_temp_drb1;
    my %h_temp_drbx;

    if($race ne "CAU"){
        my $dr_cnt = 0;
        my @a_geno_gl;
        foreach my $glidlist (keys %GLID_pile) {

            my $rhp = xr_id($race, $glidlist);
            my $n_restricted_genos_DRBX = scalar @{$$rhp{genos_DRBX_restricted}};
            my $n_restricted_genos_DRB1 = scalar @{$$rhp{genos_DRB1_restricted}};

            if ($n_restricted_genos_DRB1 > 0) {
                foreach my $geno (@{$$rhp{genos_DRBX_restricted}}) {
                    my $geno_glid;
                    if(!defined $h_geno_glid{$geno}){
                        $geno_glid             = $dr_cnt;
                        $h_geno_glid{$geno}    = $geno_glid;
                        push(@a_geno_gl,$geno);
                        $dr_cnt++;
                    }else{
                        $geno_glid = $h_geno_glid{$geno};
                    }
                    $h_temp_drbx{$glidlist}{$geno_glid}++;
                }
                foreach my $geno (@{$$rhp{genos_DRB1_restricted}}) {
                    my $geno_glid;
                    if(!defined $h_geno_glid{$geno}){
                        $geno_glid             = $dr_cnt;
                        $h_geno_glid{$geno}    = $geno_glid;
                        push(@a_geno_gl,$geno);
                        $dr_cnt++;
                    }else{
                        $geno_glid = $h_geno_glid{$geno};
                    }
                    $h_temp_drb1{$glidlist}{$geno_glid}++;
                }
            }
            else {
                foreach my $geno (@{$$rhp{genos_DRBX_all}}) {
                    my $geno_glid;
                    if(!defined $h_geno_glid{$geno}){
                        $geno_glid             = $dr_cnt;
                        $h_geno_glid{$geno}    = $geno_glid;
                        push(@a_geno_gl,$geno);
                        $dr_cnt++;
                    }else{
                        $geno_glid = $h_geno_glid{$geno};
                    }
                    $h_temp_drbx{$glidlist}{$geno_glid}++;
                }
                foreach my $geno (@{$$rhp{genos_DRB1_all}}) {
                    my $geno_glid;
                    if(!defined $h_geno_glid{$geno}){
                        $geno_glid             = $dr_cnt;
                        $h_geno_glid{$geno}    = $geno_glid;
                        push(@a_geno_gl,$geno);
                        $dr_cnt++;
                    }else{
                        $geno_glid = $h_geno_glid{$geno};
                    }
                    $h_temp_drb1{$glidlist}{$geno_glid}++;
                }
            }
        }

        %h_geno_glid = ();
        foreach my $glidlist (keys %GLID_pile) {
            foreach my $id (@{$GLID_pile{$glidlist}}) {
                $genos_id_DRBX{$id} = join "|",sort map{ $a_geno_gl[$_] } (keys %{$h_temp_drbx{$glidlist}});
                $genos_id_DRB1{$id} = join "|",sort map{ $a_geno_gl[$_] } (keys %{$h_temp_drb1{$glidlist}});
            }
        }
        @a_geno_gl   = ();
        %h_temp_drb1 = ();
        %h_temp_drbx = ();

    }else{
        my $dr_cnt = 0;
        my @a_geno_gl;
        my $num_glidlists = (scalar keys %GLID_pile);
        print STDERR "Num glidlist: $num_glidlists\n";
        foreach my $glidlist (keys %GLID_pile) {

            my $rhp = xr_id($race, $glidlist);
            my $n_restricted_genos_DRB1 = scalar @{$$rhp{genos_DRB1_restricted}};

            if ($n_restricted_genos_DRB1 > 0) {
                foreach my $geno (@{$$rhp{genos_DRBX_restricted}}) {
                    my $geno_glid;
                    if(!defined $h_geno_glid{$geno}){
                        $geno_glid             = $dr_cnt;
                        $h_geno_glid{$geno}    = $geno_glid;
                        push(@a_geno_gl,$geno);
                        $dr_cnt++;
                    }else{
                        $geno_glid = $h_geno_glid{$geno};
                    }
                    $h_temp_drbx{$glidlist}{$geno_glid}++;
                }
            }
            else {
                foreach my $geno (@{$$rhp{genos_DRBX_all}}) {
                    my $geno_glid;
                    if(!defined $h_geno_glid{$geno}){
                        $geno_glid             = $dr_cnt;
                        $h_geno_glid{$geno}    = $geno_glid;
                        push(@a_geno_gl,$geno);
                        $dr_cnt++;
                    }else{
                        $geno_glid = $h_geno_glid{$geno};
                    }
                    $h_temp_drbx{$glidlist}{$geno_glid}++;
                }
            }
        }
        print STDERR "Adding DRBX..\n";

        %h_geno_glid = ();
        foreach my $glidlist (keys %GLID_pile) {
            foreach my $id (@{$GLID_pile{$glidlist}}) {
                $genos_id_DRBX{$id} = join "|",sort map{ $a_geno_gl[$_] } (keys %{$h_temp_drbx{$glidlist}});
            }
        }
        @a_geno_gl   = ();
        %h_temp_drbx = ();
        print STDERR "Finished DRBX! On to DRB1..\n";

        $dr_cnt = 0;
        foreach my $glidlist (keys %GLID_pile) {

            my $rhp = xr_id($race, $glidlist);
            my $n_restricted_genos_DRBX = scalar @{$$rhp{genos_DRBX_restricted}};
            my $n_restricted_genos_DRB1 = scalar @{$$rhp{genos_DRB1_restricted}};

            if ($n_restricted_genos_DRB1 > 0) {
                foreach my $geno (@{$$rhp{genos_DRB1_restricted}}) {
                    my $geno_glid;
                    if(!defined $h_geno_glid{$geno}){
                        $geno_glid             = $dr_cnt;
                        $h_geno_glid{$geno}    = $geno_glid;
                        push(@a_geno_gl,$geno);
                        $dr_cnt++;
                    }else{
                        $geno_glid = $h_geno_glid{$geno};
                    }
                    $h_temp_drb1{$glidlist}{$geno_glid}++;
                }
            }
            else {
                foreach my $geno (@{$$rhp{genos_DRB1_all}}) {
                    my $geno_glid;
                    if(!defined $h_geno_glid{$geno}){
                        $geno_glid             = $dr_cnt;
                        $h_geno_glid{$geno}    = $geno_glid;
                        push(@a_geno_gl,$geno);
                        $dr_cnt++;
                    }else{
                        $geno_glid = $h_geno_glid{$geno};
                    }
                    $h_temp_drb1{$glidlist}{$geno_glid}++;
                }
            }
        }

        print STDERR "Adding DRB1..\n";

        %h_geno_glid = ();
        foreach my $glidlist (keys %GLID_pile) {
            foreach my $id (@{$GLID_pile{$glidlist}}) {
                foreach my $geno_glid (keys %{$h_temp_drb1{$glidlist}}){
                    $genos_id_DRB1{$id} = join "|",sort map{ $a_geno_gl[$_] } (keys %{$h_temp_drb1{$glidlist}});
                }
            }
        }
        %h_temp_drb1 = ();
        @a_geno_gl   = ();
        print STDERR "Finished with DRB1!\n";


    }


}

################################################################################
# DRBX~DRB1
################################################################################

sub xr_id {
    my ($race, $glidlist) = @_;
    my ($glid_DRBX,$glid_DRB1) = split /:/,$glidlist;
    my $glstring_DRBX = $glstring_DRBX{$glid_DRBX};
    my $glstring_DRB1 = $glstring_DRB1{$glid_DRB1};
    #print STDERR "GLIDS: $glidlist\n";
    my @genos_DRBX = split /\|/,$glstring_DRBX;
	if(!defined $glstring_DRBX || $glstring_DRBX !~ /\S/){
		#print STDERR "Err DRBX: $glid_DRBX\n";
	}
	if(!defined $glstring_DRB1 || $glstring_DRB1 !~ /\S/){
		#print STDERR "Err DRB1: $glid_DRB1\n";
	}
    my @genos_DRB1 = split /\|/,$glstring_DRB1;
    my @genos_DRBX_restricted;
    my @genos_DRB1_restricted;
    my @genos_DRBX_all;
    my @genos_DRB1_all;
    my @genos_DRBXDRB1_all;
    foreach my $geno_DRBX (@genos_DRBX) {
	my ($drbx_1, $drbx_2) = split /,/,$geno_DRBX;
    next if $geno_DRBX =~ /\*9999/;
	my $drbx_loc1 = substr($drbx_1,0,4);
	my $drbx_loc2 = substr($drbx_2,0,4);
	foreach my $geno_DRB1 (@genos_DRB1) {
	    my ($drb1_1,$drb1_2) = split /,/,$geno_DRB1;
	    my ($drb1_loc1,$drb1_typ1) = split /\*/,$drb1_1;
	    my ($drb1_loc2,$drb1_typ2) = split /\*/,$drb1_2;
	    my ($drb1_1_2dig) = substr ($drb1_typ1,0,2);
	    my ($drb1_2_2dig) = substr ($drb1_typ2,0,2);

# phase 1
	    my $geno_DRBXDRB1_phase1_1 = "$drbx_1~$drb1_1";
	    my $geno_DRBXDRB1_phase1_2 = "$drbx_2~$drb1_2";

	    my $geno_phase1 = $geno_DRBXDRB1_phase1_1 . "," .
		$geno_DRBXDRB1_phase1_2;
	    push @genos_DRBX_all,$geno_DRBX;
	    push @genos_DRB1_all,$geno_DRB1;
	    push @genos_DRBXDRB1_all,$geno_phase1;

# check for assoc rules
# print STDERR "Phase 1: DRB1 $drb1_1_2dig DRBX $drbx_loc1 DRB1 $drb1_2_2dig DRBX $drbx_loc2\n";
	    if (($DRBX_DRB1_assoc{$drb1_1_2dig} eq $drbx_loc1) &&
		    ($DRBX_DRB1_assoc{$drb1_2_2dig} eq $drbx_loc2)) {
		push @genos_DRBX_restricted,$geno_DRBX;
		push @genos_DRB1_restricted,$geno_DRB1;
# print STDERR "Accepted!\n";
	    }


# phase 2
	    my $geno_DRBXDRB1_phase2_1 = "$drbx_1~$drb1_2";
	    my $geno_DRBXDRB1_phase2_2 = "$drbx_2~$drb1_1";
	    my $geno_phase2 = $geno_DRBXDRB1_phase2_1 . "," .
		$geno_DRBXDRB1_phase2_2;
	    push @genos_DRBX_all,$geno_DRBX;
	    push @genos_DRB1_all,$geno_DRB1;
	    push @genos_DRBXDRB1_all,$geno_phase2;

# check for assoc rules
# print STDERR "Phase 2: DRB1 $drb1_2_2dig DRBX $drbx_loc1 DRB1 $drb1_1_2dig DRBX $drbx_loc2\n";
	    if (($DRBX_DRB1_assoc{$drb1_2_2dig} eq $drbx_loc1) &&
		    ($DRBX_DRB1_assoc{$drb1_1_2dig} eq $drbx_loc2)) {
		push @genos_DRBX_restricted,$geno_DRBX;
		push @genos_DRB1_restricted,$geno_DRB1;
# print STDERR "Accepted!\n";
	    }


	} # end foreach geno_DRB1DQB1
    } # end foreach geno_DRBX

my %h;
$h{genos_DRBX_all} = \@genos_DRBX_all;
$h{genos_DRB1_all} = \@genos_DRB1_all;
$h{genos_DRBXDRB1_all} = \@genos_DRBXDRB1_all;
$h{genos_DRBX_restricted} = \@genos_DRBX_restricted;
$h{genos_DRB1_restricted} = \@genos_DRB1_restricted;
return \%h;
}

sub xr_block {

    my %GLstrings_DRBX;
    my %GLstrings_DRB1;
    my $GLID_counter_DRBX = 2000000;
    my $GLID_counter_DRB1 = 3000000;
    # create new GLIDs for DRBX genotype lists
    foreach my $id (keys %genos_id_DRBX) {

      my $glstring = $genos_id_DRBX{$id};

      # if GLstring is new, create new GLID
      if (!exists $GLstrings_DRBX{$glstring}) {
        $glstring_DRBX_new{$GLID_counter_DRBX} = $glstring;
        $GLID_DRBX{$id} = $GLID_counter_DRBX;
        $GLstrings_DRBX{$glstring} = $GLID_counter_DRBX;
        $GLID_counter_DRBX++;
      }
      else {
        $GLID_DRBX{$id} = $GLstrings_DRBX{$glstring};
      }

      delete $genos_id_DRBX{$id};

    } # end foreach id genos_id_DRBX

    # create new GLIDs for DRB1 genotype lists
    foreach my $id (keys %genos_id_DRB1) {
      my $glstring = $genos_id_DRB1{$id};


      # if GLstring is new, create new GLID
      if (!exists $GLstrings_DRB1{$glstring}) {
        $glstring_DRB1_new{$GLID_counter_DRB1} = $glstring;
        $GLID_DRB1{$id} = $GLID_counter_DRB1;
        $GLstrings_DRB1{$glstring} = $GLID_counter_DRB1;
        $GLID_counter_DRB1++;
      }
      else {
        $GLID_DRB1{$id} = $GLstrings_DRB1{$glstring};
      }

      delete $genos_id_DRB1{$id};

    }

}
sub new_glid_DRB {
    my $race = shift;
    # print new DRBX~DRB1 GLID file
    my $new_glid_file = $o_dir."/$em_folder/$race.glid";
    print STDERR ">new_glid_file = $new_glid_file\n";
    open (NEWGLID,">$new_glid_file") or die "$new_glid_file";
    foreach my $glid (sort keys %glstring_DRB1_new) {
      print NEWGLID "$glid;$glstring_DRB1_new{$glid}\n";
    }
    foreach my $glid (sort keys %glstring_DRBX_new) {
      print NEWGLID "$glid;$glstring_DRBX_new{$glid}\n";
    }
    close NEWGLID;
}


    # print new donor pull file
sub new_pull_DRB {
    my $race = shift;
    #my %GLIDlist_count;
    %GLIDlist_count = ();
    my $new_pull_file = $o_dir."/$em_folder/$race.pull";
    print STDERR ">new_pull_file = $new_pull_file\n";
    open (NEWPULL,">$new_pull_file") or die "$new_pull_file";
    foreach my $id (sort keys %GLID_DRB1) {
      print NEWPULL "$id:$GLID_DRBX{$id}:$GLID_DRB1{$id}\n";
      my $glidlist = "$GLID_DRBX{$id}:$GLID_DRB1{$id}";
      $GLIDlist_count{$glidlist}++;
    }
    close NEWPULL;
}

__END__
