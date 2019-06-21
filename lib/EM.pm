#!/usr/local/bin/perl
##############################################################################
# PACKAGE NAME:	EM.pm
# DESCRIPTION:	
# PARAMETERS:	none
# OUTPUT:	
# TABLES:       
#
# DATE WRITTEN: 2013-05-02
# WRITTEN BY:    John L Freeman
#
# REVISION HISTORY: 
# REVISION DATE		REVISED BY	DESCRIPTION 
# ------- ----------	--------------	-------------------------------------
#
#       COPYRIGHT (C) 2013
#               ALL RIGHTS RESERVED        
##############################################################################
package EM;
use strict;    # always
use warnings;  # or else

#XRSQD^P XRSQP^D XRSDP^Q XRQDP^S XSQDP^R RSQDP^X

#LOCI_COMBOS:=RQ XR DP QS XRQS XRQSDP
our %h_loci = (
	ACB => "A", 
	CB => "C,B", 
	DP => "DPB1,DPA1",
	QS => "DQB1,DQA1",
	XR => "DRBX,DRB1",
	ACBXR => "A,C,B,DRBX,DRB1",
	XRQS => "DRBX,DRB1,DQB1,DQA1",
	XRQSDP => "DRBX,DRB1,DQB1,DQA1,DPB1,DPA1",
        XRQ    => "DQB1",
	RQS   => "DRB1,DQB1,DQA1",
	XD    => "DRBX,DPB1",       
	RD    => "DRB1,DPB1",
	QD    => "DQB1,DPA1",
	SD    => "DQA1,DPB1",
	XDP   => "DRBX,DPB1,DPA1",
	RDP   => "DRB1,DPB1,DPA1",
	QDP   => "DQB1,DPB1,DPA1",
	SDP   => "DQA1,DPB1,DPA1",
	XRD   => "DRBX,DRB1,DPB1",
	QSD   => "DQB1,DQA1,DPB1",
	XRDP  => "DRBX,DRB1,DPB1,DPA1",
	QSDP  => "DQB1,DQA1,DPB1,DPA1",
	RQSDP   => "DRB1,DQB1,DQA1,DPB1,DPA1",
	ACBXRQ  => "A,C,B,DRBX,DRB1,DQB1",
	ACBXRQP  => "DPB1",
	ACBXRQPD  => "DPA1",
	ACBXRQPDS  => "DQA1",
	ACBXRS    => "DQB1",
	ACBXRSP   => "DPB1",
	ACBXRSPD  => "DPA1",
	ACBXRSPDQ => "DQA1",
);


# our %h_loci_pos = (
# 	A => 0,
# 	C => 1,
# 	B => 2,
# 	DRBX => 3,
# 	DRB1 => 4,
# 	DQB1 => 5
# );


our %h_loci_pos = (
	A => 0,
	C => 1,
	B => 2,
	DRBX => 3,
	DRB1 => 4,
	DQA1 => 5,
	DQB1 => 6,
	DPA1 => 7,
	DPB1 => 8,
);

our %h_block = (
	ACBXR   => "XR",
	ACBXRQ   => "XRQ",
	XRACB 	=> "ACB",
	XRQACB  => "ACB",
	XRQS 	=> "QS",
	QSXR 	=> "XR",
	DPXRQS 	=> "XRQS",
	XRQSDP	=> "DP",
	DP 		=> "",
	QS 		=> "",
	XR 		=> "",
	ACB => "CB",
	XRQ => "XR",
	QXR => "",
        RD => "",
	ACBXRQP => "ACBXRQ",
	ACBXRQPD => "ACBXRQP",
	ACBXRQPDS => "ACBXRQPD",
	PACBXRQ => "",
	CB => "",
	ACBXRS 		=> "ACBXR",
	ACBXRSP 	=> "ACBXRS",
	ACBXRSPD 	=> "ACBXRSP",
	ACBXRSPDQ 	=> "ACBXRSPD",
);

our $impute_threshold = 0.01;
