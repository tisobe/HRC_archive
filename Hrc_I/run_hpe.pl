#!/usr/bin/perl

##############################################################
#
# Run hrc_process_events to creat a new level=1 event list with
# latest CALDB products, software updates, new bad pix file, etc
#
# initialize CIAO before running
#
# JPB, 8 Feb 2007
#
#   t. isobe Mar 09, 2017
#
###############################################################

$outdir = '/data/hrc/i/';

$obsfile=$ARGV[0];

print "$obsfile\n";

$obsfile = "$outdir"."$obsfile";

open(OBS, $obsfile);
while (<OBS>){
    @tmp=split;
    $obsid=trim($tmp[0]);

    $blah=`ls $outdir$obsid/secondary/hrcf*N*evt1.fits*`;
    $file=trim($blah);

    print "\n\n$obsid\t$file\n";

    $dir="$outdir${obsid}/secondary/";

    print `ls $outdir${obsid}/primary/pcad*fits* > $outdir${obsid}/primary/pcad.lst`;

    print `punlearn hrc_process_events`;
    print `pset hrc_process_events infile=${file}`;
#    print `pset hrc_process_events infile=${dir}hrcf${obsid}_evt1_new.fits.gz`;
    print `pset hrc_process_events outfile=${dir}hrcf${obsid}_evt1_new.fits`;
    print `pset hrc_process_events badpixfile=${dir}hrcf${obsid}_new_bpix1.fits`;
    print `pset hrc_process_events acaofffile=\@$outdir${obsid}/primary/pcad.lst`;
    print `pset hrc_process_events instrume=hrc-i`;
    print `pset hrc_process_events do_amp_sf=yes`;
    print `pset hrc_process_events badfile=NONE`;
    print `pset hrc_process_events gainfile=CALDB`;
    print `pset hrc_process_events degapfile=CALDB`;
    print `pset hrc_process_events logfile=${dir}hrc_process_events.log`;
    print `pset hrc_process_events verbose=4`;
    print `pset hrc_process_events mode=h`;
    print `pset hrc_process_events clobber=yes`;
    print `hrc_process_events`;
}

close(OBS);


sub trim {
    my @out = @_;
    for (@out){
        s/^\s+//;
        s/\s+$//;
    }
    return wantarray ? @out: $out[0];
}

