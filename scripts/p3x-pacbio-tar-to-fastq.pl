#!/usr/bin/env perl
use strict;

if (!$ARGV[0] or $ARGV[0] eq '-h') {
    print 'Usage: p3x-pacbio-tar-to-fastq.pl /user@bv-brc.org/home/path-to-pacbio.tar';
}
my $filepath;
if ($ARGV[0] =~ /^\d+$/) {
    # looks like a job id, get the read file from the task-status stdout file
    # this requires being on holly
    print STDERR "Try looking in /disks/p3/task_status/$ARGV[0]/stdout\n";
    open F, "/disks/p3/task_status/$ARGV[0]/stdout" or die "cannot open job status stdout for $ARGV[0]";
   
    while (<F>) {
        $filepath = $1 if /'read' => '(\/.*brc.*tar)'/;
        #'read' => '/gourabelsahalder@patricbrc.org/home/apcolr130/barcode19apcolr130.tar'
    }
    die "Could not find read file ending in 'tar' for job $ARGV[0]" unless $filepath;
}
else {
    $filepath = $ARGV[0]
}
print "filepath = $filepath\n";

my @fields = split("/", $filepath);
my $filename = pop @fields;
my $filename_fastq = $filename;
$filename_fastq =~ s/.tar$/_fastq.gz/;
if (-f $filename) {
    print "Refusing to overwrite $filename, exiting)";
    exit(1);
}
if (-f $filename_fastq) {
    print "Refusing to overwrite $filename_fastq, exiting)";
    exit(1);
}
my $workspace_directory = join("/", @fields);
shift @fields; # get rid of first empty field
my $user_brc = shift @fields;
my ($user, $brc) = split("@", $user_brc);
die "Cannot parse the username in $filepath" unless $user;

my $ls_result = `p3-ls "$filepath"`;
chomp($ls_result);
if ($ls_result ne $filepath) {
    print "Cannot see file in workspace: $ls_result\n";
    exit(1);
}
system("p3-cp ws:'$filepath' .");
if (! -f $filename) {
    print "$filename could not be copied.\n";
    exit(1);
}
system("tar -xOf '$filename' > '$filename_fastq'");
system("p3-cp -m gz=reads '$filename_fastq' ws:'$workspace_directory'");
unlink($filename);
unlink($filename_fastq);
system("p3-ls -lT '$workspace_directory/$filename_fastq'");
