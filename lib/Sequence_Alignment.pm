package Sequence_Alignment;
use strict;
use warnings;
use List::Util qw(max);

sub new {
    my ($class, $fh) = @_;
    print(STDERR "in Sequence_Alignment::new\n");
    print(STDERR " class = $class\n");
    print(STDERR " args = ", ", ".join(@_), ".\n");
    my $self = {};
    bless $self, $class;
    $self->{_seqs} = {};
    $self->{_annot} = {};
    $self->{_ids} = [];
    $self->{_is_alinged} = 0;
    $self->{_length} = 0;
    $self->{_format} = '';
    if (defined $fh) {
        print STDERR "File handle = $fh\n";
        $self->read_file($fh);
    }
    return $self;
}

sub get_ntaxa { my $self = shift; return scalar(@{$self->{_ids}})}
sub get_length { my $self = shift; return $self->{_length}}
sub is_aligned { my $self = shift; return $self->{_is_aligned}}

sub detect_format {
    my $class = shift;
    my $fh = shift;
    print STDERR "in detect_format, class=$class, fh=$fh\n";

    $_ = readline $fh;
    seek $fh, 0, 0; # reset file to beginning

    print STDERR "first line of file is :\n", $_, "\n";
    my $format = 'unknown';
    $format = 'clustal' if (/^CLUSTAL/ || /^MUSCLE/);
    $format = 'fasta' if (/^>/);
    $format = 'phylip' if (/^(\d+)\s+(\d+)\s*$/);
    $format = 'nexus' if (/\#NEXUS/);
    print STDERR "input format detected as $format\n";
    return $format
}

sub read_file {
    my $self = shift;
    my $fh = shift;
    print STDERR "in detect_format, self=$self, fh=$fh\n";
    my $format = shift;
    if (! defined $format)
    {
        $format = $self->detect_format($fh)
    }
    $self->{_format} = $format;
    if ($format eq 'unknown') {
        return undef}

    if ($format eq 'clustal') {
        my $found = 0;
        while (<$fh>) {
            if (/^CLUSTAL/ || /^MUSCLE/) {
                $found = 1;
                last;
            }
        }
        die "Format seems to be wrong, not Clustal.\n" if (!$found); 
        while (<$fh>) {
            if (/^(\S+)\s+(\S+)/) {
                push(@{$self->{_ids}}, $1) unless exists($self->{_seqs}{$1});

                $self->{_seqs}{$1} .= $2;
            }
        }
    }
    elsif ($format eq 'phylip') {
        $_ = <$fh>;
        die "Format does not seem to be phylip\n" if (!/^\s*(\d+)\s+(\d+)\s*$/);
        my $ntaxa = $1;
        my $nchar = $2;
        for my $i (1..$ntaxa) {
            $_ = <$fh>;
            /(\S+)\s+(\S.*\S)/ or die $_;
            my $id = $1;
            my $seq = $2;
            $seq =~ s/\s//g;
            $self->{_seqs}{$id} = $seq;
            push @{$self->{_ids}}, $id;
        }
        # now if there are more lines, read in same order as first set, but without identifiers
        my $index = 0;
        while (<$fh>) {
            my $seq = $_;
            $seq =~ s/\s//g;
            if ($seq) {
                my $id = $self->{_ids}[$index % $ntaxa];
                $self->{_seqs}{$id} .= $seq;
                $index += 1
            }
        }
        foreach my $id (@{$self->{_ids}}) { # phylip uses '.' as insert (unknown) character
            $self->{_seqs}{$id} =~ s/\./\-/g; # replace dot as gap char with '-'
        }
    }
    elsif ($format eq 'fasta') {
        my $id;
        while (<$fh>) {
            chomp;
            if (/^>(\S+)/) {
                $id = $1;
                push @{$self->{_ids}}, $id;
                #$self->{_annot}{$id} = $2 if $2;
                $self->{_seqs}{$id} = '';
            }
            else  {
                $_ =~ tr/\\s//d;
                $self->{_seqs}{$id} .= $_;
            }
        }
    }
    elsif ($format eq 'nexus')
    {
        $_ = <$fh>;
        die "Format does not seem to be NEXUS" if (!/\#NEXUS/);
        my $data;
        my $ntax;
        my $nchar;
        my $matrix;
        while (<$fh>) {
            chomp;
            s/\[[^\]]*\]//g;
            $data = 1 if (/^begin data/i);
            if ($data and !$matrix and /^dimensions/i) {
                $ntax = $1 if (/ntax=(\d+)/i);
                $nchar = $1 if (/nchar=(\d+)/i);
            }
            $matrix = 1 if ($data and /^matrix/i);
            if ($matrix) {	    
                if (/^(\S+)\s+(\S+)/) {
                    push @{$self->{_ids}}, $1 unless $self->{_seqs}{$1};
                    $self->{_seqs}{$1} .= $2;
                }
                last if (/;/);
            }
        }
    }
    my $first_id = $self->{_ids}[0];
    $self->{_length} = length($self->{_seqs}{$first_id});
    $self->{_is_aligned} = 1;
    print STDERR "now review sequences, length = $self->{_length}:\n";
    for my $id (@{$self->{_ids}}) {
        if (length($self->{_seqs}{$id}) != $self->{_length}) {
            $self->{_is_aligned} = 0;
            $self->{_length} = max($self->{_length}, length($self->{_seqs}{$id}))
        }
        #print STDERR "id $id ; len ", length($self->{_seqs}{$id}), " ; is_al=$self->{_is_aligned} \n";
    }
}

sub write_fasta {
  # write out in fasta format
    my $self = shift;
    my $out = shift;
    print STDERR "in write_fasta, ref(out) = ", ref($out), "\n"; 
    my $FH = $out;
    unless (ref($out) and ref($out) eq "GLOB") {
        print STDERR "opening $out for fasta output.\n";
        open(my $TEMP, ">", $out);
        $FH = $TEMP;
    }
    foreach my $id (@{$self->{_ids}}) {
        print $FH ">",$id, "\n", $self->{_seqs}->{$id}, "\n";
    }
    if ($out ne $FH) {
        print "closing $FH\n";
        close $FH  #because we opened it
    } 
}

sub write_phylip {
  # write out in phylip format
    my $self = shift;
    my $out = shift;
    print STDERR "In write_phylip\n";
    print STDERR "self = $self\n";
    print STDERR "out = $out\n";
    warn "alignment status is $self->{_is_aligned} in write_phylip()" unless $self->{_is_aligned};
    my $FH = $out;
    unless (ref($out) and ref($out) eq "GLOB") {
        print STDERR "opening $out for fasta output.\n";
        open(my $TEMP, ">", $out);
        $FH = $TEMP;
    }
    print $FH $self->get_ntaxa(), "  ", $self->get_length(), "\n";
    my $maxIdLength = 0;
    foreach my $id (@{$self->{_ids}}) {
        $maxIdLength = length($id) if length($id) > $maxIdLength;
    }
    foreach my $id (@{$self->{_ids}}) {
        my $seq = $self->{_seqs}->{$id};
        next unless $seq;
        #$seq = ambiguateEndGaps($seq) if ($opt_e);
	    printf($FH "%-${maxIdLength}s %s\n", $id, $seq);
    }
    if ($out ne $FH) {
        print "closing $FH\n";
        close $FH  #because we opened it
    } 
}

sub calc_column_gap_count {
    my $self = shift;
    print STDERR "In calc_column_gap_count\n";
    my @gap_count;
    $#gap_count = $self->{_length}-1;
    for my $id (@{$self->{_ids}}) {
        my @str_as_array = split('', $self->{_seqs}->{$id});
        for my $i (0 .. $#str_as_array) {
            $gap_count[$i] += $str_as_array[$i] eq '-';
            #$gap_count[$i] += substr($self->{_seqs}->{$id},$i, 1) eq '-';
        }
    }
    return \@gap_count;
}


sub end_trim {
    # trim gappy ends inward to a minimum occupancy threshold (proportion of non-gap chars)
    my $self = shift;
    my $threshold = shift;
    print STDERR "In end_trim($threshold)\n";
    ($threshold <= 1.0 and $threshold > 0) or die "threshold must be between 0 and 1";
    my $gap_count = $self->calc_column_gap_count();
    my $max_gaps = (1.0 - $threshold) * $self->get_ntaxa();
    my $start = 0;
    my $vis = '';
    my $vis2 = '';
    #print "Length of \@gap_count = ", scalar(@$gap_count), ", vs self->length = $self->{_length}\n";
    for my $i (0 .. $self->{_length}-1) {
        my $prop10 = int(10 * $gap_count->[$i] / $self->get_ntaxa());
        $prop10 = 9 if $prop10 > 9;
        $vis .= $prop10;
        $vis2 .= $i % 10;
    }

    print "$vis\n$vis2\n";
    $start++ while ($start < $self->{_length}-1 and $gap_count->[$start] > $max_gaps);
    my $end = $self->{_length}-1;
    $end-- while ($end and $gap_count->[$end] > $max_gaps);
    my $len = $end - $start + 1;
    print "Trim up to $start and after $end\n";
    print substr($vis, $start, $len), "\n";
    print substr($vis2, $start, $len), "\n";
    for my $id (@{$self->{_ids}}) {
            $self->{_seqs}->{$id} = substr($self->{_seqs}->{$id}, $start, $len);
    }
    $self->{_length} = $len;
}

sub calc_row_gap_count
{
    my $self = shift;
    print STDERR "In calc_row_gap_count\n";
    my %gap_count;
    for my $id (@{$self->{_ids}}) {
        $gap_count{$id} = $self->{_seqs}->{$id} =~ tr/-/-/;
    }
    return \%gap_count;
}

sub delete_gappy_seqs {
    # remove gappy sequences below minimum occupancy threshold (proportion of non-gap chars)
    my $self = shift;
    my $threshold = shift;
    print STDERR "In delete_gappy_seqs($threshold)\n";
    ($threshold <= 1.0 and $threshold > 0) or die "threshold must be between 0 and 1";
    my $gap_count = $self->calc_row_gap_count();
    my $max_gaps = (1.0 - $threshold)*$self->{_length};
    my $index = 0;
    for my $id (@{$self->{_ids}}) {
        if ($gap_count->{$id} > $max_gaps) {
            # remove id from list and seq from hash
            delete($self->{_seqs}->{$id});
            splice @{$self->{_ids}}, $index, 1;
            print STDERR "seq $id has $gap_count->{$id} gaps, deleting.\n";
        }
        else {
            $index++;
        }
    }
}    

return 1
