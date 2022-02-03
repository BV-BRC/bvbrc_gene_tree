package Sequence_Alignment;
use strict;
use warnings;
use List::Util qw(max);

our $debug = 0;

sub new {
    my ($class, $input) = @_;
    if ($debug) {
        print(STDERR "in Sequence_Alignment::new\n");
        print(STDERR " class = $class\n");
        print(STDERR " input = $input\n");
    }
    my $self = {};
    bless $self, $class;
    $self->{_seqs} = {};
    $self->{_annot} = {};
    $self->{_ids} = [];
    $self->{_is_alinged} = 0;
    $self->{_length} = 0;
    $self->{_format} = '';
    if ($input) {
        $self->read_file($input)
    }
    print( STDERR "new Sequence_Alignment: ids=", $self->{_ids}, "\n") if $debug;
    return $self;
}

sub set_debug {my $onoff = shift; $debug = $onoff ? 1 : 0}
sub get_ntaxa { my $self = shift; return scalar(@{$self->{_ids}})}
sub get_length { my $self = shift; return $self->{_length}}
sub is_aligned { my $self = shift; return $self->{_is_aligned}}
sub get_ids { my $self = shift; return $self->{_ids}}

sub instantiate_from_hash { 
    my $self = shift;
    my $in_hash = shift;
    $self->{_seqs} = %{$in_hash};
    $self->{_ids} = keys %{$in_hash};
    $self->{_is_aligned} = 1;
    for my $id ($self->{_ids}) {
        my $seqlen = length($self->{_seqs}{$id});
        $self->{_is_aligned} = 0 if $self->{_length} and $seqlen != $self->{_length};
        $self->{_length} = $seqlen if $seqlen > $self->{_length};
    }
    return $self;
}

sub detect_format {
    my $class = shift;
    my $fh = shift;
    print STDERR "in detect_format, class=$class, fh=$fh\n" if $debug;

    $_ = readline $fh;
    seek $fh, 0, 0; # reset file to beginning

    print STDERR "first line of file is :\n", $_, "\n" if $debug;
    my $format = 'unknown';
    $format = 'clustal' if (/^CLUSTAL/ || /^MUSCLE/);
    $format = 'fasta' if (/^>/);
    $format = 'phylip' if (/^(\d+)\s+(\d+)\s*$/);
    $format = 'nexus' if (/\#NEXUS/);
    print STDERR "input format detected as $format\n" if $debug;
    return $format
}

sub read_file {
    my $self = shift;
    my $fh = shift;
    if ( ! ref($fh) ) {
        print STDERR "in read_file, not a file handle, open file $fh\n" if $debug;
        my $temp = undef;
        open $temp, $fh;
        $fh = $temp;
    } 
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
                my $temp = $id;
                my $suffix = 1;
                while (exists $self->{_seqs}{$temp}) {
                    print STDERR "sequence id $temp exists\n" if $debug;
                    $suffix++;
                    $temp = "${id}_$suffix";
                    print STDERR "incrementing to $temp\n" if $debug;
                }
                $id = $temp;
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
    print STDERR "now review sequences, length = $self->{_length}:\n" if $debug;
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
    my $write_unaligned = shift;
    print STDERR "in write_fasta($out), ref(out) = ", ref($out), "\n" if $debug; 
    my $FH = $out;
    unless (ref($out) and ref($out) eq "GLOB") {
        print STDERR "opening $out for fasta output.\n" if $debug;
        open(TEMP, ">", $out);
        $FH = *TEMP;
        die "Cannot open file for writing at $out\n" unless $FH;
    }
    foreach my $id (@{$self->{_ids}}) {
        print $FH ">",$id, "\n";
        die "id from _id is not in _seqs" unless exists $self->{_seqs}->{$id};
        my $seq = $self->{_seqs}->{$id};
        $seq =~ tr/-//d if $write_unaligned;
        die "seq length is zero for $id" unless $seq;
        print $FH "$seq\n";
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
    if ($debug) {
        print STDERR "In write_phylip\n";
        print STDERR "self = $self\n";
        print STDERR "out = $out\n";
    }
    warn "alignment status is $self->{_is_aligned} in write_phylip()" unless $self->{_is_aligned};
    my $FH;
    if (ref($out) and ref($out) eq "GLOB") {
        $out = $FH
    }
    else {
        print STDERR "opening $out for phylip output.\n" if $debug;
        open($FH, ">$out");
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

sub write_fasta_for_raxml {
    my $self = shift;
    my $out = shift;
    print STDERR "in write_fasta_for_raxml, ref(out) = ", ref($out), "\n" if $debug; 
    my $FH = $out;
    unless (ref($out) and ref($out) eq "GLOB") {
        print STDERR "opening $out for fasta output.\n" if $debug;
        open(my $TEMP, ">", $out);
        $FH = $TEMP;
    }
    #my $raxml_illegal_chars = ":()[]";
    my $need_changing = 0;
    for my $id (@{$self->{_ids}}) {
        if ($id =~ tr/:()[]/:()[]/) { # counts but doesn't change
            $need_changing = 1;
            $self->{_raxml_to_orignal_id} = {};
            last
        }
    }
    for my $id (@{$self->{_ids}}) {
        my $seq = $self->{_seqs}->{$id};
        if ($need_changing) {
            my $orig = $id;
            my $changed = $id =~ tr/:()[]/_____/; #replace with underscores
            if ($changed) {
                print STDERR "in write_fasta_for_raxml: original=$orig, changed=$id\n" if $debug;
                $self->{_raxml_to_original_id}{$id} = $orig;
            }
        }
        print $FH ">$id\n$seq\n";
    }
    if ($out ne $FH) {
        print "closing $FH\n";
        close $FH  #because we opened it
    } 
    return $need_changing;
}

sub restore_original_ids_in_raxml_tree {
    my ($self, $newick);
    return $newick unless exists $self->{_raxml_to_orignal_id};
    for my $raxml_id (keys %{$self->{_raxml_to_orignal_id}}) {
        $newick =~ s/$raxml_id/$self->{_raxml_to_orignal_id}{$raxml_id}/;
    }
    return $newick
}


sub calc_column_gap_count {
    my $self = shift;
    #print STDERR "In calc_column_gap_count\n";
    my @gap_count;
    $#gap_count = $self->{_length}-1;
    for my $id (@{$self->{_ids}}) {
        my @str_as_array = split('', $self->{_seqs}->{$id});
        for my $i (0 .. $#str_as_array) {
            $gap_count[$i] += $str_as_array[$i] eq '-';
        }
    }
    return \@gap_count;
}

sub calculate_entropy_per_column {
    my $self = shift;
    my @column_entropy;
    for my $column_index (0 .. $self->get_length()) {
        my %letter_count = ();
        my $num_valid = 0;
        for my $id (@{$self->{_ids}}) {
            my $letter = substr($self->{_seqs}->{$id}, $column_index, 1);
            unless ($letter eq "-") {
                $letter_count{$letter}++;
                $num_valid++
            }
        }
        my $entropy = 0;
        for my $letter (keys %letter_count) {
            my $frequency = $letter_count{$letter} / $num_valid;
            if ($frequency) {
                $entropy += $frequency * log( 1/$frequency )
            }
        }
        push @column_entropy, $entropy;
    }
    return \@column_entropy;
}

sub write_stats {
    my $self = shift;
    my $gaps_per_seq = $self->calc_row_gap_count();
    my $worst_seq = undef;
    my $worst_seq_gaps = 0;
    my $avg_gaps_per_seq = 0;
    my $total_gaps = 0;
    for my $id (keys %$gaps_per_seq) {
        $avg_gaps_per_seq += $gaps_per_seq->{$id};
        $total_gaps += $gaps_per_seq->{$id};
        if ($gaps_per_seq->{$id} > $worst_seq_gaps) {
            $worst_seq_gaps = $gaps_per_seq->{$id};
            $worst_seq = $id
        }
    }
    $avg_gaps_per_seq /= $self->get_ntaxa();

    my $per_col_entropy = $self->calculate_entropy_per_column();
    my $avg_entropy = 0;
    for my $e (@$per_col_entropy) {
        $avg_entropy += $e
    }
    $avg_entropy /= $self->get_length();
    my $retval = "";
    $retval .= "Alignment Statistics\n";
    $retval .= "\tNumber of sequences    = ". $self->get_ntaxa(). "\n";
    $retval .= "\tAlignment length       = ". $self->get_length(). "\n";
    $retval .= "\tProportion gaps        = ". sprintf("%.4f", $avg_gaps_per_seq/$self->get_length()). "\n";
    $retval .= "\tAverage column entropy = ". sprintf("%.3f", $avg_entropy). "\n";
    return $retval;
}

sub end_trim {
    # trim gappy ends inward to a minimum occupancy threshold (proportion of non-gap chars)
    my $self = shift;
    my $threshold = shift;
    print STDERR "In end_trim($threshold)\n" if $debug;
    ($threshold <= 1.0 and $threshold > 0) or die "threshold must be between 0 and 1";
    my $gap_count = $self->calc_column_gap_count();
    my $max_gaps = (1.0 - $threshold) * $self->get_ntaxa();
    my $vis = '';
    my $vis2 = '';
    #print "Length of \@gap_count = ", scalar(@$gap_count), ", vs self->length = $self->{_length}\n";
    for my $i (0 .. $self->{_length}-1) {
        my $prop10 = int(10 * $gap_count->[$i] / $self->get_ntaxa());
        $prop10 = 9 if $prop10 > 9;
        $vis .= $prop10;
        $vis2 .= $i % 10;
    }
    #print "$vis\n$vis2\n";

    my $start = 0;
    $start++ while ($gap_count->[$start] > $max_gaps and $start < $self->{_length}-1);
    my $end = $self->{_length}-1;
    $end-- while ($end and $gap_count->[$end] > $max_gaps);
    my $num_end_columns_trimmed = $self->{_length} - $end - 1;
    my $len = $end - $start + 1;
    print STDERR "Trim up to $start and after $end\n" if $debug;
    #print substr($vis, $start, $len), "\n";
    #print substr($vis2, $start, $len), "\n";
    for my $id (@{$self->{_ids}}) {
            $self->{_seqs}->{$id} = substr($self->{_seqs}->{$id}, $start, $len);
    }
    $self->{_length} = $len;
    return ($start, $num_end_columns_trimmed);
}

sub calc_row_gap_count
{
    my $self = shift;
    print STDERR "In calc_row_gap_count\n" if $debug;
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
    print STDERR "In delete_gappy_seqs($threshold)\n" if $debug;
    ($threshold <= 1.0 and $threshold > 0) or die "threshold must be between 0 and 1";
    my $gap_count = $self->calc_row_gap_count();
    my $max_gaps = (1.0 - $threshold)*$self->{_length};
    my $index = 0;
    my @retval = ();
    while ($index <= $#{$self->{_ids}}) {
        my $id = $self->{_ids}->[$index];
        if ($id ne $self->{_ids}->[$index]) {
            print "item at $index is not $id: ", join(", ", @{$self->{_ids}}[$index - 1, $index, $index + 1]), "\n";
        }
        if ($gap_count->{$id} > $max_gaps) {
            # remove id from list and seq from hash
            print STDERR "seq $id has $gap_count->{$id} gaps, deleting index $index.\n" if $debug;
            my $hashlen_before = scalar keys %{$self->{_seqs}};
            my $arraylen_before = scalar @{$self->{_ids}};
            my @abefore= @{$self->{_ids}}[$index-1,$index, $index+1];
            delete($self->{_seqs}->{$id});
            my $removed = splice @{$self->{_ids}}, $index, 1;
            my $hashlen_after = scalar keys %{$self->{_seqs}};
            my $arraylen_after = scalar @{$self->{_ids}};
            if ($self->{_ids}->[$index] eq $id) {
                print "list before = ",join(", ", @abefore), "\n"; 
                print "list after  = ", join(", ", @{$self->{_ids}}[$index - 1, $index, $index + 1]), "\n";
                print "Removed = $removed\n";
                die "found $id on list after deleting it, index=$index hashlen_b=$hashlen_before, hashlen_a=$hashlen_after; arraylen_b=$arraylen_before, arraylen_a=$arraylen_after"; 
            }
            push @retval, $id;
        }
        else {
            $index++;
        }
    }
    die "array vs hash mismatch after deteting gappy seqs" if (scalar @{$self->{_ids}} != scalar keys %{$self->{_seqs}});
    return \@retval
}    

sub get_sequence_lengths {
    my $self = shift;
    print STDERR "In get_sequence_lengths, num ids: ", scalar @{$self->{_ids}}, "\n" if $debug;
    my %seq_len;
    for my $id (@{$self->{_ids}}) {
        $seq_len{$id} = length($self->{_seqs}->{$id});
    }
    return \%seq_len;
}

return 1
