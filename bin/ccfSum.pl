#! usr/bin/perl
use warnings; use strict; use Cwd 'abs_path';

die "Usage: $0 fastA_file out_directory locus_tag cpu\n" unless @ARGV == 4;
my $thispath = abs_path($0); $thispath =~ s/[^\/]+$//;
my ($fasta, $outdir, $name, $cpu) = @ARGV; $fasta = abs_path($fasta);
my $ccfPath = `which CRISPRCasFinder.pl 2> /dev/null` or die "CRISPRCasFinder.pl is not in your path\n";
chomp $ccfPath; $ccfPath =~ s/CRISPRCasFinder.pl$/bin/;
my (%gff, $parent, %serials, %parents, $dna, %evidences);

mkdir $outdir; chdir $outdir;
unless (-d 'ccf' and -f "ccf/result.json") {
 system "rm -rf ccf";
 my $call = "CRISPRCasFinder.pl -cf CasFinder-2.0.2 -def General -cas -i $fasta -minSeqSize 500 -cpuM $cpu -out ccf -keep -so $ccfPath/sel392v2.so";
 print "Running from $outdir: $call\n";
 system $call;
}

open OUT, ">ccf.gff";
chdir 'ccf';
for (`cat TSV/Crisprs_REPORT.tsv`) {
 chomp; my @f = split "\t";
 next unless $f[1] and $f[5] and $f[6] and $f[5] !~ /[^0-9]/ and $f[6] !~ /[^0-9]/;
 my ($dna, $L, $R, $evidence) = ($f[1], $f[5], $f[6], $f[-1]);
 $evidences{"$dna,$L,$R"} = $evidence;
}
for my $file (glob "GFF/*gff") {
 next if $file =~ /^GFF\/annotation_/;  # Prodigal outputs
 for (`cat $file`) {
  # CP014688.1      CRISPRCasFinder CRISPR  313945  314034  .       -       .       DR=CTCGGCTCATCCCCGCACACGCGGGGAACAC;DR_length=31;Number_of_spacers=1;Name=CP014688.1_313945_314034;ID=CP014688_Crispr_2;potential_direction=Reverse;
  # CP014688.1      CRISPRCasFinder CRISPRdr        313945  313975  .       -       .       sequence=CTCGGCTCATCCCCGCACACGCGGGGAACAC;Parent=CP014688.1_313945_314034;ID=DR_313945
  # CP014688.1      CRISPRCasFinder CRISPRspacer    313976  314003  .       -       .       sequence=CCGAAGCGCAATGTTCTCCGCGTTCTCT;Name=spacer_313976_28;Parent=CP014688.1_313945_314034;ID=sp_313976
  chomp;
  my @f = split "\t";
  next unless $f[2];
  if ($f[2] eq 'CRISPR' and /Number_of_spacers=([^;]+);.*Name=([^;]+);.*ID=([^;]+)/) {%{$parents{$2}} = (dna => $f[0], L => $f[3], R => $f[4], strand => $f[6],
   spacer_count => $1, id => "$name.$3")}
  elsif ($f[2] eq 'CRISPRdr' and /sequence=([^;]+);.*Parent=([^;]+)/) {push @{$parents{$2}{parts}}, uc($1)}
  elsif ($f[2] eq 'CRISPRspacer' and /sequence=([^;]+);.*Parent=([^;]+)/) {push @{$parents{$2}{parts}}, lc($1)}
  elsif ($f[2] =~ /FLANK/) {}
  else {die "Can't parse $_\n"}
 }
}
for (keys %parents) {
 $parents{$_}{seq} = join '', @{$parents{$_}{parts}};
 $parents{$_}{seq} = Revcomp($parents{$_}{seq}) if $parents{$_}{strand} eq '-';
 my ($dna, $L, $R) = ($parents{$_}{dna}, $parents{$_}{L}, $parents{$_}{R});
 %{$gff{$dna}{$L}{$R}{CRISPR}} = (strand => $parents{$_}{strand}, ID => $parents{$_}{id}, spacers => $parents{$_}{spacer_count},
  sequence => $parents{$_}{seq}, cat => 'CRISPR', evidence => $evidences{"$dna,$L,$R"});
}

%parents = ();
open SUMS, ">sums.bed";
for (`cat TSV/Cas_REPORT.tsv`) {
 # CP014688.1_298  cse1_TypeIE     mandatory       CAS-TypeIE      CDS     278410  279909  +       ID=1_298;partial=00;start_type=TTG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.544;conf=100.00;score=98.68;cscore=98.19;sscore=0.49;rscore=6.98;uscore=0.69;tscore=-7.17;
 # ####Summary system CAS:begin=1034;end=9419:{sequenceID=ATAH01000101} : [cas2_TypeIE (9120,9419,+); cas3_TypeI (1034,3706,+); cas3_TypeI (1034,3706,+); cas5_TypeIE (6958,7620,+); cas6_TypeIE (7595,8206,+); cas7_TypeIE (5814,6956,+); cse1_TypeIE (3709,5277,+); cse2_TypeIE (5274,5804,+)]
 if (/^([^#]\S*)/) {$dna = $1; $dna =~ s/_\d+$//}
 next unless /^####Summary system ([^:]+):begin=(\d+);end=(\d+).*\[([^\]]+)/;
 my $parent = $1; $serials{item} ++; $parent .= "_$serials{item}";
 %{$parents{$parent}} = (dna => $dna, type => $1, L => $2, R => $3, cds => $4);
 print SUMS join("\t", $dna, $2-1, $3, $parent), "\n";
}
close SUMS;
system "perl ${thispath}overlapBed.pl sums.bed > sums.overlap";  # Overlapping systems on same line

for (`cat sums.overlap`) {
 chomp; my @f = split "\t";
 my ($parent, $type, %cds) = ($f[0], $parents{$f[0]}{type}); $type .= 'plus' if @f > 1; 
 my %strands = ('-' => 0, '+' => 0);
 my ($dna, $L, $R, $strand) = ($parents{$parent}{dna}, $parents{$parent}{L}, $parents{$parent}{R}, '+');
 for my $p2 (@f) {
  $L = $parents{$p2}{L} if $parents{$p2}{L} < $L;
  $R = $parents{$p2}{R} if $parents{$p2}{R} > $R;
  for (split '; ', $parents{$p2}{cds}) {
   warn "Parse $name $parent $_\n" unless /(\S+) \((\d+),(\d+),([\+\-])/;
   %{$cds{$2}} = (name => $1, L => $2, R => $3, strand => $4);
  }
 }
 for (keys %cds) {
  $strands{$cds{$_}{strand}} ++;
  %{$gff{$dna}{$cds{$_}{L}}{$cds{$_}{R}}{CDS}} = (strand => $cds{$_}{strand}, Name => $cds{$_}{name}, Parent => $parent, cat => $cds{$_}{name});
 }
 $strand = '-' if $strands{'-'} > $strands{'+'};
 %{$gff{$dna}{$L}{$R}{Cas_system}} = (strand => $strand, cds => scalar(keys %cds), type => $type, cat => $type, parentID => $parent);
}
%serials = ();
for my $dna (sort keys %gff) {
 for my $L (sort {$a <=> $b} keys %{$gff{$dna}}) {
  for my $R (sort {$b <=> $a} keys %{$gff{$dna}{$L}}) {
   for my $feature (qw/Cas_system CDS CRISPR/) {
    next unless $gff{$dna}{$L}{$R}{$feature};
    my $g = $gff{$dna}{$L}{$R}{$feature};
    $serials{$$g{cat}} ++; $$g{ID} = "${name}_$$g{cat}_$serials{$$g{cat}}";
    $parents{$$g{parentID}}{id} = $$g{ID} if $feature eq 'Cas_system';
    $$g{Parent} = $parents{$$g{Parent}}{id} if $feature eq 'CDS';
    print OUT join("\t", $dna, 'CRISPRCasFinder', $feature, $L, $R, '.', $$g{strand}, '.', '');
    for (qw/ID Name Parent type cds spacers sequence evidence/) {
     print OUT "$_=$gff{$dna}{$L}{$R}{$feature}{$_};" if $gff{$dna}{$L}{$R}{$feature}{$_};
    }
    print OUT "\n";
   }
  }
 }
}
close OUT;

sub Revcomp {my $ret = reverse $_[0]; $ret =~ tr/ACGTacgt/TGCAtgca/; return $ret}
