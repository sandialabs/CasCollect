use strict; use warnings;

die "perl $0 4-columnBedFile\n" unless @ARGV == 1;
my $file = $ARGV[0];
my (%mems, %outs, %mem2gp);
for (`bedtools intersect -wo -a $file -b $file`) {
 my @f = split "\t";
 $mems{$f[3]}{$f[7]} ++;
 $mems{$f[7]}{$f[3]} ++;
}
for my $mem (keys %mems) {
 my (%gps);
 for (keys %{$mems{$mem}}) {$gps{$mem2gp{$_}} ++ if $mem2gp{$_}}
 my $gp = $mem;
 $gp = (keys %gps)[0] if %gps;
 for (keys %{$mems{$mem}}) {$mem2gp{$_} = $gp}  # Performs merging if needed
}

for my $mem (sort keys %mem2gp) {
 $outs{$mem2gp{$mem}}{$mem} ++;
}

for my $gp (sort {scalar(keys %{$outs{$b}}) <=> scalar(keys %{$outs{$a}}) || $a cmp $b} keys %outs) {
 print join("\t", sort {length($b) <=> length($a)} keys %{$outs{$gp}}), "\n";
}
