#!/usr/bin/perl

my $repnum;

if (@ARGV == 0){
    $repnum = 00
}else{
    $repnum = $ARGV[0]
}

$ENV{'PATH'} = "../bin:$ENV{'PATH'}";
# MUST put smartpca bin directory in path for smartpca.perl to work

$command = "smartpca.perl";
$command .= " -i rep". $repnum .".geno ";
$command .= " -a rep". $repnum .".snp ";
$command .= " -b mssamples.ind " ;
$command .= " -k 2 ";
$command .= " -o rep". $repnum .".pca ";
$command .= " -p rep". $repnum .".plot ";
$command .= " -e rep". $repnum .".eval ";
$command .= " -l rep". $repnum .".log ";
$command .= " -m 5 ";
$command .= " -t 2 ";
$command .= " -s 6.0 ";
print("$command\n");
system("$command");

#$command = "smarteigenstrat.perl ";
#$command .= " -i rep". $repnum .".geno ";
#$command .= " -a rep". $repnum .".snp ";
#$command .= " -b mssamples.ind ";
#$command .= " -p rep". $repnum .".pca ";
#$command .= " -k 1 ";
#$command .= " -o rep". $repnum .".chisq ";
#$command .= " -l rep". $repnum .".2.log ";
#print("$command\n");
#system("$command");
#
#$command = "gc.perl rep". $repnum .".chisq rep". $repnum .".chisq.GC";
#print("$command\n");
#system("$command");
