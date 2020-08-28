#!/usr/bin/perl

##### parsAnnotSV.pl ####

# Author : Thomas Guignard 2020

# Description : 
# Create an User friendly Excel file from an AnnotSV annotated file. 


use strict; 
use warnings;
use Getopt::Long; 
use Switch;
use YAML::XS 'LoadFile';
use Data::Dumper;

use Sort::Key::Natural qw(rnatkeysort);
use Sort::Naturally;

#parameters
my $man = "USAGE : \nperl parsAnnotSV.pl 
\n--configFile <YAML config file for customizing output>
\n--annotSVfile <annotSV annotated file> 
\n--annotSVranking <annotSV ranking explanations file>
\n--outDir <output directory (default = current dir)> 
\n--outPrefix <output file prelifx (default = \"\")>"; 


my $configFile = ".";

my $help;
my $current_line;
my $incfile = "";
my $outDir = ".";
my $outPrefix = "";


#vcf parsing and output
my @line;
my $variantID; 
my $count = 0;

my @finalSortData;
my %hashFinalSortData;
	
my $config;
my %columnHash;
my %InColHash;
my %OutColHash;
my %NameColHash;
my %dataHash;
my %dataCommentHash;
my %dataColorHash;

my $scorePenalty;

#$arguments = GetOptions( "vcf=s" => \$incfile ) or pod2usage(-vcf => "$0: argument required\n") ;

GetOptions( 	"annotSVfile=s"				=> \$incfile,
		"configFile=s"		=> \$config,
		"outDir=s"			=> \$outDir,
		"outPrefix:s"			=> \$outPrefix,
		"help|h"			=> \$help);
				
				


#check mandatory arguments
#if(defined $help || $incfile eq ""){
#	die("$man");
#}

#add underscore to output prefix
if($outPrefix ne ""){
	$outPrefix .= "_";
}


			
			
print  STDERR "Starting a new fishing trip ... \n" ; 



open( CONFIG, "<$config" ) or die ("Cannot open vcf file $config") ;
print STDERR "Parsing config file....\n";
my @configLine;
my %configHash;

my $configHash;

$configHash = LoadFile($config);
	#@configLine = split (/=/ , $_);
	#$configHash{$configLine[0]}=$configLine[1];
	#print Dumper($configHash);


foreach my $item (keys %{$configHash}){
		#print $item."\n";
		my $field = $item;
		#print $configHash->{$field}->{POSITION}."\n";
		
		$OutColHash{$configHash->{$field}->{POSITION}} = $field;
		$NameColHash{$field} = $configHash->{$field}->{POSITION};
		
		
		if (defined $configHash->{$field}->{COMMENTLIST}){
			$dataCommentHash{$field}{'commentFieldList'}=$configHash->{$field}->{COMMENTLIST} ;
		}
		#print Dump($field)."\n";
}







#print Dumper($configHash);

#open( VCF , "<$incfile" )or die("Cannot open vcf file $incfile") ;


#TODO check if header contains required INFO
#Parse VCF header to fill the dictionnary of parameters
print STDERR "Parsing VCF header in order to get sample names and to check if required informations are present ... \n";
my %dicoParam;

#Hash of ACMG incidentalome genes
my %ACMGgene = ("ACTA2" =>1,"ACTC1" =>1,"APC" =>1,"APOB" =>1,"ATP7B" =>1,"BMPR1A" =>1,"BRCA1" =>1,"BRCA2" =>1,"CACNA1S" =>1,"COL3A1" =>1,"DSC2" =>1,"DSG2" =>1,"DSP" =>1,"FBN1" =>1,"GLA" =>1,"KCNH2" =>1,"KCNQ1" =>1,"LDLR" =>1,"LMNA" =>1,"MEN1" =>1,"MLH1" =>1,"MSH2" =>1,"MSH6" =>1,"MUTYH" =>1,"MYBPC3" =>1,"MYH11" =>1,"MYH7" =>1,"MYL2" =>1,"MYL3" =>1,"NF2" =>1,"OTC" =>1,"PCSK9" =>1,"PKP2" =>1,"PMS2" =>1,"PRKAG2" =>1,"PTEN" =>1,"RB1" =>1,"RET" =>1,"RYR1" =>1,"RYR2" =>1,"SCN5A" =>1,"SDHAF2" =>1,"SDHB" =>1,"SDHC" =>1,"SDHD" =>1,"SMAD3" =>1,"SMAD4" =>1,"STK11" =>1,"TGFBR1" =>1,"TGFBR2" =>1,"TMEM43" =>1,"TNNI3" =>1,"TNNT2" =>1,"TP53" =>1,"TPM1" =>1,"TSC1" =>1,"TSC2" =>1,"VHL" =>1,"WT1"=>1);




#Define column title order
my @columnTitles;
#foreach my $key  (sort { $dicoColumnNbr{$a} <=> $dicoColumnNbr{$b} } keys %dicoColumnNbr)  {
#	push @columnTitles,  $key;
#	#DEBUG print STDERR $key."\n";
#}


#final strings for comment
my $commentMPAscore;

#define sorted arrays with score for comment
my @CommentMPA_score = ("MPA_ranking",
						"MPA_final_score",
						'fathmm-MKL_coding_pred');

####################################################################
#############################################
##################   Start parsing VCF

open( VCF , "<$incfile" )or die("Cannot open vcf file ".$incfile) ;



while( <VCF> ){
  	$current_line = $_;
	$count++;

#############################################
##############   skip header

	chomp $current_line;
	@line = split( /\t/, $current_line );	
	
	#DEBUG print STDERR $dicoColumnNbr{'Gene.refGene'}."\n";



#############################################
##############   Treatment for First line to create header of the output

	if ( $line[0] eq "AnnotSV ID" )   {

		#initialize column order
		for( my $fieldNbr = 0 ; $fieldNbr < scalar @line; $fieldNbr++){
			
			#print STDERR "Field found:  ".$line[$fieldNbr]."\n";
			$InColHash{$fieldNbr} = $line[$fieldNbr];
			
			$dataHash{$line[$fieldNbr]} = ".";
			

			#$columnHash{$fields}{"newColPosition"}=$configHash{$line[$fields]};
			#$columnHash{$fields}{"isComment"}=0;
			#$columnHash{$fields}{"hasComment"}=0;
			#$columnHash{$fields}{"commentFields"}="";
			#$columnHash{$fields}{"finalColName"}=$configHash{$line[$fields]};

			
		}

		next;
		
#############################################
##############################
##########  start to compute variant lines	

	}else {

		#initialise final printable string
		@finalSortData = ("");

		#fill printable string
		for( my $fieldNbr = 0 ; $fieldNbr < scalar @line; $fieldNbr++){

			#if (defined $columnHash{$fieldNbr}{"newColPosition"} and $columnHash{$fieldNbr}{"newColPosition"} > 0) {
			#	$finalSortData[$columnHash{$fieldNbr}{"newColPosition"}]=$line[$fieldNbr]-1;
			#}

			$dataHash{$InColHash{$fieldNbr}} = $line[$fieldNbr];

		}

		#get all comments values
		foreach my $field (keys %dataCommentHash){
			if (defined $dataCommentHash{$field}){
				foreach my $fieldCom (@{$dataCommentHash{$field}{'commentFieldList'}}){
					$dataCommentHash{$field}{'values'} .= $fieldCom . " : ".$dataHash{$fieldCom}."\n\n";
				}
			
			}
		}

		#fill finalSortData array
		foreach my $field (keys %dataHash){
			if (defined $NameColHash{$field}){
				$finalSortData[$NameColHash{$field} - 1] = $dataHash{$field};
			}
		}



		
		#	print Dumper(\@finalSortData);

#############################################################################
######### FILL HASH STRUCTURE WITH COMMENTS 
		
		    
            $variantID = $line[0];

			#key for good sorting should be rank concatenated with annotSV_ID
			#
			#
			if ($dataHash{"AnnotSV type"} eq "full"){
				$scorePenalty = ".1";	
			}else{
				$scorePenalty = ".0";	

			}	
			$hashFinalSortData{$dataHash{"AnnotSV ranking"}.$scorePenalty."_".$variantID}{$count}{'finalArray'} = [@finalSortData] ; 
			
			foreach my $fieldCom (keys %dataCommentHash){
				
				$hashFinalSortData{$dataHash{"AnnotSV ranking"}."_".$variantID}{'hashComments'}{$NameColHash{$fieldCom}} = $dataCommentHash{$fieldCom} ; 
				#TODO check color for gene according to output position of gene field 
				$hashFinalSortData{$dataHash{"AnnotSV ranking"}."_".$variantID}{'hashColor'}{$NameColHash{$fieldCom}} = "#FF0000" ; 
			}


			#$hashFinalSortData{$finalSortData[$configHash{'AnnotSV ranking'}]}{$variantID}{'finalArray'} = [@finalSortData] ; 
			#$hashFinalSortData{$finalSortData[$configHash{'AnnotSV ranking'}]}{$variantID}{'hashComments'} = [@finalSortData] ; 
#			$hashFinalSortData{$finalSortData[$dicoColumnNbr{'AnnotSV Ranking'}]}{$variantID}{'finalArray'} = [@finalSortData] ;


#			$hashFinalSortData{$finalSortData[$dicoColumnNbr{'MPA_ranking'}]}{$variantID}{'nbSample'} = scalar keys %dicoSamples ;
#			$hashFinalSortData{$finalSortData[$dicoColumnNbr{'MPA_ranking'}]}{$variantID}{'commentGnomADexome'} = $commentGnomADExomeScore  ;

				#ACMG
				#if(defined $ACMGgene{$finalSortData[$dicoColumnNbr{'Gene.refGene'}]} )
#				if(defined $ACMGgene{$geneName} ){
#					$hashFinalSortData{$finalSortData[$dicoColumnNbr{'MPA_ranking'}]}{$variantID}{'worksheet'} .= "#ACMG";
#				}


	} #END of IF-ELSE(#CHROM)	

}#END WHILE VCF


############ THE TIME HAS COME TO FILL THE HTML OUTPUT FILE
####################################################################
######################### HTML  Initialisation ######################
#
#file header
my $htmlStart = "<!DOCTYPE html>\n<html>
\n<head>
\n<meta charset=\"utf-8\">
\n<title>".$outPrefix." AnnotSV</title>\n
\n<script type=\"text/javascript\" src='DataTables/js/jquery.js'></script>
\n<script type=\"text/javascript\" src='DataTables/js/jquery.dataTables.min.js'></script>
\n<script type=\"text/javascript\" src='achab.js'></script>
\n<link rel=\"stylesheet\" type=\"text/css\" href='DataTables/css/jquery.dataTables.min.css'>
\n<link rel=\"stylesheet\" type=\"text/css\" href='DataTables/css/W3CSS.css'>
\n</head>
\n\t<body>";

#table and columns names



my $htmlStartTable = "\n\t<table id='tab' class='display' >
        \n\t\t<thead><tr>";
		foreach my $col (sort {$a <=> $b} keys %OutColHash){
			        #print HTML "\t<th style=\"word-wrap: break-word\"   >";
					        $htmlStartTable .= "\t<th >".$OutColHash{$col}."\t</th>\n";
        }
$htmlStartTable .= "</tr>\n</thead>\n<tbody>\n";
					
					
					
my $htmlEndTable = "</tbody>\n</table>\n</div>\n";
					
					
					


my $htmlEnd = "<div class=\"tab\">\n
<button class=\"tablinks\" onclick=\"openCity(event, 'METADATA')\">METADATA</button>\n
<button class=\"tablinks\" onclick=\"openCity(event, 'ALL')\">ALL</button>\n
</div>";


$htmlEnd .= "\n</body>\n</html>";
my $htmlALL= "<div id=\"ALL\" class=\"tabcontent\">". $htmlStartTable;
my $htmlACMG="<div id=\"ACMG\" class=\"tabcontent\">". $htmlStartTable;






open(HTML, '>', $outDir."/".$outPrefix."annotSV.html") or die $!;


#########################################################################
#################### Sort by MPA ranking for the output

#create user friendly ranking score
my $kindRank=0;

foreach my $rank (rnatkeysort { "$_-$hashFinalSortData{$_}" } keys %hashFinalSortData){
	print $rank."\n";
	
	foreach my $variant ( keys %{$hashFinalSortData{$rank}}){

#		$format_pLI = $workbook->add_format(bg_color => '#FFFFFF');
		
		#increse rank number then change final array
		#$kindRank++;
		
		#print $variant ."___".  $hashFinalSortData{$rank}{$variant}{'finalArray'}[$dicoColumnNbr{'MPA_ranking'}]."\n";
		#$hashFinalSortData{$rank}{$variant}{'finalArray'}[$dicoColumnNbr{'MPA_ranking'}] = $kindRank; 

		#last finalSortData assignation
		#$hashFinalSortData{$finalSortData[$dicoColumnNbr{'MPA_ranking'}]}{$variantID}{'finalArray'} = [@finalSortData] ; 
		



		#create reference of Hashes
		#my $hashTemp = $hashFinalSortData{$rank}{$variant};
		#my $hashColumn_ref = \%dicoColumnNbr;




#FILL tab 'ALL';
$htmlALL .= "<tr>\n";

#foreach my $field ( @{$hashFinalSortData{$rank}{$variant}{'finalArray'}} ){
foreach my $field ( @{$hashFinalSortData{$rank}{$variant}{'finalArray'}} ){
		print $rank."bis\n";
        #print HTML "\t<td style=\"word-wrap: break-word\";class=\"tooltip\"; title=".$field."  >";
        #print HTML "\t<td >";
	$htmlALL .= "\t<td >";
	if (defined $field){
	#print HTML $field;
	$htmlALL .= $field;
	}else {
		#print HTML "."
        $htmlALL .= "."
	}
    #print HTML "\t</td>\n";
	$htmlALL .= "\t</td>\n";
}
	
	#print HTML "</tr>\n";
	$htmlALL .= "</tr>\n";
	
	
	




##############################################################
###########################      ALL     #####################

	
	}   # END FOREACH VARIANT
}	#END FOREACH RANK



close(VCF);




print HTML $htmlStart;
print HTML $htmlALL;
print HTML $htmlEndTable;
print HTML $htmlEnd;


close(HTMl);




print STDERR "Done!\n\n\n";

###################################################################
######################### FUNCTIONS ###############################
# Write in sheets
sub writeThisSheet {
}#END OF SUB

exit 0; 
