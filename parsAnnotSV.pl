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
#use Sort::Naturally;

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
my %SV_ID;


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


#initialize gene colore according to pLI/LOEUF
my %pLI_ColorHash;
my $pLI_Color;

$pLI_ColorHash{'0.9'} = '#FF0000';
$pLI_ColorHash{'0.8'} = '#FF3300';
$pLI_ColorHash{'0.7'} = '#FF6600';
$pLI_ColorHash{'0.6'} = '#FF9900';
$pLI_ColorHash{'0.5'} = '#FFCC00';
$pLI_ColorHash{'0.4'} = '#FFFF00';
$pLI_ColorHash{'0.3'} = '#BFFF00';
$pLI_ColorHash{'0.2'} = '#7FFF00';
$pLI_ColorHash{'0.1'} = '#3FFF00';
$pLI_ColorHash{'0.0'} = '#00FF00';
#$pLI_ColorHash{'.'} = '#FFFFFF';

#$pLI_ColorHash{'0.0'} = '#FF0000';
#$pLI_ColorHash{'1.0'} = '#FF3300';
#$pLI_ColorHash{'2.0'} = '#FF6600';
#$pLI_ColorHash{'3.0'} = '#FF9900';
#$pLI_ColorHash{'4.0'} = '#FFCC00';
#$pLI_ColorHash{'5.0'} = '#FFFF00';
#$pLI_ColorHash{'6.0'} = '#BFFF00';
#$pLI_ColorHash{'7.0'} = '#7FFF00';
#$pLI_ColorHash{'8.0'} = '#3FFF00';
#$pLI_ColorHash{'9.0'} = '#00FF00';
#$pLI_ColorHash{'.'} = '#FFFFFF';


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
        $pLI_Color=".";






        #fill nbr of SV_ID;
        if (defined $SV_ID{$line[0]} ){
            $SV_ID{$line[0]}++;
        }else{
            $SV_ID{$line[0]} = 1;
            #reinitialize Comment Values
		    foreach my $field (keys %dataCommentHash){
                undef($dataCommentHash{$field}{'values'});
            }

        }



		#fill printable string
		for( my $fieldNbr = 0 ; $fieldNbr < scalar @line; $fieldNbr++){

			#if (defined $columnHash{$fieldNbr}{"newColPosition"} and $columnHash{$fieldNbr}{"newColPosition"} > 0) {
			#	$finalSortData[$columnHash{$fieldNbr}{"newColPosition"}]=$line[$fieldNbr]-1;
			#}

			$dataHash{$InColHash{$fieldNbr}} = $line[$fieldNbr];

		}

		#get all comments values
		foreach my $field (keys %dataCommentHash){
            #print $field."\n";
			#if (defined $dataCommentHash{$field}){
			#	foreach my $fieldCom (@{$dataCommentHash{$field}{'commentFieldList'}}){
			if (defined $dataCommentHash{$field}){
                if (! defined $dataCommentHash{$field}{'values'}){
				    foreach my $fieldCom (@{$dataCommentHash{$field}{'commentFieldList'}}){
					    $dataCommentHash{$field}{'values'} .= "<br><br><b>".$fieldCom . " :</b> ".$dataHash{$fieldCom};
                        #print $dataCommentHash{$field}{'values'}."\n";
				    }
			    }
			}
		}

		#fill finalSortData array
		foreach my $field (keys %dataHash){
			if (defined $NameColHash{$field}){
				$finalSortData[$NameColHash{$field} - 1] = $dataHash{$field};
			}
		}


		#fill color pLI 
		foreach my $field (keys %dataHash){
			if (defined $dataHash{'pLI_ExAC'} ){
                foreach my $pli (sort {$b <=> $a} keys %pLI_ColorHash){
                    if ( $dataHash{'pLI_ExAC'} eq ""){
                        $pLI_Color = '#FFFFFF';
                        last;
                    }
                    if ( $dataHash{'pLI_ExAC'}  >= $pli){
                        $pLI_Color = $pLI_ColorHash{$pli};
                        last;
                    }
                }
			}else{
                $pLI_Color = '#FFFFFF';            
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
				
				$hashFinalSortData{$dataHash{"AnnotSV ranking"}.$scorePenalty."_".$variantID}{$count}{'hashComments'}{$NameColHash{$fieldCom}} = $dataCommentHash{$fieldCom}{'values'} ; 
				#TODO check color for gene according to output position of gene field 
			}

            #assign color to Gene Name
			$hashFinalSortData{$dataHash{"AnnotSV ranking"}.$scorePenalty."_".$variantID}{$count}{'hashColor'}{$NameColHash{'Gene name'}-1} = $pLI_Color ;
            #print $NameColHash{'Gene name'}."\n";

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
\n<body>\n\n";

#table and columns names


my $htmlALL= "<div id=\"FULL+SPLIT\" class=\"tabcontent\">";
$htmlALL .= "\n\t<table id='tabFULLSPLIT' class='display' >
        \n\t\t<thead><tr>";


my $htmlFULL="<div id=\"FULL\" class=\"tabcontent\">";
$htmlFULL .= "\n\t<table id='tabFULL' class='display' >
        \n\t\t<thead><tr>";

foreach my $col (sort {$a <=> $b} keys %OutColHash){
			#print HTML "\t<th style=\"word-wrap: break-word\"   >";
			$htmlALL .= "\t<th >".$OutColHash{$col}."\t</th>\n";
			$htmlFULL .= "\t<th >".$OutColHash{$col}."\t</th>\n";
}


$htmlALL .= "</tr>\n</thead>\n<tbody>\n";
$htmlFULL .= "</tr>\n</thead>\n<tbody>\n";
					
					

my $htmlEndTable = "</tbody>\n</table>\n</div>\n\n\n";
					
					
					


my $htmlEnd = "<div class=\"tab\">\n
<button class=\"tablinks\" onclick=\"openCity(event, 'FULL')\">FULL</button>\n
<button class=\"tablinks\" onclick=\"openCity(event, 'FULL+SPLIT')\">FULL+SPLIT</button>\n
</div>";


$htmlEnd .= "\n</body>\n</html>";






open(HTML, '>', $outDir."/".$outPrefix."annotSV.html") or die $!;


#########################################################################
#################### Sort by MPA ranking for the output

#create user friendly ranking score
my $kindRank=0;

# check if last line was "FULL" to finish the raw with </tr>
my $FULLboolean=0;


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


#Once for "ALL"  = FULL+SPLIT
for( my $fieldNbr = 0 ; $fieldNbr < scalar @{$hashFinalSortData{$rank}{$variant}{'finalArray'}} ; $fieldNbr++){
    #print $fieldNbr."\n";
    #print $hashFinalSortData{$rank}{$variant}{'finalArray'}[$fieldNbr]."\n";



#foreach my $field ( @{$hashFinalSortData{$rank}{$variant}{'finalArray'}} ){
		#print $rank."bis\n";
        #print HTML "\t<td style=\"word-wrap: break-word\";class=\"tooltip\"; title=".$field."  >";
        #print HTML "\t<td >";
	
    #implementing rowspan only if line is "full" type
    # TODO loop on finalArray with counter to skip first td SVID for split line
    #if (     $hashFinalSortData{$rank}{$variant}{'finalArray'}[$NameColHash{'AnnotSV type'} - 1] eq "full") {
       
		
       if (defined $hashFinalSortData{$rank}{$variant}{'hashColor'}{$fieldNbr}){
            $htmlALL .= "\t<td bgcolor=\"".$hashFinalSortData{$rank}{$variant}{'hashColor'}{$fieldNbr}."\">";
        
       }else{
            $htmlALL .= "\t<td >";
       }
       #$htmlALL .= "\t<td rowspan=\".$hashFinalSortData{$rank}{$variant}{SVIDNbr}."\">";
    #}




	#if (defined $field){
	if (defined $hashFinalSortData{$rank}{$variant}{'finalArray'}[$fieldNbr]){
	#print HTML $field;
	$htmlALL .= "<div class=\"tooltip\">".$hashFinalSortData{$rank}{$variant}{'finalArray'}[$fieldNbr];
    $htmlALL .= "<span class=\"tooltiptext tooltip-bottom\">".$hashFinalSortData{$rank}{$variant}{'finalArray'}[$fieldNbr];
    # add comments
        if (defined $hashFinalSortData{$rank}{$variant}{'hashComments'}{$fieldNbr}){
            $htmlALL .= $hashFinalSortData{$rank}{$variant}{'hashComments'}{$fieldNbr}."</span></div>";   
        }else{
            $htmlALL .= "</span></div>";
        }

    #	$hashFinalSortData{$rank}{$variant}{'hashComments'}{$NameColHash{$fieldCom}} 

    }else {
		#print HTML "."
        $htmlALL .= "."
	}
    #print HTML "\t</td>\n";
	$htmlALL .= "\t</td>\n";
}

#Once for "FULL"  = FULL only
for( my $fieldNbr = 0 ; $fieldNbr < scalar @{$hashFinalSortData{$rank}{$variant}{'finalArray'}} ; $fieldNbr++){

    if ($hashFinalSortData{$rank}{$variant}{'finalArray'}[$NameColHash{'AnnotSV type'} - 1] eq "split") {
		$FULLboolean = 0;
		next;
	}elsif ($FULLboolean != 1){
			$htmlFULL .= "<tr>\n";
			$FULLboolean = 1;
	} 


       
	#assign color	
	if (defined $hashFinalSortData{$rank}{$variant}{'hashColor'}{$fieldNbr}){
		$htmlFULL .= "\t<td bgcolor=\"".$hashFinalSortData{$rank}{$variant}{'hashColor'}{$fieldNbr}."\">";    
	}else{
		$htmlFULL .= "\t<td >";
	}

	#assign values and tooltip value
	if (defined $hashFinalSortData{$rank}{$variant}{'finalArray'}[$fieldNbr]){
		$htmlFULL .= "<div class=\"tooltip\">".$hashFinalSortData{$rank}{$variant}{'finalArray'}[$fieldNbr];
		$htmlFULL .= "<span class=\"tooltiptext tooltip-bottom\">".$hashFinalSortData{$rank}{$variant}{'finalArray'}[$fieldNbr];
    # add comments
        if (defined $hashFinalSortData{$rank}{$variant}{'hashComments'}{$fieldNbr}){
            $htmlFULL .= $hashFinalSortData{$rank}{$variant}{'hashComments'}{$fieldNbr}."</span></div>";   
        }else{
            $htmlFULL .= "</span></div>";
        }


    }else {
        $htmlFULL .= "."
	}
	$htmlFULL .= "\t</td>\n";
}




	#end of table line
	$htmlALL .= "</tr>\n";
	if ($FULLboolean == 1){
		$htmlFULL .= "</tr>\n";
		$FULLboolean = 0;
	}
	
	




##############################################################
###########################      ALL     #####################

	
	}   # END FOREACH VARIANT
}	#END FOREACH RANK



close(VCF);




print HTML $htmlStart;
print HTML $htmlFULL;
print HTML $htmlEndTable;
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
