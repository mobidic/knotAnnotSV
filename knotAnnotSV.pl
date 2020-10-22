#!/usr/bin/perl

##############################################################################################
# knotAnnotSV 1.0                                                                            #	
#                                                                                            #
# knotAnnotSV: Creation of a customizable html file to visualize, filter                   	 # 
#                   and analyze an AnnotSV output                                            #
#                                                                                            #
# Author: Thomas Guignard 2020                                                               #
#                                                                                            #
# Copyright (C) 2020 Thomas Guignard (t-guignard@chu-montpellier.fr)                         #
#                                                                                            #
# This is part of knotAnnotSV source code.        	                                         #
#                                                                                            #
# This program is free software; you can redistribute it and/or                              #
# modify it under the terms of the GNU General Public License                                #
# as published by the Free Software Foundation; either version 3                             #
# of the License, or (at your option) any later version.                                     #
#                                                                                            #
# This program is distributed in the hope that it will be useful,                            #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                             #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                              #
# GNU General Public License for more details.                                               #
#                                                                                            #
# You should have received a copy of the GNU General Public License                          #
# along with this program; If not, see <http://www.gnu.org/licenses/>.                       #
##############################################################################################


use strict; 
use warnings;
use Getopt::Long; 
use YAML::XS 'LoadFile';
use Sort::Key::Natural qw(rnatkeysort);

#use Switch;
#use Data::Dumper;
#use Sort::Naturally;

#parameters
my $man = "USAGE : \nperl knotAnnotSV.pl 
\n--configFile <YAML config file for customizing output>
\n--annotSVfile <annotSV annotated file> 
\n--annotSVranking <annotSV ranking explanations file>
\n--outDir <output directory (default = current dir)> 
\n--outPrefix <output file prefix (default = \"\")> 
\n--genomeBuild <Genome Assembly reference (default = hg19)> 
\n--datatableDir <Local Path to dataTables directory containing css and js files (default = \"\", requires web connection)>"; 


my $configFile = ".";

my $help;
my $current_line;
my $incfile = "";
my $outDir = ".";
my $outPrefix = "";
my $annotSVranking = "";
my $datatableDir = "";


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

my $url2UCSC="";
my $url2OMIM="";
my $url2DECIPHER="";
my $url2ensembl="";

my $genomeBuild="";

GetOptions( "annotSVfile=s"		=> \$incfile,
			"configFile=s"		=> \$config,
			"annotSVranking=s"	=> \$annotSVranking,
			"outDir=s"			=> \$outDir,
			"outPrefix:s"		=> \$outPrefix,
			"datatableDir=s"	=> \$datatableDir,
			"genomeBuild=s"		=> \$genomeBuild,
			"help|h"			=> \$help);
				
				


#check mandatory arguments
if(defined $help || $incfile eq ""){
	die("$man");
}

#add underscore to output prefix
if($outPrefix ne ""){
	$outPrefix .= "_";
}

#check datatable pathname
if ($datatableDir ne "" && ! -d $datatableDir){
    print "\nError with datatableDir argument: ".$datatableDir." is not an existing directory.\n\n";
    exit 1;    
}
			
#check if genomeBuild exists (hg19 grch37 hg38, mm etc..) and put hg19 as default build 
if($genomeBuild eq ""){
	$genomeBuild = "hg19";
}


#open( CONFIG, "<$config" ) or die ("Cannot open vcf file $config") ;

#CONFIG FILE PARSING
print STDERR "Parsing config file....\n";
my %configHash;
my $configHash;

$configHash = LoadFile($config);
	#print Dumper($configHash);


foreach my $item (keys %{$configHash}){
		#print $item."\n";
		my $field = $item;
		#print $configHash->{$field}->{POSITION}."\n";
	
        if ($configHash->{$field}->{POSITION} != 0 && defined $OutColHash{$configHash->{$field}->{POSITION}}){
            print "\nError in config file: 'POSITION : ".$configHash->{$field}->{POSITION}."' is assigned twice.\nPlease correct config file to avoid 2 identical positions.\n\n";
            exit 1;
        }else {
		    if ($configHash->{$field}->{POSITION} != 0){
                $OutColHash{$configHash->{$field}->{POSITION}} = $field;
            }
		    $NameColHash{$field} = $configHash->{$field}->{POSITION};
        }
		
		if (defined $configHash->{$field}->{COMMENTLIST}){
			$dataCommentHash{$field}{'commentFieldList'}=$configHash->{$field}->{COMMENTLIST} ;
		}
		#print Dump($field)."\n";
}


#ANNOTSV RANKING PARSING to get ranking decision
my @SVrankLine;
my %SVrankHash;

if ($annotSVranking ne ""){
	open( SVRANK , "<$annotSVranking" )or die("Cannot open ranking file ".$annotSVranking) ;
	
	
	while( <SVRANK> ){

		chomp $_;
		@SVrankLine = split( /\t/, $_ );	

		if (defined $SVrankHash{$SVrankLine[0]."_".$SVrankLine[2]}){
			print "Doublon detected in ranking file : ".$_."\n\n";
		}else{
			$SVrankHash{$SVrankLine[0]."_".$SVrankLine[2]} = $SVrankLine[4];
		}	
	}

	#add comment with rank to AnnotSV ranking column
	if (! defined $dataCommentHash{'AnnotSV ranking'}{'commentFieldList'}){
		$dataCommentHash{'AnnotSV ranking'}{'SVrank'} = "OK";
	}

}





#Hash of ACMG incidentalome genes  => grey color #808080
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



####################################################################
#############################################
##################   Start parsing VCF

open( VCF , "<$incfile" )or die("Cannot open annotSV file ".$incfile) ;



while( <VCF> ){
  	$current_line = $_;
	$count++;


	chomp $current_line;
	@line = split( /\t/, $current_line );	
	

#############################################
##############   Treatment for First line to create header of the output

	if ( $line[0] eq "AnnotSV ID" )   {
		
		#TODO check if header contains required INFO
		print STDERR "Parsing AnnotSV header in order to get column names ... \n";

		#initialize column order
		for( my $fieldNbr = 0 ; $fieldNbr < scalar @line; $fieldNbr++){
			
			#print STDERR "Field found:  ".$line[$fieldNbr]."\n";
			$InColHash{$fieldNbr} = $line[$fieldNbr];
			
			$dataHash{$line[$fieldNbr]} = ".";
			
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
        }
        
		#reinitialize Comment Values
		foreach my $field (keys %dataCommentHash){
			delete $dataCommentHash{$field}{'values'};
		}


		#fill printable string
		for( my $fieldNbr = 0 ; $fieldNbr < scalar @line; $fieldNbr++){

			$dataHash{$InColHash{$fieldNbr}} = $line[$fieldNbr];
		}

		

		#get all comments values
		foreach my $field (keys %dataCommentHash){
            #print $field."\n";
			#if (defined $dataCommentHash{$field})
			#	foreach my $fieldCom (@{$dataCommentHash{$field}{'commentFieldList'}})
			if (defined $dataCommentHash{$field}){
                if (! defined $dataCommentHash{$field}{'values'}){
				    foreach my $fieldCom (@{$dataCommentHash{$field}{'commentFieldList'}}){
					    if (defined $dataHash{$fieldCom}){
                            $dataCommentHash{$field}{'values'} .= "<br><br><b>".$fieldCom . " :</b> ".$dataHash{$fieldCom};
                        }else{
                            print $field."\n";
                            $dataCommentHash{$field}{'values'} .= "<br><br><b>".$fieldCom . " :</b> Absent in file";
                        }
                        #print $field.":\t".$fieldCom.":\t".$dataCommentHash{'Gene name'}{'values'}."\n";
				    }
			    }
				
				# add ranking decision as comment of AnnotSV ranking field
				if (defined $dataCommentHash{$field}{'SVrank'}){
	  				if (defined $SVrankHash{$dataHash{'AnnotSV ID'}."_".$dataHash{'Gene name'}}){
                		if (! defined $dataCommentHash{$field}{'values'}){
							$dataCommentHash{$field}{'values'}= "<br><br><b>Ranking :</b> " . $SVrankHash{$dataHash{'AnnotSV ID'}."_".$dataHash{'Gene name'}} ;
						}else{
							$dataCommentHash{$field}{'values'} = "<br><br><b>Ranking :</b> " .$SVrankHash{$dataHash{'AnnotSV ID'}."_".$dataHash{'Gene name'}} . $dataCommentHash{$field}{'values'};  	
						}
					}
				}
			}
		}


		#fill finalSortData array   (try to invert foreach with NameColHash for external name)
		foreach my $field (keys %dataHash){
			if (defined $NameColHash{$field} && $NameColHash{$field} != 0){
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
        
		#change pli color if ACMG gene
		if (defined $ACMGgene{$dataHash{'Gene name'}}){
			$pLI_Color = '#808080';
		}



		# check if column names and data have the same size
		if (scalar @finalSortData != keys %OutColHash){
			print "\nError in config file: it seems that some POSITION are missing.\nPlease correct config file with continuous 'POSITION' number.\n\n"; 
			exit 1;
		}

		
		#	print Dumper(\@finalSortData);

#############################################################################
######### FILL HASH STRUCTURE WITH COMMENTS 
		
		    
            $variantID = $line[0];

			#key for good sorting should be rank concatenated with annotSV_ID
			#
			#grant full line to be first
			if ($dataHash{"AnnotSV type"} eq "full"){
				$scorePenalty = ".1";	
			}else{
				$scorePenalty = ".0";	
			}


			#finalsortData assigment
			$hashFinalSortData{$dataHash{"AnnotSV ranking"}.$scorePenalty."_".$variantID}{$count}{'finalArray'} = [@finalSortData] ; 

			
			#url to UCSC for SV , highlight in blue and zoomout x3
			$hashFinalSortData{$dataHash{"AnnotSV ranking"}.$scorePenalty."_".$variantID}{$count}{'url2UCSC'} = "http://genome.ucsc.edu/cgi-bin/hgTracks?db=".$genomeBuild."&position=chr".$dataHash{"SV chrom"}.":".$dataHash{"SV start"}."-".$dataHash{"SV end"}."&hgt.out1=submit&highlight=".$genomeBuild.".chr".$dataHash{"SV chrom"}.":".$dataHash{"SV start"}."-".$dataHash{"SV end"}."#aaedff\" target=\"_blank\" rel=\"noopener noreferrer\"" ; 

			#url to OMIM
			$hashFinalSortData{$dataHash{"AnnotSV ranking"}.$scorePenalty."_".$variantID}{$count}{'url2OMIM'} = "https://www.omim.org/entry/".$dataHash{"Mim Number"}  ; 



			#comment assigment
			foreach my $fieldCom (keys %dataCommentHash){
				$hashFinalSortData{$dataHash{"AnnotSV ranking"}.$scorePenalty."_".$variantID}{$count}{'hashComments'}{$NameColHash{$fieldCom} -1} = $dataCommentHash{$fieldCom}{'values'} ; 
				#print $dataCommentHash{'Gene name'}{'values'}."\n\n";
				#TODO check color for gene according to output position of gene field 
			}


            #assign color to Gene Name
			$hashFinalSortData{$dataHash{"AnnotSV ranking"}.$scorePenalty."_".$variantID}{$count}{'hashColor'}{$NameColHash{'Gene name'}-1} = $pLI_Color ;

			#print $NameColHash{'Gene name'}."\n";

				#ACMG
				#if(defined $ACMGgene{$finalSortData[$dicoColumnNbr{'Gene.refGene'}]} )
#				if(defined $ACMGgene{$geneName} ){
#					$hashFinalSortData{$finalSortData[$dicoColumnNbr{'MPA_ranking'}]}{$variantID}{'worksheet'} .= "#ACMG";
#				}


	} #END of IF-ELSE(#CHROM)	

}#END WHILE VCF


#prepare path to lib css and js


#original url used
#'https://code.jquery.com/jquery-3.5.1.js'
#'https://cdn.datatables.net/1.10.22/js/jquery.dataTables.min.js'
#'https://cdn.datatables.net/fixedheader/3.1.7/js/dataTables.fixedHeader.min.js'
#'https://cdn.datatables.net/1.10.22/css/jquery.dataTables.min.css'
#'https://cdn.datatables.net/fixedheader/3.1.7/css/fixedHeader.dataTables.min.css'



my $path2jquery = "https://code.jquery.com/";
my $path2jqueryDT = "https://cdn.datatables.net/1.10.22/js/";
my $path2jsFHDT = "https://cdn.datatables.net/fixedheader/3.1.7/js/";

my $path2css = "https://cdn.datatables.net/1.10.22/css/"; 
my $path2FHcss = "https://cdn.datatables.net/fixedheader/3.1.7/css/";

if ($datatableDir ne ""){

	#remove trailing slash from path
	$datatableDir =~ s/\/$//;

	$path2jquery = $datatableDir."/js/";
	$path2jqueryDT = $datatableDir."/js/";
	$path2jsFHDT = $datatableDir."/js/";

	$path2css = $datatableDir."/css/"; 
	$path2FHcss = $datatableDir."/css/";
}





############ THE TIME HAS COME TO FILL THE HTML OUTPUT FILE
####################################################################
######################### HTML  Initialisation ######################
#
#file header
my $htmlStart = "<!DOCTYPE html>\n<html>
\n<head>
\n<meta charset=\"utf-8\">
\n<title>".$outPrefix." knotAnnotSV</title>\n
\n<script type=\"text/javascript\" language=\"javascript\" src='".$path2jquery."jquery-3.5.1.js'></script>
\n<script type=\"text/javascript\" language=\"javascript\" src='".$path2jqueryDT."jquery.dataTables.min.js'></script>
\n<script type=\"text/javascript\" language=\"javascript\" src='".$path2jsFHDT."dataTables.fixedHeader.min.js'></script>
\n<script> 
\$(document).ready(function () {
 
 	\$('#tabFULL thead tr').clone(true).appendTo( '#tabFULL thead' );
            \$('#tabFULL thead tr:eq(1) th').each( function (i) {
             var title = \$(this).text();
              \$(this).html( '<input type=\"text\" placeholder=\"Search '+title+'\" />' );
                                                 
              \$( 'input', this ).on( 'keyup change', function () {
                  if ( table.column(i).search() !== this.value ) {
                        table
                       .column(i)
                       .search( this.value,true,false )
                       .draw();
                  }
            } );
        } );

	 \$('#tabFULLSPLIT thead tr').clone(true).appendTo( '#tabFULLSPLIT thead' );
            \$('#tabFULLSPLIT thead tr:eq(1) th').each( function (i) {
             var title = \$(this).text();
              \$(this).html( '<input type=\"text\" placeholder=\"Search '+title+'\" />' );
                                                 
              \$( 'input', this ).on( 'keyup change', function () {
                  if ( tableFULLSPLIT.column(i).search() !== this.value ) {
                        tableFULLSPLIT
                       .column(i)
                       .search( this.value,true,false )
                       .draw();
                  }
            } );
        } );



	var table = \$('#tabFULL').DataTable(        {\"order\": [] ,\"lengthMenu\":[ [ 50, 100, -1 ],[ 50, 100, \"All\" ]], \"fixedHeader\": true, \"orderCellsTop\": true, \"oLanguage\": { \"sLengthMenu\": \"Show _MENU_ lines\",\"sInfo\": \"Showing _START_ to _END_ of _TOTAL_ lines\" } } );
	var tableFULLSPLIT = \$('#tabFULLSPLIT').DataTable(   {\"order\": [] ,\"lengthMenu\":[ [ 50, 100, -1 ],[ 50, 100, \"All\" ]], \"fixedHeader\": true, \"orderCellsTop\": true, \"oLanguage\": { \"sLengthMenu\": \"Show _MENU_ lines\",\"sInfo\": \"Showing _START_ to _END_ of _TOTAL_ lines\" } } );

	window.onload = document.getElementById('focusFULLfirst').focus();

});

function openCity(evt, cityName) {
	var i, tabcontent, tabcontentFULL, tablinks;
	
	tabcontentFULL = document.getElementsByClassName(\"tabcontentFULL\");
	for (i = 0; i < tabcontentFULL.length; i++) {
  		tabcontentFULL[i].style.display = \"none\";	
  	}

	tabcontent = document.getElementsByClassName(\"tabcontent\");
	for (i = 0; i < tabcontent.length; i++) {
		tabcontent[i].style.display = \"none\";
	}
	tablinks = document.getElementsByClassName(\"tablinks\");
	for (i = 0; i < tablinks.length; i++) {
		tablinks[i].className = tablinks[i].className.replace(\" active\", \"\");
	}
	document.getElementById(cityName).style.display = \"block\";
	evt.currentTarget.className += \" active\";

	\$('#tabFULLSPLIT').DataTable().fixedHeader.adjust();
}





</script>
\n<link rel=\"stylesheet\" type=\"text/css\" href='".$path2css."jquery.dataTables.min.css'>
\n<link rel=\"stylesheet\" type=\"text/css\" href='".$path2FHcss."fixedHeader.dataTables.min.css'>

\n<style>
/* Style the tab */
	.tab {
		overflow: hidden;
		border: 1px solid #ccc;
		background-color: #f1f1f1;
		left: 0;
		bottom: 0;
		position: fixed;
		position: -webkit-sticky;
		position: sticky;
		}

/* Style the buttons inside the tab */
	.tab button {
		background-color: inherit;
		float: left;
		border: none;
		outline: none;
		cursor: pointer;
		padding: 14px 16px;
		transition: 0.3s;
		font-size: 17px;
		}
/* Change background color of buttons on hover */
	.tab button:hover {
		background-color: #ddd;
		}

/* Create an active/current tablink class */
	.tab button.active {
		background-color: #ccc;
		}
	.tab button:focus {
		background-color: #ccc;
		}

	thead input {
		width: 100%;
		}

/* Style the tab content */
	.tabcontent {
		display: none;
		padding: 6px 12px;
		/*border: 1px solid #ccc;*/
		border-top: none;
		}

/* Style the tab content FULL */
	.tabcontentFULL {
		display: block;
		padding: 6px 12px;
		/*border: 1px solid #ccc;*/
		border-top: none;
		}

	.tooltip {
		position: relative;
		display: inline-block;
		border-bottom: 1px dotted black;
		max-width: 300px;
		word-wrap: break-word;
		}
	
	.tooltip .tooltiptext {
		visibility: hidden;
		width: auto;
		min-width: 300px;
		max-width: 600px;
		height: auto;
		background-color: #555;
		color: #fff;
		text-align: left;
		border-radius: 6px;
		padding: 5px 0;
		position: absolute;
		z-index: 1;
		top: 100%;
		left: 20%;
		margin-left: -60px;
		opacity: 0;
		transition: opacity 0.3s;
		overflow: hidden;
		text-overflow: ellipsis;
		white-space: normal;
		overflow-wrap: break-word;
		}

	.tooltip:hover .tooltiptext {
		visibility: visible;
		opacity: 1;
		}



\n</style>

\n</head>
\n<body>\n\n";

#table and columns names


my $htmlALL= "<div id=\"FULL+SPLIT\" class=\"tabcontent\">";
$htmlALL .= "\n\t<table id='tabFULLSPLIT' class='display' >
        \n\t\t<thead><tr>";


my $htmlFULL="<div id=\"FULL\" class=\"tabcontentFULL\">";
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
<button id=\"focusFULLfirst\" class=\"tablinks\" onclick=\"openCity(event, 'FULL')\">FULL</button>\n
<button class=\"tablinks\" onclick=\"openCity(event, 'FULL+SPLIT')\">FULL+SPLIT</button>\n
</div>";


$htmlEnd .= "\n</body>\n</html>";



open(HTML, '>', $outDir."/".$outPrefix."annotSV.html") or die $!;


#########################################################################
#################### Sort by AnnotSV ranking for the output

#create user friendly ranking score
my $kindRank=0;

# check if last line was "FULL" to finish the raw with </tr>
my $FULLboolean=0;


foreach my $rank (rnatkeysort { "$_-$hashFinalSortData{$_}" } keys %hashFinalSortData){
	print $rank."\n";
	
	foreach my $variant ( keys %{$hashFinalSortData{$rank}}){

		#increse rank number then change final array
		#$kindRank++;
		
		#print $variant ."___".  $hashFinalSortData{$rank}{$variant}{'finalArray'}[$dicoColumnNbr{'MPA_ranking'}]."\n";
		#$hashFinalSortData{$rank}{$variant}{'finalArray'}[$dicoColumnNbr{'MPA_ranking'}] = $kindRank; 

		#last finalSortData assignation
		#$hashFinalSortData{$finalSortData[$dicoColumnNbr{'MPA_ranking'}]}{$variantID}{'finalArray'} = [@finalSortData] ; 
		



		#FILL tab 'ALL';
		# change bgcolor for FULL row
		if (     $hashFinalSortData{$rank}{$variant}{'finalArray'}[$NameColHash{'AnnotSV type'} - 1] eq "full") {
			$htmlALL .= "<tr style=\"background-color:teal\">\n";
		}else{
			$htmlALL .= "<tr>\n";
		}

		#Once for "ALL"  = FULL+SPLIT
		for( my $fieldNbr = 0 ; $fieldNbr < scalar @{$hashFinalSortData{$rank}{$variant}{'finalArray'}} ; $fieldNbr++){
			#print $fieldNbr."\n";
				#print $hashFinalSortData{$rank}{$variant}{'finalArray'}[$fieldNbr]."\n";

	
			  #implementing rowspan only if line is "full" type
			# TODO loop on finalArray with counter to skip first td SVID for split line
			#if (     $hashFinalSortData{$rank}{$variant}{'finalArray'}[$NameColHash{'AnnotSV type'} - 1] eq "full") {
       
					
			if (defined $hashFinalSortData{$rank}{$variant}{'hashColor'}{$fieldNbr}){
				if ($hashFinalSortData{$rank}{$variant}{'finalArray'}[$NameColHash{'location'} - 1] ne "txStart-txEnd" && $hashFinalSortData{$rank}{$variant}{'finalArray'}[$NameColHash{'location'} - 1] ne "") {
					$htmlALL .= "\t<td style=\"background: linear-gradient(-45deg,".$hashFinalSortData{$rank}{$variant}{'hashColor'}{$fieldNbr}." 50%, white 50% )\">";
				}else{
					$htmlALL .= "\t<td style=\"background-color:".$hashFinalSortData{$rank}{$variant}{'hashColor'}{$fieldNbr}."\">";
				}
        
			}else{
				$htmlALL .= "\t<td >";
			}
				#$htmlALL .= "\t<td rowspan=\".$hashFinalSortData{$rank}{$variant}{SVIDNbr}."\">";
			#}




			#if (defined $field){
			if (defined $hashFinalSortData{$rank}{$variant}{'finalArray'}[$fieldNbr]){
				#print HTML $field;

				#Add url to UCSC browser
				if ($fieldNbr eq $NameColHash{'AnnotSV ID'} - 1){
					$htmlALL .= "<div class=\"tooltip\"><a href=\"".$hashFinalSortData{$rank}{$variant}{'url2UCSC'}."\">".$hashFinalSortData{$rank}{$variant}{'finalArray'}[$fieldNbr]."</a>";
				}else{
					$htmlALL .= "<div class=\"tooltip\">".$hashFinalSortData{$rank}{$variant}{'finalArray'}[$fieldNbr];
				}
				$htmlALL .= "<span class=\"tooltiptext tooltip-bottom\"><b>".$OutColHash{$fieldNbr + 1}. " :</b> ".$hashFinalSortData{$rank}{$variant}{'finalArray'}[$fieldNbr];
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
		} #END FOR FULL+SPLIT tab


		#Once for "FULL" tab  = FULL only
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
				$htmlFULL .= "\t<td style=\"background-color:".$hashFinalSortData{$rank}{$variant}{'hashColor'}{$fieldNbr}."\">";    
			}else{
				$htmlFULL .= "\t<td >";
			}

			#assign values and tooltip value
			if (defined $hashFinalSortData{$rank}{$variant}{'finalArray'}[$fieldNbr]){
				#Add url to UCSC browser
				if ($fieldNbr eq $NameColHash{'AnnotSV ID'} - 1){
					$htmlFULL .= "<div class=\"tooltip\"><a href=\"".$hashFinalSortData{$rank}{$variant}{'url2UCSC'}."\">".$hashFinalSortData{$rank}{$variant}{'finalArray'}[$fieldNbr]."</a>";
				}else{
					$htmlFULL .= "<div class=\"tooltip\">".$hashFinalSortData{$rank}{$variant}{'finalArray'}[$fieldNbr];
				}	
				$htmlFULL .= "<span class=\"tooltiptext tooltip-bottom\"><b>".$OutColHash{$fieldNbr + 1}. " :</b> ".$hashFinalSortData{$rank}{$variant}{'finalArray'}[$fieldNbr];
				
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
		} #END FULL only tab




		#end of table line
		$htmlALL .= "</tr>\n";
		if ($FULLboolean == 1){
			$htmlFULL .= "</tr>\n";
			$FULLboolean = 0;
		}
	
	
	
	}   # END FOREACH VARIANT
}	#END FOREACH RANK



close(VCF);



#Write in HTML file
print HTML $htmlStart;
print HTML $htmlFULL;
print HTML $htmlEndTable;
print HTML $htmlALL;
print HTML $htmlEndTable;
print HTML $htmlEnd;


close(HTML);




print STDERR "\n\n\nDone!\n\n\n";

###################################################################
######################### FUNCTIONS ###############################
# Write in sheets
sub toto {
}#END OF SUB




exit 0; 


