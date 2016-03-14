package gene_library;
use Moose;
use strict;
use warnings;
use KEGG_Pathway;
use LWP::Simple;
use gene;
use JSON;
use network;
use Data::Dumper qw(Dumper);
use go;
has 'genes' => (
	is => 'rw',
	isa => 'HashRef[gene]',
	required => 0
);
has 'networks' => (
	is => 'rw',
	isa => 'HashRef[network]'
);
has 'number_of_networks' => (
	is => 'rw',
	isa => 'Int'
);
sub record_gene_info {
	if (@_) {
		my ($self, $filein) = @_ ;
		if (open FILE, ">$filein") {
			print "File opened";
			print FILE "UNIPROT\tKEGG\n";
			foreach my $line (sort keys %{$self->genes}) {
				print FILE $self->genes->{$line}->return_uniprot, "\t";
				print FILE $self->genes->{$line}->return_keggid, "\n";
			}

		} else {
			die "ERROR : Couldn't get 2 file\n";
		}
	}
}
sub open_gen_info {
	if (@_) {
		my ($self, $filein) = @_;
		if (open FILE, $filein) {
			my @file = <FILE>;
			my %gene_hash	;
			print "getting information from $filein\n";
			shift(@file)	;
			foreach my $line (@file) {
				my ($uniprot_id, $kegg_id) = split("\t", $line);
				$gene_hash{$uniprot_id} = gene -> new(
					'UNIPROT_ID' => $uniprot_id,
					'KEGG_ID' => $kegg_id
				);
			}
			$self -> genes(\%gene_hash);
		}
	}
}
sub input_genes_from_kegg {
	#
	#	This function allows us to create gene instances inside our gene
	#	library, which will be used for develeoping gene-interaction
	#	networks.
	#
	if (@_) {
		#
		#	input genes loads KEGG Identifiersfrom file, and looks for
		#	it through KEGG Rest Interface.
		#
		my ($self, $path) = @_;
		if ($path) {
			if (open(INFILE, $path)) {
				print "Reading information from KEGG-identifiers file\n";
				my $page					;
				my @kegg_file = <INFILE> 	;
				my $json = JSON->new		;
				my $ref_content				;
				my $refs					;
				my $uniprotreff				;
				my %hashgene				;
				my $seq						;
				my $kegg_pathway_name		;
				my $kegg_pathway_id			;
				my $iter = 0				;
				foreach my $keggref (@kegg_file) {
					#
					#
					#
					my @kegg_pathway_object_array;
					$iter ++;
					chomp($keggref);
					print "[$iter]\tworking in : fetching information for $keggref\n";
					$page = get ('http://togows.dbcls.jp/entry/genes/ath:'.$keggref.'/dblinks.json')
						or next;
					$ref_content	=	$json->decode($page);
					$refs			=	$ref_content -> [0]	;
					($uniprotreff	=	$refs->{UniProt}->[0]) or next	;
					$page = get ('http://togows.dbcls.jp/entry/genes/ath:'.$keggref.'/aaseq.json')
						or next;
					$ref_content	=	$json->decode($page);
					$seq			=	$ref_content->[0];
					$page = get('http://togows.dbcls.jp/entry/genes/ath:'.$keggref.'/pathways.json')
						or next;
					$ref_content	=	$json->decode($page);
					$refs			=	$ref_content -> [0];
					foreach my $kegg_accesion (keys %{$refs}) {
						#
						#
						#
						push (@kegg_pathway_object_array, KEGG_Pathway->new(
								'KEGG_id' => $kegg_accesion,
								'KEGG_name' => $refs->{$kegg_accesion}
							)
						) or print "ERROR : Unable to create KEGG anotation for $kegg_accesion\n";
					}


					$hashgene{$uniprotreff} = gene-> new(
						'KEGG_ID' => $keggref,
						'UNIPROT_ID' => $uniprotreff,
						'SEQ' => $seq,
						'PATHWAY' => \@kegg_pathway_object_array


					) or (print $uniprotreff." not avalaible\n" and next);
				}
				$self -> genes(\%hashgene) ;
				print "saving data to gene_information_file.txt\n";
				$self -> record_gene_info('gene_information_file.txt');
			}
		} else {
			die "ERROR : No path\n";
		}

	} else {
		die "ERROR : No input\n";
	}
}
sub look4goanotation {
	if (@_) {
		#
		#	This function pourpose is to get process antotation
		#	from the GO
		#
		my ($self) = @_ ;
		my $json = JSON->new;
		if ($self->networks) {
			print "working in : looking for GO anotation\n";
			foreach my $netname (keys %{$self -> networks}) {
				my %go_network_asociated;
				foreach my $interaction_name (keys %{$self->networks->{$netname}->interaction}) {

					my @interaction_terms = split (":",$interaction_name);
					while (@interaction_terms) {
						my %go_gene_asociated;
						my $current_gene = shift(@interaction_terms);
						if ($self->genes->{$current_gene}->GO) {
							%go_gene_asociated = %{$self->genes->{$current_gene}->GO};
						}
						print "\tworking in : looking gor GO anotation in $netname:$current_gene\n";
						my $page = get ('http://togows.dbcls.jp/entry/uniprot/'.$current_gene.'/dr.json');
						if ($page) {
							my $ref_content = $json->decode($page);
							my $refs = $ref_content->[0] ;
							foreach my $go (@{$refs->{'GO'}}){
								foreach my $go_entry_ref ($go) {
									if (exists $go_gene_asociated{$go_entry_ref->[0]}) {
										next;
									} else {
										if ($go_entry_ref->[1] =~ /^P:(.*)/) {
											print "\t\tadding a new ontology term to network\n";
											$go_gene_asociated{$go_entry_ref->[0]} = go -> new(
												'GO_id' => $go_entry_ref->[0],
												'GO_category' => 'P',
												'GO_term' => $1
											);
											$go_network_asociated{$go_entry_ref->[0]} = $go_gene_asociated{$go_entry_ref->[0]};
										}
									}
								}
							}
						}
						$self->genes->{$current_gene}->GO(\%go_gene_asociated);
					}

				}
				$self->networks->{$netname}->go_anotation(\%go_network_asociated);
			}
		} else {
			print "ERROR : No networks\n";
		}
	}
}
sub createnetworks {
	#
	#	This method provides the genes not present in the original
	#	dataset, and calls the subrutines in charge of linking
	#	interactions.
	#
	if (@_) {
		my ($self) = @_ ;
		my $interaction_result;
		my @interactors;
		my %interactions;
		my $deepness = 1;
		my $json = JSON -> new;
		foreach my $gene (sort keys %{$self->genes}) {
			print "working in : looking for interactions involving $gene\n";
			$interaction_result = $self-> look4interactions ($gene, $deepness);

			if (defined $interaction_result) {
				print "\tinteraction found\n";
				@interactors = split("\t", $interaction_result);
				print '|'.$interaction_result.'|', "\n";
				my $error_control = 0;
				foreach my $interactor (@interactors) {

					if (not exists $self->genes->{$interactor}) {
						my %hashgenes;
						if (%{$self->genes}) {
							%hashgenes = %{$self->genes};
						}
						print "\tworking in : adding intermediate node\n";
						print "\t\tworking in : looking for sequence\n";
						my $page;
						if (not $page= get ('http://togows.dbcls.jp/entry/uniprot/'.$interactor.'/seq.json') ) {
							print "ERROR : Couldn't get to sequence\n";
							$error_control=1;
							next;
						}
						my $ref_content	=	$json->decode($page);
						my $seq			=	$ref_content->[0];
						print "\t\tworking in : looking for KEGG id\n";
						if (not $page = get ('http://togows.dbcls.jp/entry/uniprot/'.$interactor.'/dr.json') ) {
							$error_control=1;
							print "ERROR : Couldn't get to togo rest interface\n";
							next;
						}
						$ref_content	=	$json->decode($page);
						my $ref			=	$ref_content->[0];
						my $keggref;
						if (not $keggref = $ref->{KEGG}->[0]->[0] ) {
							$error_control=1;
							print "ERROR : Unable to find KEGG id\n";
							next;
						}
						print "\t\t\tkegg id = $keggref\n";
						($keggref) and ($keggref =~ s/^ath\://);
						if (not $page = get('http://togows.dbcls.jp/entry/genes/ath:'.$keggref.'/pathways.json')) {
							$error_control=1;
							print "ERROR : Unable to fing KEGG pathway id\n";
							next;
						}
						$ref_content	=	$json->decode($page);
						$ref			=	$ref_content -> [0];
						my @kegg_pathway_object_array;
						foreach my $kegg_accesion (keys %{$ref}) {
							push (@kegg_pathway_object_array, KEGG_Pathway->new(
									'KEGG_id' => $kegg_accesion,
									'KEGG_name' => $ref->{$kegg_accesion}
								)
							) or print "ERROR : Unable to create KEGG anotation for $kegg_accesion\n";
						}
						if ($error_control==0) {

							$hashgenes{$interactor} = gene-> new(
								'UNIPROT_ID' => $interactor,
								'SEQ' => $seq,
								'KEGG_ID' => $keggref,
								'PATHWAY' => \@kegg_pathway_object_array
							) or die "ERROR : Couldn't create object gene for $interactor\n";
							$self -> genes(\%hashgenes);
						}
					}
				}
				print "\t\t\terror control value = $error_control";
				if ($error_control == 0) {


					while (scalar(@interactors) > 1) {
						my $gen1 = shift @interactors ;
						$interactions{$gen1.':'.$interactors[0]} = interaction -> new(
							'gene1' =>	$self->genes->{$gen1},
							'gene2' =>	$self->genes->{$interactors[0]},
							'id' =>		$gen1.':'.$interactors[0]
						) and print "interaction added\n";
					}

				} else {
					print "\tinteraction not added.\n";
				}
			}
		}

		$self -> joininteractions(%interactions);
		if ($self->networks) {
			foreach my $network (keys %{$self->networks}) {
				print $self->networks->{$network}->network_id, "\n";
				foreach my $interaction (keys %{$self->networks->{$network}->interaction}) {
					print "\t",$self->networks->{$network}->interaction->{$interaction}->id, "\n" ;
				}
			}
		}
	}
}

sub joininteractions {
	#
	#	This method calls recursively itself to link all the diferent
	#	interactions in a network-
	#
	my ($self, %interaction_hash) = @_ 	;

	if (%interaction_hash) {
		print "working in : identifying possible gene interaction network\n";

		my $door = 1		;
		my @temp_list = keys(%interaction_hash)		;
		my $network_tupple = $temp_list[0]			;
		my %network_final_terms 					;

		$network_final_terms{$network_tupple}	=	$interaction_hash{$network_tupple};
		delete $interaction_hash{$network_tupple}	;
		while ($door == 1) {
			my @network_terms = split(":", $network_tupple)	;
			$door = 0;
			foreach my $current_interaction (sort keys %interaction_hash) {

				if (defined $interaction_hash{$current_interaction}->id) {
					foreach my $current_term (@network_terms) {
						my @current_interaction_identifiers = split(":", $current_interaction);
						if ($current_term eq $current_interaction_identifiers[0]){
							$network_tupple = $network_tupple.":".$current_interaction_identifiers[1];
							$network_final_terms{$current_interaction} = $interaction_hash{$current_interaction};

							delete $interaction_hash{$current_interaction};
							$door = 1;
						} elsif ($current_term eq $current_interaction_identifiers[1]) {
							$network_tupple = $network_tupple.":".$current_interaction_identifiers[0];
							$network_final_terms{$current_interaction} = $interaction_hash{$current_interaction};

							delete $interaction_hash{$current_interaction};
							$door = 1;
						}
					}
				} else {
					delete $interaction_hash{$current_interaction};
				}
			}
		}
		print "\tnew network found\n";
		my $iter=0;
		foreach my $line (sort keys %network_final_terms) {
			if (not defined $network_final_terms{$line}) {
				delete $network_final_terms{$line};
			}
		}
		my %built_networks;
		if (defined $self->networks) {
			%built_networks = %{$self->networks};
		}
		my $number			=	$self->number_of_networks + 1;
		my $network_name	=	"N$number";

		print "\t\twhose name will be $network_name\n";
		$built_networks{$network_name} = network -> new (
			'network_id' => $network_name,
			'interaction' => \%network_final_terms
		);

		$self -> networks(\%built_networks);
		$self -> number_of_networks($number);
		if (%interaction_hash) {
			$self-> joininteractions(%interaction_hash);
		}
		return;
	}
}
sub look4interactions {
	if (@_) {
		#
		#	Recursive Function
		#

		my ($self, $currentgene, $deepness)	=	@_;

		if ($deepness == 3) {
			if (exists $self->genes-> {$currentgene}) {

				return $currentgene ;
			} else {
				return	;
			}
		} else {
			my $ebi_url =  'http://www.ebi.ac.uk/';
			my $psicquic_tool_url = 'Tools/webservices/psicquic/intact/webservices/current/search/';
			my $method	=	'interactor/';
			my $format	=	'format=tab27';
			my $rest_url				;
				$rest_url = $rest_url = $ebi_url.
							$psicquic_tool_url.
							$method.
							$currentgene.'?'.
							$format	;
			my $rest_page	=	get($rest_url) or
				return;
			my @interaction_list	;
			my $focusgene			;
			my $recursive_answer	;
			my @instances			;
				@instances	=	split("\n", $rest_page) ;
			my @fields				;
			my $answer				;

			foreach my $inst (@instances) {

				@fields = split("\t", $inst);
				my ($idA, $idB, $altIdA, $altIdB, $aliasA, $aliasB,
				$detMethod, $author, $pub, $orgA, $orgB, $intType,
				$sourceDb, $intId, $conf) = @fields ;
				$idA =~ s/uniprotkb://g;
				$idB =~ s/uniprotkb://g;
				($idA eq $idB)
					and next ;
				($idA  eq $currentgene)
					and ($focusgene = $idB);
				($idB  eq $currentgene)
					and ($focusgene = $idA);
				(defined $focusgene)
					or die "ERROR : Couldn't get focus gene value\n";
				#print $indenter."FOCUS GENE = ",$focusgene, "\n";
				if (exists $self->genes->{$focusgene}) {
					return $currentgene."\t".$focusgene ;
				} else {
					$recursive_answer = $self -> look4interactions ($focusgene, $deepness + 1);
					if (defined $recursive_answer) {
						($recursive_answer =~ m/$currentgene/)
							and next;
						($recursive_answer =~ m/$currentgene/)
							or $answer = $currentgene."\t".$recursive_answer;
						defined $answer
							and return $answer ;
					}
				}

			}
		}
	} else {
		die "ERROR : No information to work with\n";
	}
}
# sub look4interactions {
# 	if (@_) {
# 		#
# 		#	This method links genes inside the dataset using, if necessary,
# 		#	genes not present in the dataset, through recursive calling.
# 		#
#
# 		my ($self, $currentgene, $deepness)	=	@_;
#
# 		if ($deepness == 3) {
# 			if (exists $self->genes-> {$currentgene}) {
#
# 				return $currentgene ;
# 			} else {
# 				return	;
# 			}
# 		} else {
# 			my $ebi_url =  'http://www.ebi.ac.uk/';
# 			my $psicquic_tool_url = 'Tools/webservices/psicquic/intact/webservices/current/search/';
# 			my $method	=	'interactor/';
# 			my $format	=	'format=tab27';
# 			my $rest_url				;
# 				$rest_url = $rest_url = $ebi_url.
# 							$psicquic_tool_url.
# 							$method.
# 							$currentgene.'?'.
# 							$format	;
# 			my $rest_page	=	get($rest_url) or
# 				return;
# 			my @interaction_list	;
# 			my $focusgene			;
# 			my $recursive_answer	;
# 			my @instances			;
# 				@instances	=	split("\n", $rest_page) ;
# 			my @fields				;
# 			my $answer				;
#
# 			foreach my $inst (@instances) {
#
# 				@fields = split("\t", $inst);
# 				my ($idA, $idB, $altIdA, $altIdB, $aliasA, $aliasB,
# 				$detMethod, $author, $pub, $orgA, $orgB, $intType,
# 				$sourceDb, $intId, $conf) = @fields ;
# 				$idA =~ s/uniprotkb://g;
# 				$idB =~ s/uniprotkb://g;
# 				($idA eq $idB)
# 					and next ;
# 				($idA  eq $currentgene)
# 					and ($focusgene = $idB);
# 				($idB  eq $currentgene)
# 					and ($focusgene = $idA);
# 				(defined $focusgene)
# 					or die "ERROR : Couldn't get focus gene value\n";
# 				#print $indenter."FOCUS GENE = ",$focusgene, "\n";
# 				if (exists $self->genes->{$focusgene}) {
# 					return $currentgene."\t".$focusgene ;
# 				} else {
# 					$recursive_answer = $self -> look4interactions ($focusgene, $deepness + 1);
# 					if (defined $recursive_answer) {
# 						($recursive_answer =~ m/$currentgene/)
# 							and next;
# 						($recursive_answer =~ m/$currentgene/)
# 							or $answer = $currentgene."\t".$recursive_answer;
# 						defined $answer
# 							and return $answer ;
# 					}
# 				}
#
# 			}
# 		}
# 	} else {
# 		die "ERROR : No information to work with\n";
# 	}
# }
sub report_interactions {
	#
	#	Creates a report.
	#
	if (@_) {
		my ($self, $fileout) = @_ ;
		if (open OUTILE, '>'.$fileout) {

			if ($self -> networks) {


				print OUTILE "INTERACTION REPORT\n";

				foreach my $network (sort keys %{$self->networks}) {
					print OUTILE "\tNETWORK : $network\n";
					foreach my $interaction (sort keys %{$self->networks->{$network}->interaction}) {
						print OUTILE "\t\tINTERACTION : $interaction\n";
						print OUTILE "\t\t\tGENE ID : ".$self->networks->{$network}->interaction->{$interaction}->gene1->UNIPROT_ID."\n";
						if ($self->networks->{$network}->interaction->{$interaction}->gene1->PATHWAY){
							print OUTILE "\t\t\t\tKEGG PATHWAY ANOTATION : \n";
							foreach my $keggpathway (@{$self->networks->{$network}->interaction->{$interaction}->gene1->PATHWAY}){
								print OUTILE "\t\t\t\t\t", $keggpathway->KEGG_name, "\n";
							}
						}
						if ($self->networks->{$network}->interaction->{$interaction}->gene1->GO) {
							print OUTILE "\t\t\t\tGO ONTOLOGY ANOTATION : \n";
							foreach my $goterm (sort keys %{$self->networks->{$network}->interaction->{$interaction}->gene1->GO}) {
								print OUTILE "\t\t\t\t\t", $self->networks->{$network}->interaction->{$interaction}->gene1->GO->{$goterm}->GO_id, " : ";
								print OUTILE  $self->networks->{$network}->interaction->{$interaction}->gene1->GO->{$goterm}->GO_term, "\n";
							}
						}
						print OUTILE "\t\t\tGENE ID : ".$self->networks->{$network}->interaction->{$interaction}->gene2->UNIPROT_ID."\n";
						if ($self->networks->{$network}->interaction->{$interaction}->gene2->PATHWAY){
							foreach my $keggpathway (@{$self->networks->{$network}->interaction->{$interaction}->gene2->PATHWAY}){
								print OUTILE "\t\t\t\t", $keggpathway->KEGG_name, "\n";
							}
						}
						if ($self->networks->{$network}->interaction->{$interaction}->gene2->GO) {
							print OUTILE "\t\t\t\tGO ONTOLOGY ANOTATION : \n";
							foreach my $goterm (sort keys %{$self->networks->{$network}->interaction->{$interaction}->gene2->GO}) {
								print OUTILE "\t\t\t\t\t", $self->networks->{$network}->interaction->{$interaction}->gene2->GO->{$goterm}->GO_id, " : ";
								print OUTILE  $self->networks->{$network}->interaction->{$interaction}->gene2->GO->{$goterm}->GO_term, "\n";
							}
						}

					}
					print OUTILE "\t\tGO TERMS :\n";
					foreach my $goterm (sort keys %{$self -> networks->{$network}->go_anotation}) {
						print OUTILE "\t\t\tGO : ".$self->networks->{$network}->go_anotation->{$goterm}->GO_id, "\n";
						print OUTILE "\t\t\t\tProces : ".$self->networks->{$network}->go_anotation->{$goterm}->GO_term, "\n";
					}
				}
			}
		}
	}
}
1;
