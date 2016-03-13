package gene_library;
use KEGG_Pathway;
use Moose;
use strict;
use warnings;
use LWP::Simple;
use gene;
use JSON;
use network;
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

sub input_genes_from_kegg {

	if (@_) {
		my ($self, $path) = @_;
		if ($path) {
			if (open(INFILE, $path)) {
				#
				#
				#
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
					print "[$iter]\tworking in $keggref\n";
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
					foreach my $kegg_accesion (keys $refs) {
						#
						#
						#
						push (@kegg_pathway_object_array, KEGG_Pathway->new(
								'KEGG_id' => $kegg_accesion,
								'KEGG_name' => $refs->{$kegg_accesion}
							)
						);
					}


					$hashgene{$uniprotreff} = gene-> new(
						'KEGG_ID' => $keggref,
						'UNIPROT_ID' => $uniprotreff,
						'SEQ' => $seq,
						'PATHWAY' => \@kegg_pathway_object_array


					) or (print $uniprotreff." not avalaible\n" and next);
				}
				$self -> genes(\%hashgene) ;
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
		my ($self) = @_ ;
		my $json = JSON->new;
		foreach my $netname (keys %{$self -> networks}) {
			foreach my $interaction_name (keys %{$self->networks->{$netname}}) {
				my @interaction_terms = split (":",$interaction_name);
				while (@interaction_terms) {
					my $page = get ('http://togows.dbcls.jp/entry/uniprot/'.shift(@interaction_terms).'/dr.json');
					my $ref_content = $json->decode($page);
					my $refs = $ref_content -> [0] ;
					foreach my $go_ref (@{$refs->{GO}}) {
						print $go_ref, "\n";
					}
				}

			}
		}
	}
}
sub createnetworks {
	if (@_) {
		my ($self) = @_ ;
		my $interaction_result;
		my @interactors;
		my %interactions;
		my $deepness = 1;
		my $json = JSON -> new;
		foreach my $gene (sort keys %{$self->genes}) {
			$interaction_result = $self-> look4interactions ($gene, $deepness);
			if (defined $interaction_result) {
				@interactors = split("\t", $interaction_result);
				print join ("\t", @interactors), "\n";
				foreach my $interactor (@interactors) {
					if (exists $self->genes->{$interactor}) {
						my %hashgenes;
						if ($self->genes) {
							%hashgenes = %{$self->genes};
						}
						my $page = get ('http://togows.dbcls.jp/entry/uniprot/'.$interactor.'/seq.json')
							or next;
						my $ref_content	=	$json->decode($page);
						my $seq			=	$ref_content->[0];

						$hashgenes{$interactor} = gene-> new(
							'UNIPROT_ID' => $interactor,
							'SEQ' => $seq
						) ;

						$self -> genes(\%hashgenes);
					}
				}
				while (scalar(@interactors) > 1) {
					my $gen1 = shift @interactors ;
					#print scalar @interactors, "\n";
					$interactions{$gen1.':'.$interactors[0]} = interaction -> new(
						'gen1' =>	$gen1,
						'gen2' =>	$interactors[0],
						'id' =>		$gen1.':'.$interactors[0]
					);
					print "Interaccion añadida\n";
				}
			}
		}
		$self -> joininteractions(\%interactions);
		foreach my $network (keys %{$self->networks}) {
			print $self->networks->{$network}->network_id, "\n";
			foreach my $interaction (keys %{$self->networks->{$network}->interaction}) {
				print "\t".$self->networks->{$network}->interaction->{$interaction}->id, "\n" ;
			}
		}
	}
}
sub joininteractions {

	my ($self, %interaction_hash) = @_ 	;
	if (%interaction_hash) {
		print "Geetting into\n";
		my $door = 1		;
		my @temp_list = keys(%interaction_hash)		;
		my $network_tupple = $temp_list[0]			;
		my %network_final_terms 					;
		$network_final_terms{$network_tupple}	=	$interaction_hash{$network_tupple};
		delete $interaction_hash{$network_tupple}	;
		while ($door == 1) {
			print "In\n";
			my @network_terms = split(":", $network_tupple)	;
			print "\t$network_tupple\n";
			$door = 0;
			foreach my $current_interaction (sort keys %interaction_hash) {
				#
				#
				#
				print "\tworking in $current_interaction\n";
				foreach my $current_term (@network_terms) {
					my @current_interaction_identifiers = split(":", $current_interaction);
					print "\t\t\t$current_term\n";
					if ($current_term eq $current_interaction_identifiers[0]){
						print "Network expanded. Term = $current_interaction_identifiers[1]\n";
						$network_tupple = $network_tupple.":".$current_interaction_identifiers[1];
						$network_final_terms{$current_interaction} = $interaction_hash{$current_interaction};
						delete $interaction_hash{$current_interaction};
						$door = 1;
					} elsif ($current_term eq $current_interaction_identifiers[1]) {
						print "Network expanded. Term = $current_interaction_identifiers[0]\n";
						$network_tupple = $network_tupple.":".$current_interaction_identifiers[0];
						$network_final_terms{$current_interaction} = $interaction_hash{$current_interaction};
						delete $interaction_hash{$current_interaction};
						$door = 1;
					}
				}
			}
		}
		print "NEW NETWORK\n";
		my $iter=0;
		foreach my $line (sort keys %network_final_terms) {
			print "[$iter]\t";
			print $line;
			print "\n";
			$iter ++;
		}
		my %built_networks;
		if (defined $self->networks) {
			#print "Coming to preexistent networks\n";
			%built_networks = %{$self->networks};
		}
		my $number			=	$self->number_of_networks + 1;
		my $network_name	=	"N$number";

		$built_networks{$network_name} = network -> new (
			'network_id' => $network_name
			#'interaction' => %network_final_terms
		);
		$built_networks{$network_name} -> interaction(\%built_networks);
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
1;
