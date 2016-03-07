#!/usr/bin/perl
#7th March 2016
#
=head1 'network class'

=head1 'DESCRIPTION'

This class will allow us to store interaction networks for the
second assignment

=head1 'SYNOPSIS'

=cut

package network	;
use warnings	;
use strict		;
use Moose		;

has 'network_id' = (
	is => 'rw',
	isa => 'ArrayRef[Gene]'
);
