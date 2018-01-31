package MaranasTools::MaranasToolsClient;

use JSON::RPC::Client;
use POSIX;
use strict;
use Data::Dumper;
use URI;
use Bio::KBase::Exceptions;
my $get_time = sub { time, 0 };
eval {
    require Time::HiRes;
    $get_time = sub { Time::HiRes::gettimeofday() };
};

use Bio::KBase::AuthToken;

# Client version should match Impl version
# This is a Semantic Version number,
# http://semver.org
our $VERSION = "0.1.0";

=head1 NAME

MaranasTools::MaranasToolsClient

=head1 DESCRIPTION


A KBase module: MaranasTools
This sample module contains one small method - filter_contigs.


=cut

sub new
{
    my($class, $url, @args) = @_;
    

    my $self = {
	client => MaranasTools::MaranasToolsClient::RpcClient->new,
	url => $url,
	headers => [],
    };

    chomp($self->{hostname} = `hostname`);
    $self->{hostname} ||= 'unknown-host';

    #
    # Set up for propagating KBRPC_TAG and KBRPC_METADATA environment variables through
    # to invoked services. If these values are not set, we create a new tag
    # and a metadata field with basic information about the invoking script.
    #
    if ($ENV{KBRPC_TAG})
    {
	$self->{kbrpc_tag} = $ENV{KBRPC_TAG};
    }
    else
    {
	my ($t, $us) = &$get_time();
	$us = sprintf("%06d", $us);
	my $ts = strftime("%Y-%m-%dT%H:%M:%S.${us}Z", gmtime $t);
	$self->{kbrpc_tag} = "C:$0:$self->{hostname}:$$:$ts";
    }
    push(@{$self->{headers}}, 'Kbrpc-Tag', $self->{kbrpc_tag});

    if ($ENV{KBRPC_METADATA})
    {
	$self->{kbrpc_metadata} = $ENV{KBRPC_METADATA};
	push(@{$self->{headers}}, 'Kbrpc-Metadata', $self->{kbrpc_metadata});
    }

    if ($ENV{KBRPC_ERROR_DEST})
    {
	$self->{kbrpc_error_dest} = $ENV{KBRPC_ERROR_DEST};
	push(@{$self->{headers}}, 'Kbrpc-Errordest', $self->{kbrpc_error_dest});
    }

    #
    # This module requires authentication.
    #
    # We create an auth token, passing through the arguments that we were (hopefully) given.

    {
	my %arg_hash2 = @args;
	if (exists $arg_hash2{"token"}) {
	    $self->{token} = $arg_hash2{"token"};
	} elsif (exists $arg_hash2{"user_id"}) {
	    my $token = Bio::KBase::AuthToken->new(@args);
	    if (!$token->error_message) {
	        $self->{token} = $token->token;
	    }
	}
	
	if (exists $self->{token})
	{
	    $self->{client}->{token} = $self->{token};
	}
    }

    my $ua = $self->{client}->ua;	 
    my $timeout = $ENV{CDMI_TIMEOUT} || (30 * 60);	 
    $ua->timeout($timeout);
    bless $self, $class;
    #    $self->_validate_version();
    return $self;
}




=head2 run_optstoic

  $output = $obj->run_optstoic($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a MaranasTools.OptStoicParams
$output is a MaranasTools.OptStoicOutput
OptStoicParams is a reference to a hash where the following keys are defined:
	model has a value which is a MaranasTools.model_upa
	start_compound has a value which is a MaranasTools.compound_id
	target_compound has a value which is a MaranasTools.compound_id
	max_steps has a value which is an int
	use_heterologous_steps has a value which is a MaranasTools.boolean
	dG_threshold has a value which is a float
	workspace_name has a value which is a string
model_upa is a string
compound_id is a string
boolean is an int
OptStoicOutput is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string

</pre>

=end html

=begin text

$params is a MaranasTools.OptStoicParams
$output is a MaranasTools.OptStoicOutput
OptStoicParams is a reference to a hash where the following keys are defined:
	model has a value which is a MaranasTools.model_upa
	start_compound has a value which is a MaranasTools.compound_id
	target_compound has a value which is a MaranasTools.compound_id
	max_steps has a value which is an int
	use_heterologous_steps has a value which is a MaranasTools.boolean
	dG_threshold has a value which is a float
	workspace_name has a value which is a string
model_upa is a string
compound_id is a string
boolean is an int
OptStoicOutput is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string


=end text

=item Description



=back

=cut

 sub run_optstoic
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function run_optstoic (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to run_optstoic:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'run_optstoic');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "MaranasTools.run_optstoic",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'run_optstoic',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method run_optstoic",
					    status_line => $self->{client}->status_line,
					    method_name => 'run_optstoic',
				       );
    }
}
 


=head2 run_steadycom

  $return = $obj->run_steadycom($params)

=over 4

=item Parameter and return types

=begin html

<pre>
$params is a MaranasTools.SteadyComParams
$return is a MaranasTools.SteadyComOutput
SteadyComParams is a reference to a hash where the following keys are defined:
	model_inputs has a value which is a reference to a list where each element is a MaranasTools.ModelInput
	medium_upa has a value which is a string
	flux_output has a value which is a string
	workspace_name has a value which is a string
ModelInput is a reference to a hash where the following keys are defined:
	model_upa has a value which is a string
	fixed_gr has a value which is a float
SteadyComOutput is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string
	flux_output has a value which is a string

</pre>

=end html

=begin text

$params is a MaranasTools.SteadyComParams
$return is a MaranasTools.SteadyComOutput
SteadyComParams is a reference to a hash where the following keys are defined:
	model_inputs has a value which is a reference to a list where each element is a MaranasTools.ModelInput
	medium_upa has a value which is a string
	flux_output has a value which is a string
	workspace_name has a value which is a string
ModelInput is a reference to a hash where the following keys are defined:
	model_upa has a value which is a string
	fixed_gr has a value which is a float
SteadyComOutput is a reference to a hash where the following keys are defined:
	report_name has a value which is a string
	report_ref has a value which is a string
	flux_output has a value which is a string


=end text

=item Description



=back

=cut

 sub run_steadycom
{
    my($self, @args) = @_;

# Authentication: required

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function run_steadycom (received $n, expecting 1)");
    }
    {
	my($params) = @args;

	my @_bad_arguments;
        (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"params\" (value was \"$params\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to run_steadycom:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'run_steadycom');
	}
    }

    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
	    method => "MaranasTools.run_steadycom",
	    params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{error}->{code},
					       method_name => 'run_steadycom',
					       data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method run_steadycom",
					    status_line => $self->{client}->status_line,
					    method_name => 'run_steadycom',
				       );
    }
}
 
  
sub status
{
    my($self, @args) = @_;
    if ((my $n = @args) != 0) {
        Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
                                   "Invalid argument count for function status (received $n, expecting 0)");
    }
    my $url = $self->{url};
    my $result = $self->{client}->call($url, $self->{headers}, {
        method => "MaranasTools.status",
        params => \@args,
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
                           code => $result->content->{error}->{code},
                           method_name => 'status',
                           data => $result->content->{error}->{error} # JSON::RPC::ReturnObject only supports JSONRPC 1.1 or 1.O
                          );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method status",
                        status_line => $self->{client}->status_line,
                        method_name => 'status',
                       );
    }
}
   

sub version {
    my ($self) = @_;
    my $result = $self->{client}->call($self->{url}, $self->{headers}, {
        method => "MaranasTools.version",
        params => [],
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(
                error => $result->error_message,
                code => $result->content->{code},
                method_name => 'run_steadycom',
            );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(
            error => "Error invoking method run_steadycom",
            status_line => $self->{client}->status_line,
            method_name => 'run_steadycom',
        );
    }
}

sub _validate_version {
    my ($self) = @_;
    my $svr_version = $self->version();
    my $client_version = $VERSION;
    my ($cMajor, $cMinor) = split(/\./, $client_version);
    my ($sMajor, $sMinor) = split(/\./, $svr_version);
    if ($sMajor != $cMajor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Major version numbers differ.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor < $cMinor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Client minor version greater than Server minor version.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor > $cMinor) {
        warn "New client version available for MaranasTools::MaranasToolsClient\n";
    }
    if ($sMajor == 0) {
        warn "MaranasTools::MaranasToolsClient version is $svr_version. API subject to change.\n";
    }
}

=head1 TYPES



=head2 boolean

=over 4



=item Description

A boolean - 0=false, 1=true
@range (0, 1)


=item Definition

=begin html

<pre>
an int
</pre>

=end html

=begin text

an int

=end text

=back



=head2 model_upa

=over 4



=item Description

An X/Y/Z style reference to an FBA model.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 compound_id

=over 4



=item Description

The id of a compound that exists either in the model or in the biochemistry.


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 OptStoicParams

=over 4



=item Description

model - the FBA model to use as a basis for modification
start_compound - the initial compound to be used as a source for the pathway
target_compound - the target compound to maximize yield for in the pathway
max_steps - the maximum number of steps to allow in the optimized pathway - any pathway
            created that has more than this number of steps is disqualified
use_heterologous_steps - allows adding
dG_threshold - a threshold free energy value to further constrain the path optimization


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
model has a value which is a MaranasTools.model_upa
start_compound has a value which is a MaranasTools.compound_id
target_compound has a value which is a MaranasTools.compound_id
max_steps has a value which is an int
use_heterologous_steps has a value which is a MaranasTools.boolean
dG_threshold has a value which is a float
workspace_name has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
model has a value which is a MaranasTools.model_upa
start_compound has a value which is a MaranasTools.compound_id
target_compound has a value which is a MaranasTools.compound_id
max_steps has a value which is an int
use_heterologous_steps has a value which is a MaranasTools.boolean
dG_threshold has a value which is a float
workspace_name has a value which is a string


=end text

=back



=head2 OptStoicOutput

=over 4



=item Description

report_name - name of the report object that gets generated.
report_ref - UPA of the report object that gets generated.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string


=end text

=back



=head2 ModelInput

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
model_upa has a value which is a string
fixed_gr has a value which is a float

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
model_upa has a value which is a string
fixed_gr has a value which is a float


=end text

=back



=head2 SteadyComParams

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
model_inputs has a value which is a reference to a list where each element is a MaranasTools.ModelInput
medium_upa has a value which is a string
flux_output has a value which is a string
workspace_name has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
model_inputs has a value which is a reference to a list where each element is a MaranasTools.ModelInput
medium_upa has a value which is a string
flux_output has a value which is a string
workspace_name has a value which is a string


=end text

=back



=head2 SteadyComOutput

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string
flux_output has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
report_name has a value which is a string
report_ref has a value which is a string
flux_output has a value which is a string


=end text

=back



=cut

package MaranasTools::MaranasToolsClient::RpcClient;
use base 'JSON::RPC::Client';
use POSIX;
use strict;

#
# Override JSON::RPC::Client::call because it doesn't handle error returns properly.
#

sub call {
    my ($self, $uri, $headers, $obj) = @_;
    my $result;


    {
	if ($uri =~ /\?/) {
	    $result = $self->_get($uri);
	}
	else {
	    Carp::croak "not hashref." unless (ref $obj eq 'HASH');
	    $result = $self->_post($uri, $headers, $obj);
	}

    }

    my $service = $obj->{method} =~ /^system\./ if ( $obj );

    $self->status_line($result->status_line);

    if ($result->is_success) {

        return unless($result->content); # notification?

        if ($service) {
            return JSON::RPC::ServiceObject->new($result, $self->json);
        }

        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    elsif ($result->content_type eq 'application/json')
    {
        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    else {
        return;
    }
}


sub _post {
    my ($self, $uri, $headers, $obj) = @_;
    my $json = $self->json;

    $obj->{version} ||= $self->{version} || '1.1';

    if ($obj->{version} eq '1.0') {
        delete $obj->{version};
        if (exists $obj->{id}) {
            $self->id($obj->{id}) if ($obj->{id}); # if undef, it is notification.
        }
        else {
            $obj->{id} = $self->id || ($self->id('JSON::RPC::Client'));
        }
    }
    else {
        # $obj->{id} = $self->id if (defined $self->id);
	# Assign a random number to the id if one hasn't been set
	$obj->{id} = (defined $self->id) ? $self->id : substr(rand(),2);
    }

    my $content = $json->encode($obj);

    $self->ua->post(
        $uri,
        Content_Type   => $self->{content_type},
        Content        => $content,
        Accept         => 'application/json',
	@$headers,
	($self->{token} ? (Authorization => $self->{token}) : ()),
    );
}



1;
