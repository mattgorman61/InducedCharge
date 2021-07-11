%
%  Test Laplace particle FMMs in R^3
%
clear;
clc;

nsource = 10000;

source = zeros(3,nsource);

source(1,:)=rand(1,nsource);
source(2,:)=rand(1,nsource);
source(3,:)=rand(1,nsource);


% figure(1);
% scatter3(source(1,:),source(2,:),source(3,:));

%
%  timings
%

ifcharge=1;
charge = rand(1,nsource);
ifdipole=1;
dipstr = rand(1,nsource);
dipvec = rand(3,nsource);

ifpot = 1;
iffld = 1;


ntarget = 0;
target = source(:,1:ntarget);
target(1,:) = target(1,:) + 10;
[ndim,ntarget] = size(target);
%%%ntarget = 0;

%%%ntarget = 0;

ntarget;
ifpottarg = 1;
iffldtarg = 1;

if( ntarget == 0 )
ifpottarg = 0;
iffldtarg = 0;
end


% disp('')
% 'Laplace particle target FMM in R^3'

tic
iprec=2;
[U]=lfmm3dpart(iprec,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarget,target,ifpottarg,iffldtarg);
fmm_time=toc

% 'Laplace particle direct evaluation in R^3'

tic
[F]=l3dpartdirect(nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld,ntarget,target,ifpottarg,iffldtarg);
direct_time=toc


if( ifpot ), U.pot=U.pot/(4*pi); end
if( iffld ), U.fld=U.fld/(4*pi); end

if( ifpot ), F.pot=F.pot/(4*pi); end
if( iffld ), F.fld=F.fld/(4*pi); end

if( ifpot )
%rms_pot = norm((F.pot),2)/sqrt(nsource)
%rms_error_pot = norm((U.pot - F.pot),2)/sqrt(nsource)
rel_error_pot = norm((U.pot - F.pot),2)/norm((F.pot),2);
end

if( iffld )
%rms_fld = norm(F.fld,2)/sqrt(nsource)
%rms_error_fld = norm(U.fld - F.fld,2)/sqrt(nsource)
rel_error_fld = norm(U.fld - F.fld,2)/norm(F.fld,2);
end
%%%break;

if( ifpottarg ), U.pottarg=U.pottarg/(4*pi); end
if( iffldtarg ), U.fldtarg=U.fldtarg/(4*pi); end

if( ifpottarg ), F.pottarg=F.pottarg/(4*pi); end
if( iffldtarg ), F.fldtarg=F.fldtarg/(4*pi); end

if( ifpottarg )
%rms_pottarg = norm((F.pottarg),2)/sqrt(nsource)
%rms_error_pottarg = norm((U.pottarg - F.pottarg),2)/sqrt(ntarget)
%norm_pottarg = norm((F.pottarg),2)
rel_error_pottarg = norm((U.pottarg - F.pottarg),2)/norm((F.pottarg),2);
end

if( iffldtarg )
%rms_fldtarg = norm(F.fldtarg,2)/sqrt(ntarget)
%rms_error_fldtarg = ...
%    norm(U.fldtarg - F.fldtarg,2)/sqrt(ntarget)
rel_error_fldtarg = ...
    norm(U.fldtarg - F.fldtarg,2)/ ...
    norm(F.fldtarg,2);
end
%%%break;

sum(imag(U.fld)~=0,'all')
max(max(imag(U.fld)))