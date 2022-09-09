function [e,se,A,B]=sampen(y,m,r,sflag,cflag,vflag)
%function e=sampen(y,m,r);
%
%Input Parameters
%
%y  input signal vector
%m  maximum template length (default m=5)
%r  matching threshold (default r=.2)
%
%Output Parameters
%
%e sample entropy estimates for m=0,1,...,m-1
%
%Full usage:
%
%[e,se,A,B]=sampen(y,m,r,sflag,cflag,vflag)
%
%Input Parameters
%
%sflag    flag to standardize signal(default yes/sflag=1) 
%cflag    flag to use fast C code (default yes/cflag=1) 
%vflag    flag to calculate standard errors (default no/vflag=0) 
%
%Output Parameters
%
%se standard error estimates for m=0,1,...,m-1
%A number of matches for m=1,...,m
%B number of matches for m=0,...,m-1
%  (excluding last point in Matlab version)

if ~exist('m')|isempty(m),m=5;end
if ~exist('r')|isempty(r),r=.2;end
if ~exist('sflag')|isempty(sflag),sflag=1;end
if ~exist('cflag')|isempty(cflag),cflag=1;end
if ~exist('vflag')|isempty(cflag),vflag=0;end
y=y(:);
n=length(y);
if sflag>0
   y=y-mean(y);
   s=sqrt(mean(y.^2));   
   y=y/s;
end
if nargout>1
    if vflag>0
        se=sampense(y,m,r);
    else
        se=[];
    end
end    
if cflag>0
   [match,R]=cmatches(y,n,r);
   match=double(match);
else   
   [e,A,B]=sampenc(y,m,r);
   return
end
k=length(match);
if k<m
   match((k+1):m)=0;
end
N=n*(n-1)/2;
A=match(1:m);
B=[N;A(1:(m-1))];
N=n*(n-1)/2;
p=A./B;
e=-log(p);
