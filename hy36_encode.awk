#! /bin/awk -f
#
#  jiffy script for re-numbering and filtering large PDB files in differnt conventions
#
#
BEGIN{
}

! /^ATOM|HETAT/{print;next}


# parse the atom line
{
    if(debug) print "DEBUG" substr($0,6);
    pre=substr($0,1,22);
    post=substr($0,31);
    resnum=substr($0,23,8);
    if(resnum~/[A-Za-z]/) {
      gsub(" ","",resnum);
      resnum=hy36decode(resnum);
    }
    code=hy36encode(resnum+0);
    printf("%s%4s    %s       %d\n",pre,code,post,resnum);
}

# these are not yet fully functional
function hy36decode( code )
{
  if(! match(code,/[A-Za-z]/)) return code;
  udigits="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  ldigits=tolower(udigits);
  n=length(code);
  value=0;
  ucode=toupper(code);
  for(i=1;i<=n;++i) {
    c = substr(ucode,i,1);
    value+=(index(udigits,c)-1)*(36**(n-i));
  }
  return value - 10*36**(n-1) + 10**n;
}

function hy36encode( value )
{
  if(value < 10000) return value;

  n=4;code="";
  value = value + 10*36**3 - 10**n;
  udigits="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  ldigits=tolower(udigits);
  for(i=1;i<=n;++i) {
    q = int( value / 36**(n-i) )
    c = substr(udigits,q+1,1);
    code = code c;
    value = value - q*36**(n-i)
  }
  return code;
}

