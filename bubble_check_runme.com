#! /bin/tcsh -f
#
#   Look for vacuum bubbles big enough to hold water molecules in a pdb file      -James Holton 1-5-24
#
#
set pdbfile = ""

set Vwater = 96
set minvoid = 5

set quiet = 0
set debug = 0


set tempfile = /dev/shm/${USER}/tempfile_void_$$_

foreach Arg ( $* )
    set arg = `echo $Arg | awk '{print tolower($0)}'`
    set assign = `echo $arg | awk '{print ( /=/ )}'`
    set Key = `echo $Arg | awk -F "=" '{print $1}'`
    set Val = `echo $Arg | awk '{print substr($0,index($0,"=")+1)}'`
    set Csv = `echo $Val | awk 'BEGIN{RS=","} {print}'`
    set key = `echo $Key | awk '{print tolower($1)}'`
    set num = `echo $Val | awk '{print $1+0}'`
    set int = `echo $Val | awk '{print int($1+0)}'`

    if( $assign ) then
      # re-set any existing variables
      set test = `set | awk -F "\t" '{print $1}' | egrep "^${Key}"'$' | wc -l`
      if ( $test ) then
          set $Key = $Val
          echo "$Key = $Val"
          continue
      endif
      # synonyms
      if("$key" == "output") set outfile = "$Val"
    else
      # no equal sign
      if("$Arg" =~ *.pdb ) then
         set pdbfile = $Arg
      endif
      if("$Arg" =~ *.rst7 ) then
         set rstfile = $Arg
      endif
    endif
    if("$key" == "debug") set debug = "$Val"
end

if( $debug && $tempfile =~ /dev/shm/* ) set tempfile = ./tempfile_debug_

if(! -e "$pdbfile" ) then
    set BAD = "no coordinates provided"
    goto exit
endif


set t = ${tempfile}



gemmi mask ${pdbfile} ${t}void.msk
gemmi map2sf ${t}void.msk ${t}void.mtz Fbulk PHIbulk
gemmi blobs ${t}void.mtz -f Fbulk -p PHIbulk $pdbfile  >! ${t}gemmi_blobs.txt
# voids only count if more than $minvoid waters fit
awk -v Vwater=$Vwater '$3=="el"{print $5/Vwater}'  ${t}gemmi_blobs.txt |\
awk '$1>1' >!  ${t}void_list.txt
set nbulk_sum = `awk -v minvoid=$minvoid '$1>=minvoid{sum+=$1} END{print sum+0}' ${t}void_list.txt`
set nbulk_max = `awk 'NR==1{print $1}' ${t}void_list.txt`
if("$nbulk_max" == "") set nbulk_max = 0


echo "$nbulk_max waters can fit in biggest void"
echo "$nbulk_sum total waters can fit in all voids bigger than $minvoid waters"

set test = `echo $nbulk_max $minvoid | awk '{print ( $1 < $2/2 )}'`
if( $test ) then
   cat << EOF
recommend dropping some waters. Try:

remap_waters_runme.com amber.rst7 diffmap.mtz mtzlabel=FOFCWT
dehydrate_amber_runme.com xtal.prmtop remapped.rst7 diffmap.mtz mtzlabel=FOFCWT maxreject=10

EOF
endif
  set test = `echo $nbulk_max $minvoid | awk '{print ( $1 > 2*$2 )}'`
  if( $test ) then
    set nadd = `echo $nbulk_max | awk '{print int($1)}'`
    cat << EOF
recommend adding $nadd waters

hydrate_runme.com amber.rst7 outrst=wetter.rst7 paddedparm=padded.parm7 nadd=$nadd
endif


exit:

if( $?BAD ) then
   echo "ERROR: $BAD"
   exit 9   
endif

if( $debug ) exit
if( "$t" != "" && "$t" != "./") then
   rm -f ${t}* > /dev/null
endif

exit



