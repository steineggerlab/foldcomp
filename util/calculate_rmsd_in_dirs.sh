dir1=$1
dir2=$2
foldcomp_path=/Users/hbk/Projects/Lab/01_FoldU/foldcomp/build/foldcomp
for file in $dir1/*.pdb
do
    filename=$(basename $file)
    $foldcomp_path rmsd $file $dir2/$filename
done
