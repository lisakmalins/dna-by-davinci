#SUFFIX=$(echo $1 | grep -o '[^.]*$')
SUFFIX=$(echo $1 | grep -o '\..*')
BASE=$(basename $1 $SUFFIX)
echo $1
echo $SUFFIX
echo $BASE
echo $BASE'-1'$SUFFIX
cat $1 \
    | paste - - - - - - - - \
    | tee >(cut -f 1-4 | tr '\t' '\n' | gzip > $BASE'-1'$SUFFIX) \
    | cut -f 5-8 | tr '\t' '\n' | gzip > $BASE'-2'$SUFFIX 
