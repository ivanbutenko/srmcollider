for file in $*
do
    epstopdf ${file}
    rm ${file}
done

