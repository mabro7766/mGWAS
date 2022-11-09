#!/bin/bash

echo "adding filename ..."

for file in output_trait*.txt
   do sed -i "s/$/\t$file/" $file
done

head -n 1 output_trait1.txt > output.txt && tail -q -n +2 output_trait*.txt >> output.txt

## note: the following code works only when all headers are the same
## due to the sed loop, each header also contains filename, resulting in different headers

# awk '
#     FNR==1 && NR!=1 { while (/^<header>/) getline; }
#     1 {print}
# ' output_trait*.txt > all.txt