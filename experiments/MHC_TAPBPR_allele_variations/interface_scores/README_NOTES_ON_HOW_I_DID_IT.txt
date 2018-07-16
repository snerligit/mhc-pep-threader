/home/snerli/rosetta/Rosetta_ref2015/main/source/bin/InterfaceAnalyzer.linuxgccrelease -database /home/snerli/rosetta/Rosetta_ref2015/main/database -in:file:silent *.out -cutpoint 378 -out:file:score_only interface_score.sc

where 387 is the residue that separates the MHC/beta2m from TAPBPR
where *.out is each of the “silent” files generated from the Rosetta relaxes

The interface_score.sc file will store the interface scores, which can be extracted using awk.
You might be interested in "dG_separated score" column which is in column 7.
To extract it use the command:

awk '{print $7,$NF}' interface_score.sc > results_interface_score.txt

To calculate the average of each 10th row (because there are ten structures for each) use script awk.sh

#!/bin/bash

awk -v count=10 '
    {
        if ( NF > tot_col )
            tot_col = NF;

        cur = 1;
        while ( cur <= NF )
        {
            sums[cur] += $cur;
            cur++;
        }

        if ( ( NR % count ) == 0 )
        {
            cur = 1;
            while ( cur <= tot_col )
            {
                printf("%0.2f ", sums[cur] / count);
                cur++;
            }

            print "";

            delete sums;
            tot_col = 0;
        }
    }' "$@"

Then use the following commands:

chmod 777 awk.sh

./awk.sh results_interface_score.txt > test1.txt

awk 'NR == 1 || NR % 10 == 0' results_interface_score.txt | awk '{print $2}' > test2.txt

paste test1.txt test2.txt > combined_average_score.txt

awk '{print $1, $3}' combined_average_score.txt > combined_average_score_final.txt

rm combined_average_score.txt






