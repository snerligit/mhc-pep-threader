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
