#!/bin/zsh

taxa="AA1 AA2 AL1 AL2 AH1 AH2 AH3"
scale="within between all"
depth="05 10 20"

for name in $taxa
    for item in $scale
        if [ $taxa == "AA1" ] || [ $taxa == "AL1" ] || [ $taxa == "AL2" ]
            then
                depth=20
                do Rscript distance_plots.R $name $depth $scale
        elif  [ $taxa == "AH2" ] || [ $taxa == "AH3" ]
            then
                depth=05
                do Rscript distance_plots.R $name $depth $scale
        else
            for number in $depth
                then
                    do Rscript distance_plots.R $name $depth $scale
        fi
done done

