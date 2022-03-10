
for date in "20181110" "20181124" "20190106" "20190130"
do
    for exp in "default" "colMix2_Akernel"
    do
        echo $exp
        python prec_field.py --date $date -int 0 -exp $exp
        #convert -delay 20 -loop 0 /home/mkarrer/Dokumente/plots/ICON/precip_field/${date}/*${exp}* /home/mkarrer/Dokumente/plotsICON/precip_field/${date}/${date}_obs.gif
    done
done
