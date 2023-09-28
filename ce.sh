S=/cygdrive/c/Python3/python.exe
SR=/cygdrive/c/Desktop/R/R-3.4.4/bin/Rscript.exe
P=C:/Desktop/juicebox/dump/Anopheles/kr.oe/25000
O=C:/Desktop/Arbeit/R/ABCD/HC
T=C:/Desktop/Arbeit/R/ABCD/GC
R=25000
F1="kr.oe.5.ce.prc.prs"
F2="genedensity.bedGraph"

A=AatrHC
B=X
D="31 339 536 807"
C="1 4 14 30"
L1=0 
L2=16400000
#$S contrast_enhancing.py -i $P/$A.$B.$B.$R.kr.oe -o $O -l $L -r $R -d $D -c $C
L1=0 
L2=16400000
$SR eig_framed.r $O/$A.$B.$B.$R.$F1 $T/$A.$B.$R.$F2 $R 15000000 $B $L1 $L2
B=2R
D="115 914 1055 1234 1354 1682"
C="1 4 6 11 20 30"
L="0 40350000"
#$S contrast_enhancing.py -i $P/$A.$B.$B.$R.kr.oe -o $O -l $L -r $R -d $D -c $C
L1=0
L2=40350000
$SR eig_framed.r $O/$A.$B.$B.$R.$F1 $T/$A.$B.$R.$F2 $R 15000000 $B $L1 $L2


B=2L
D="127 733 1118 1493"
C="1 5 14 24"
L="6625000 45750000"
#$S contrast_enhancing.py -i $P/$A.$B.$B.$R.kr.oe -o $O -l $L -r $R -d $D -c $C
L1=6625000
L2=42750000
$SR eig_framed.r $O/$A.$B.$B.$R.$F1 $T/$A.$B.$R.$F2 $R 15000000 $B $L1 $L2

B=3R
D="146 1016"
C="1 5"
L="0 0"
#$S contrast_enhancing.py -i $P/$A.$B.$B.$R.kr.oe -o $O -l $L -r $R -d $D -c $C
L1=0
L2=54100000
$SR eig_framed.r $O/$A.$B.$B.$R.$F1 $T/$A.$B.$R.$F2 $R 15000000 $B $L1 $L2

B=3L
D="77 647 1208"
C="1 4 10"
L="6350000 41800000"
#$S contrast_enhancing.py -i $P/$A.$B.$B.$R.kr.oe -o $O -l $L -r $R -d $D -c $C
L1=6350000
L2=41800000
$SR eig_framed.r $O/$A.$B.$B.$R.$F1 $T/$A.$B.$R.$F2 $R 15000000 $B $L1 $L2

A=AfreeHC
B=X
D="31 366"
C="1 3"
L="0 0"
#$S contrast_enhancing.py -i $P/$A.$B.$B.$R.kr.oe -o $O -l $L -r $R -d $D -c $C
L1=0
L2=17250000
$SR eig_framed.r $O/$A.$B.$B.$R.$F1 $T/$A.$B.$R.$F2 $R 15000000 $B $L1 $L2

B=2R
D="49 532 1316"
C="1 3 7"
L="0 0"
#$S contrast_enhancing.py -i $P/$A.$B.$B.$R.kr.oe -o $O -l $L -r $R -d $D -c $C
L1=0
L2=54475000
$SR eig_framed.r $O/$A.$B.$B.$R.$F1 $T/$A.$B.$R.$F2 $R 15000000 $B $L1 $L2

B=2L
D="47 1190"
C="2 8"
L="0 0"
#$S contrast_enhancing.py -i $P/$A.$B.$B.$R.kr.oe -o $O -l $L -r $R -d $D -c $C
L1=0
L2=35275000
$SR eig_framed.r $O/$A.$B.$B.$R.$F1 $T/$A.$B.$R.$F2 $R 15000000 $B $L1 $L2

B=3R
D="52 546 1206"
C="1 3 7"
L="0 0"
#$S contrast_enhancing.py -i $P/$A.$B.$B.$R.kr.oe -o $O -l $L -r $R -d $D -c $C
L1=0
L2=41525000
$SR eig_framed.r $O/$A.$B.$B.$R.$F1 $T/$A.$B.$R.$F2 $R 15000000 $B $L1 $L2

B=3L
D="48 439 1050"
C="1 3 7"
L="3600000 42125000"
#$S contrast_enhancing.py -i $P/$A.$B.$B.$R.kr.oe -o $O -l $L -r $R -d $D -c $C
L1=3600000
L2=42125000
$SR eig_framed.r $O/$A.$B.$B.$R.$F1 $T/$A.$B.$R.$F2 $R 15000000 $B $L1 $L2