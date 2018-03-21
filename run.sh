CALIB="/Users/alexandertuna/Downloads/CALIBRATIONS/2017-10-26-hybrid"
PDO="${CALIB}/hybrid_PDO.root"
TDO="${CALIB}/AllBoards_2017-10-26_Calib_TDO.root"
# ALI="alignment_3545-3548_tranXZ.root"

RUN="3545-3549"
DATA="/Users/alexandertuna/Downloads/octuplet_analysis_3545/Run${RUN}.root"

make TPCAnalysis.x
./TPCAnalysis.x -i ${DATA} -o test_${RUN}.root -p ${PDO} -t ${TDO} -g # -a ${ALI}
python plots_for_talk.py -r 3545-3549

say caw
