#!/bash/bin                                                                       

# drawn signal samples                                                             
SIGNAL=(SS_onestepCC)
#SIGNAL=(SS_direct GG_direct GG_onestepCC SS_onestepCC SM_GG_N2)
# luminosity (please change it by your hand here)                                 
LUMI=36.1ifb
# draw observed limit or not                                                      
UNBLIND=1
# draw each SR's contours                                                         
ALLSR=0

OUTPUTDIR=Outputs
PLOTDIR=plots

mkdir -v $OUTPUTDIR
mkdir -v $OUTPUTDIR/${LUMI}

#:<<'#_COMMENT_'                                                                  
#for TYPE in ${SIGNAL[@]}
#do
  # run all processes (merging, make contours, choose the best SR(Oring))         
#  ./makeContours_Run2.py --all --grid ${TYPE} --outputDir "$OUTPUTDIR/$LUMI/$TYPE" 2>&1 | tee log-MakeContours_Run2_${TYPE}.out
  # making contours & choosing the best SR (Oring)                                
#  ./makeContours_Run2.py -c -o --grid ${TYPE} --outputDir "$OUTPUTDIR/$LUMI/$TYPE" 2>&1 | tee log-MakeContours_Run2_${TYPE}.out                                    
  # choosing the best SR (Oring)                                                  
  #./makeContours_Run2.py -o --grid ${TYPE} --outputDir "$OUTPUTDIR/$LUMI/$TYPE" 2>&1 | tee log-MakeContours_Run2_${TYPE}.out                                       
#done
#_COMMENT_                                                                        

#:<<'#_COMMENT_'                                                                  
mkdir -v $PLOTDIR
for TYPE in ${SIGNAL[@]}
do
 # cp -v $OUTPUTDIR/$LUMI/$TYPE/summary_harvest_tree_description.h summary_harvest_tree_description.h
  root -l -q -b "makecontourplots_CLs.C(\"$TYPE\",\"$OUTPUTDIR/$LUMI/$TYPE\",\"$PLOTDIR\",0,$UNBLIND,$ALLSR)" 2>&1 | tee log-makecontourplots_CLs_${TYPE}.out
  root -l -q -b "makecontourplots_CLs.C(\"$TYPE\",\"$OUTPUTDIR/$LUMI/$TYPE\",\"$PLOTDIR\",1,$UNBLIND,$ALLSR)" 2>&1 | tee log-makecontourplots_CLs_${TYPE}.out
done
#_COMMENT_                                                                        