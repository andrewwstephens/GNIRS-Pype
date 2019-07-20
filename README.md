# GNIRS_PIPELINE_PLAN


————————————————————————————— 
JUL 5 - AUG 20  ::  GNIRS DATA REDUCTION PIPELINE PLAN
—————————————————————————————

JUL 05    ::    gnirsPipeline.py    YES

JUL 08 - JUL 10    ::    gnirsPipeline.py    YES
                                    gnirsSort.py    NO
                                    gnirsSort.py Test Run*    YES

JUL 10    ::    Revise plan (if required)    YES

JUL 11 - JUL 12    ::    gnirsSort.py    NO
                                    gnirsSort.py Test Run*    YES
                                    gnirsBaselineCalibration.py    NO

**JUL 17 - JUL 19    ::    gnirsBaselineCalibration.py Test Run    NO

JUL 19    ::    Revise plan (if required)    YES

JUL 22 - JUL 26    ::    gnirsBaselineCalibration.py    ____
                                    gnirsReduce.py    ____

JUL 26    ::    Revise plan (if required)    ____

JUL 29  - AUG 02    ::    gnirsTelluric.py    ____
                                       gnirsFluxCalibrate.py    ____

AUG 05 - AUG 09    ::    gnirsCombineOrdersXD.py    ____
                                       gnirsCalcSNR.py    ____
                                       writeDataSheet.py    ____

AUG 07    ::    Revise plan (if required)    ____

AUG 12 - AUG 15***    ::    gnirsSort.py    ____
                                           remaining documentation    ____

AUG 19 - AUG 20    ::    Miscellaneous (wrap up) __


*  Add gnirsReadHeaders.py as required while working on gnirsSort.py    [YES]

**  JUL 15 - JUL 16:  O’ahu Visit    [YES]

***  AUG 16:  Holiday
