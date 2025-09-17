"""FitParameters for glimepiride fitting."""

from sbmlsim.fit import FitParameter


parameters_all = [
    # tissue distribution
    FitParameter(
        pid="ftissue_gli",
        lower_bound=1E-4,
        start_value=0.1,
        upper_bound=1,
        unit="l/min",
    ),
    FitParameter(
        pid="Kp_gli",
        lower_bound=10,
        start_value=660,
        upper_bound=1000,
        unit="dimensionless",
    ),

    # absorption rate
    FitParameter(
        pid="GU__GLIABS_k",
        lower_bound=1E-4,
        start_value=0.02,
        upper_bound=1,
        unit="1/min",
    ),
    # reabsorption rate
    FitParameter(
        pid="GU__MREABS_k",
        lower_bound=1E-4,
        start_value=0.05,
        upper_bound=1,
        unit="1/min",
    ),
    # feces excretion
    FitParameter(
        pid="GU__MEXC_k",
        lower_bound=1E-4,
        start_value=0.1,
        upper_bound=1,
        unit="1/min",
    ),


    # hepatic metabolism
    # FIXME: import set fast
    # FitParameter(
    #     pid="LI__GLIIM_k",
    #     lower_bound=1E-3,
    #     start_value=1.0,
    #     upper_bound=100,
    #     unit="1/min",
    # ),
    FitParameter(
        pid="LI__GLI2M1_Vmax",
        lower_bound=1E-7,
        start_value=0.0001,
        upper_bound=0.01,
        unit="mmole/min/l",
    ),
    FitParameter(
        pid="LI__M1EX_k",
        lower_bound=1E-3,
        start_value=1.0,
        upper_bound=100,
        unit="1/min",
    ),
    FitParameter(
        pid="LI__M12M2_k",
        lower_bound=1E-3,
        start_value=1.0,
        upper_bound=100,
        unit="1/min",
    ),

    FitParameter(
        pid="LI__M2EX_k",
        lower_bound=1E-3,
        start_value=1.0,
        upper_bound=100,
        unit="1/min",
    ),

    # kidney removal
    FitParameter(
        pid="KI__M1EX_k",
        lower_bound=1E-4,
        start_value=1E-1,
        upper_bound=1,
        unit="1/min",
    ),
    FitParameter(
        pid="KI__M2EX_k",
        lower_bound=1E-4,
        start_value=1E-1,
        upper_bound=1,
        unit="1/min",
    ),
]

