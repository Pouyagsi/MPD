
import pandas as pd
import math
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


# Read rheological data
PV = pd.read_excel("Data.xlsx", sheet_name='Mud Rheological Data', usecols="B", nrows=1).iloc[0, 0]
YV = pd.read_excel("Data.xlsx", sheet_name='Mud Rheological Data', usecols="B", nrows=2, skiprows=1).iloc[0, 0]
Theta_600 = pd.read_excel("Data.xlsx", sheet_name='Mud Rheological Data', usecols="B", nrows=3, skiprows=2).iloc[0, 0]
Theta_300 = pd.read_excel("Data.xlsx", sheet_name='Mud Rheological Data', usecols="B", nrows=4, skiprows=3).iloc[0, 0]
Theta_6 = pd.read_excel("Data.xlsx", sheet_name='Mud Rheological Data', usecols="B", nrows=5, skiprows=4).iloc[0, 0]
Theta_3 = pd.read_excel("Data.xlsx", sheet_name='Mud Rheological Data', usecols="B", nrows=6, skiprows=5).iloc[0, 0]
n = math.log(Theta_600 / Theta_300) / math.log(600 / 300)
K = 511 * Theta_300 / (511 ** n)


# Read pump data
Q = pd.read_excel("Data.xlsx", sheet_name='Pump', usecols="B", nrows=1).iloc[0, 0]
MW = pd.read_excel("Data.xlsx", sheet_name='Pump', usecols="B", nrows=2, skiprows=1).iloc[0, 0]

# Read drilling parameters
Hole_Diameter = pd.read_excel("Data.xlsx", sheet_name='Drilling Parameters', usecols="B", nrows=1).iloc[0, 0]
OD_Dc = pd.read_excel("Data.xlsx", sheet_name='Drilling Parameters', usecols="B", nrows=2, skiprows=1).iloc[0, 0]
OD_Dp = pd.read_excel("Data.xlsx", sheet_name='Drilling Parameters', usecols="B", nrows=3, skiprows=2).iloc[0, 0]
ID_Casing = pd.read_excel("Data.xlsx", sheet_name='Drilling Parameters', usecols="B", nrows=4, skiprows=3).iloc[0, 0]
Mud_Window = pd.read_excel("Data.xlsx", sheet_name='Required Mud Window')
try:
    # Read additional drilling parameters
    Dp_Length = pd.read_excel("Data.xlsx", sheet_name='Drilling Parameters', usecols="B", nrows=5, skiprows=4).iloc[0, 0]
    DC_Length = pd.read_excel("Data.xlsx", sheet_name='Drilling Parameters', usecols="B", nrows=6, skiprows=5).iloc[0, 0]
    Casing_shoe_Depth = pd.read_excel("Data.xlsx", sheet_name='Drilling Parameters', usecols="B", nrows=7, skiprows=6).iloc[0, 0]
    Start_depth = pd.read_excel("Data.xlsx", sheet_name='Drilling Parameters', usecols="B", nrows=9, skiprows=8).iloc[0, 0]
    End_depth = pd.read_excel("Data.xlsx", sheet_name='Drilling Parameters', usecols="B", nrows=10, skiprows=9).iloc[0, 0]
    Deviation_start_depth = pd.read_excel("Data.xlsx", sheet_name='Drilling Parameters', usecols="B", nrows=11, skiprows=10).iloc[0, 0]
    Deviation_end_depth = pd.read_excel("Data.xlsx", sheet_name='Drilling Parameters', usecols="B", nrows=12, skiprows=11).iloc[0, 0]
    Horizontal_start_depth = pd.read_excel("Data.xlsx", sheet_name='Drilling Parameters', usecols="B", nrows=13, skiprows=12).iloc[0, 0]
    Horizontal_length = pd.read_excel("Data.xlsx", sheet_name='Drilling Parameters', usecols="B", nrows=14, skiprows=13).iloc[0, 0]
    Incl = pd.read_excel("Data.xlsx", sheet_name='Survey Data', usecols="B", nrows=3, skiprows=2).iloc[0, 0]
    
    # Converting Degree to Radian
    DTR = Incl * (math.pi)/180   
    # Calculate deviated and horizontal sections
    Deviation_length = (Deviation_end_depth - Deviation_start_depth) / (math.cos(DTR))
    
    # Calculate total depth
    MD = Horizontal_length + Deviation_length + Deviation_start_depth 
    
    # Calculate Drilled_depth
    Step = 1
    Drilled_depth = [i * Step + Start_depth for i in range((int(MD) - int(Start_depth)) // Step + 1)]
    
    # Calculate Dp_Length
    Dp_Length = [i * Step + Dp_Length for i in range((int(MD) - int(Start_depth)) // Step + 1)]
    
    L1_list = []
    L2_list = []
    L3_list = []
    L4_list = []

    # Initialize lists for Power Law Model
    P1_list = []
    P2_list = []
    P3_list = []
    P4_list = []
    AF_PL_list = []
    
    # Initialize lists for Herschel-Bulkley Model
    P1_HB_list = []
    P2_HB_list = []
    P3_HB_list = []
    P4_HB_list = []
    AF_HB_list = []

    # Initialize section identifiers
    Section = []

    for depth in Drilled_depth:
        
        if depth < Deviation_start_depth:
           section = 'Vertical'
           
        elif Deviation_start_depth <= depth < (Deviation_start_depth + (Deviation_end_depth-Deviation_start_depth)/(math.cos(DTR))) :
             section = 'Deviated'
            
        elif (Deviation_start_depth + (Deviation_end_depth-Deviation_start_depth)/(math.cos(DTR))) <= depth:
             section = 'Horizontal'
        Section.append(section)  

        
    for depth, dp_length in zip(Drilled_depth, Dp_Length):
        if depth > Casing_shoe_Depth:
           L1 = min(DC_Length, depth - Casing_shoe_Depth)
           L2 = max(0, depth - Casing_shoe_Depth - L1)
           L3 = max(0, DC_Length - L1)
           L4 = max(0, dp_length - L2)
        else:
           L1, L2, L3, L4 = 0, 0, DC_Length, dp_length
                 
        L1_list.append(L1)
        L2_list.append(L2)
        L3_list.append(L3)
        L4_list.append(L4)
        
 
        
    # Calculate V_Annulus1
    V_Annulus1 = (1.28342246 * Q) / (math.pi * ((Hole_Diameter) ** 2 - (OD_Dc) ** 2))
        
    # Calculate Re_PL1
    Re_PL1 = 109000 * ((MW * V_Annulus1 ** (2 - n)) / (K)) * ((0.0208 * (Hole_Diameter - OD_Dc)) / (2 + (1 / n))) ** n
        
    # Calculate V_Annulus2
    V_Annulus2 = (1.28342246 * Q) / (math.pi * ((Hole_Diameter) ** 2 - (OD_Dp) ** 2))
        
    # Calculate Re_PL2
    Re_PL2 = 109000 * ((MW * V_Annulus2 ** (2 - n)) / (K)) * ((0.0208 * (Hole_Diameter - OD_Dp)) / (2 + (1 / n))) ** n
    
    # Calculate V_Annulus3
    V_Annulus3 = (1.28342246 * Q) / (math.pi * ((ID_Casing) ** 2 - (OD_Dc) ** 2))
        
    # Calculate Re_PL3
    Re_PL3 = 109000 * ((MW * V_Annulus3 ** (2 - n)) / (K)) * ((0.0208 * (ID_Casing - OD_Dc)) / (2 + (1 / n))) ** n
        
    # Calculate V_Annulus4
    V_Annulus4 = (1.28342246 * Q) / (math.pi * ((ID_Casing) ** 2 - (OD_Dp) ** 2))
        
    # Calculate Re_PL4
    Re_PL4 = 109000 * ((MW * V_Annulus4 ** (2 - n)) / (K)) * ((0.0208 * (ID_Casing - OD_Dp)) / (2 + (1 / n))) ** n

    for i in range(len(L1_list)):

    # Calculate Annular Friction loss (Power Law Model)
       if Re_PL1 <= 3470 - 1370 * n:
          Laminar_Delta_P1 = (((144 * V_Annulus1 * ((2 * n) + 1)) / (((Hole_Diameter - OD_Dc) * 3 * n))) ** n) * (
                    K * L1_list[i] / (300 * (Hole_Diameter - OD_Dc)))
          Turbolant_Delta_P1 = 0
          
       elif Re_PL1 >= 4270 - 1370 * n:
            f1 = (0.0791 / (Re_PL1 ** 0.25))
            Turbolant_Delta_P1 = (f1 * MW * (V_Annulus1 ** 2) * L1_list[i]) / (21.1 * (Hole_Diameter - OD_Dc))
            Laminar_Delta_P1 = 0
        
       # Append P1 to the list
       P1 = Laminar_Delta_P1 + Turbolant_Delta_P1
       P1_list.append(P1)

       if Re_PL2 <= 3470 - 1370 * n:
          Laminar_Delta_P2 = (((144 * V_Annulus2 * ((2 * n) + 1)) / (((Hole_Diameter - OD_Dp) * 3 * n))) ** n) * (
                    K * L2_list[i] / (300 * (Hole_Diameter - OD_Dp)))
          Turbolant_Delta_P2 = 0
          
       elif Re_PL2 >= 4270 - 1370 * n:
            f2 = (0.0791 / (Re_PL2 ** 0.25))
            Turbolant_Delta_P2 = (f2 * MW * (V_Annulus2 ** 2) * L2_list[i]) / (21.1 * (Hole_Diameter - OD_Dp))
            Laminar_Delta_P2 = 0
        
       # Append P2 to the list
       P2 = Laminar_Delta_P2 + Turbolant_Delta_P2
       P2_list.append(P2)

       if Re_PL3 <= 3470 - 1370 * n:
          Laminar_Delta_P3 = (((144 * V_Annulus3 * ((2 * n) + 1)) / (((ID_Casing - OD_Dc) * 3 * n))) ** n) * (
                    K * L3_list[i] / (300 * (ID_Casing - OD_Dc)))
          Turbolant_Delta_P3 = 0
          
       elif Re_PL3 >= 4270 - 1370 * n:
            f3 = (0.0791 / (Re_PL3 ** 0.25))
            Turbolant_Delta_P3 = (f3 * MW * (V_Annulus3 ** 2) * L3_list[i]) / (21.1 * (ID_Casing - OD_Dc))
            Laminar_Delta_P3 = 0
        
       # Append P3 to the list
       P3 = Laminar_Delta_P3 + Turbolant_Delta_P3
       P3_list.append(P3)

       if Re_PL4 <= 3470 - 1370 * n:
          Laminar_Delta_P4 = (((144 * V_Annulus4 * ((2 * n) + 1)) / (((ID_Casing - OD_Dp) * 3 * n))) ** n) * (
                    K * L4_list[i] / (300 * (ID_Casing - OD_Dp)))
          Turbolant_Delta_P4 = 0
          
       elif Re_PL4 >= 4270 - 1370 * n:
            f4 = (0.0791 / (Re_PL4 ** 0.25))
            Turbolant_Delta_P4 = (f4 * MW * (V_Annulus4 ** 2) * L4_list[i]) / (21.1 * (ID_Casing - OD_Dp))
            Laminar_Delta_P4 = 0
        
       # Append P4 to the list
       P4 = Laminar_Delta_P4 + Turbolant_Delta_P4
       P4_list.append(P4)



  # Calculate AF_PL for each iteration
    AF_PL_list = [P1 + P2 + P3 + P4 for P1, P2, P3, P4 in zip(P1_list, P2_list, P3_list, P4_list)]

  # Annular Friction loss (Herschel-Bulkley Model)
 
    YV_HB = 0.0020885 * (2*Theta_3 - Theta_6) 
    MW_HB = 7.48 * MW
    q = 0.0022280093 * Q
    n_HB = math.log((Theta_600 - YV_HB)/ (Theta_300 - YV_HB)) / math.log(600 / 300)
    K_HB = 0.0020885 * (500 * (Theta_300 - YV) / (511**n_HB))
    
    if Q == 0:
       AF_HB_list  = AF_PL_list
       
    else:    
         Ca_1 = 1 - ((1/(1 + n_HB))*(YV_HB/(YV_HB+ K_HB*((2*q*(2*n+1))/((n_HB * math.pi)*((Hole_Diameter-OD_Dc)/24)*(((Hole_Diameter**2)-(OD_Dc**2))/48)))**n_HB)))
         Re_HB1 = (4*(2*n_HB+1)/(n_HB))*((MW_HB * (V_Annulus1 **(2-n_HB))*((Hole_Diameter-OD_Dc)/24) ** (n_HB))/(YV_HB*((Hole_Diameter-OD_Dc)/(24* V_Annulus1))**(n_HB) + K_HB *((4*n_HB + 2)/(n_HB * Ca_1))**(n_HB)))

         Ca_2 = 1 - ((1/(1 + n_HB))*(YV_HB/(YV_HB+ K_HB*((2*q*(2*n+1))/((n_HB * math.pi)*((Hole_Diameter-OD_Dp)/24)*(((Hole_Diameter**2)-(OD_Dp**2))/48)))**n_HB)))
         Re_HB2 = (4*(2*n_HB+1)/(n_HB))*((MW_HB * (V_Annulus2 **(2-n_HB))*((Hole_Diameter-OD_Dp)/24) ** (n_HB))/(YV_HB*((Hole_Diameter-OD_Dp)/(24* V_Annulus2))**(n_HB) + K_HB *((4*n_HB + 2)/(n_HB * Ca_2))**(n_HB)))

         Ca_3 = 1 - ((1/(1 + n_HB))*(YV_HB/(YV_HB+ K_HB*((2*q*(2*n+1))/((n_HB * math.pi)*((ID_Casing-OD_Dc)/24)*(((ID_Casing**2)-(OD_Dc**2))/48)))**n_HB)))
         Re_HB3 = (4*(2*n_HB+1)/(n_HB))*((MW_HB * (V_Annulus3 **(2-n_HB))*((ID_Casing-OD_Dc)/24) ** (n_HB))/(YV_HB*((ID_Casing-OD_Dc)/(24* V_Annulus3))**(n_HB) + K_HB *((4*n_HB + 2)/(n_HB * Ca_3))**(n_HB)))

         Ca_4 = 1 - ((1/(1 + n_HB))*(YV_HB/(YV_HB+ K_HB*((2*q*(2*n+1))/((n_HB * math.pi)*((ID_Casing-OD_Dp)/24)*(((ID_Casing**2)-(OD_Dp**2))/48)))**n_HB)))
         Re_HB4 = (4*(2*n_HB+1)/(n_HB))*((MW_HB * (V_Annulus4 **(2-n_HB))*((ID_Casing-OD_Dp)/24) ** (n_HB))/(YV_HB*((ID_Casing-OD_Dp)/(24* V_Annulus4))**(n_HB) + K_HB *((4*n_HB + 2)/(n_HB * Ca_4))**(n_HB)))

         Y = (math.log(n_HB) + 3.93 )/50
         Z = (1.75 - math.log(n_HB))/7
         Re_HBc = ((8*(2*n+1)/(n_HB*Y)))**(1/(1-Z))
    

         for i in range(len(L1_list)):
        
              if Re_HB1 < Re_HBc:
                 Laminar_Delta_P1_HB = ((4*K_HB*L1_list[i])/(14400*(Hole_Diameter-OD_Dc)/12))*((YV_HB/K_HB)+(((32*n_HB+16)/(n_HB*Ca_1*(Hole_Diameter-OD_Dc)/12))*(q/(math.pi*((Hole_Diameter**2)-(OD_Dc**2))/144)))**n_HB)
                 Turbolant_Delta_P1_HB = 0
           
              elif Re_HB1 > Re_HBc:
                   Fa1 = Y * (Ca_1 * Re_HB1)**(-Z)
                   Turbolant_Delta_P1_HB = (L1_list[i] * Fa1 * MW_HB * q**2)/(1421.22*((Hole_Diameter-OD_Dc)/12)*(((Hole_Diameter**2)-(OD_Dc**2))/144)**2)
                   Laminar_Delta_P1_HB = 0
          
              P1_HB = Laminar_Delta_P1_HB + Turbolant_Delta_P1_HB
              P1_HB_list.append(P1_HB)
 
              if Re_HB2 < Re_HBc:
                 Laminar_Delta_P2_HB = ((4*K_HB*L2_list[i])/(14400*(Hole_Diameter-OD_Dp)/12))*((YV_HB/K_HB)+(((32*n_HB+16)/(n_HB*Ca_2*(Hole_Diameter-OD_Dp)/12))*(q/(math.pi*((Hole_Diameter**2)-(OD_Dp**2))/144)))**n_HB)
                 Turbolant_Delta_P2_HB = 0
              
              elif Re_HB2 > Re_HBc :
                   Fa2 = Y * (Ca_2 * Re_HB2)**(-Z)
                   Turbolant_Delta_P2_HB = (L2_list[i] * Fa2 * MW_HB * q**2)/(1421.22*((Hole_Diameter-OD_Dp)/12)*(((Hole_Diameter**2)-(OD_Dp**2))/144)**2)
                   Laminar_Delta_P2_HB = 0
         
              P2_HB = Laminar_Delta_P2_HB + Turbolant_Delta_P2_HB
              P2_HB_list.append(P2_HB) 
              
              if Re_HB3 < Re_HBc:
                 Laminar_Delta_P3_HB = ((4*K_HB*L3_list[i])/(14400*(ID_Casing-OD_Dc)/12))*((YV_HB/K_HB)+(((32*n_HB+16)/(n_HB*Ca_3*(ID_Casing-OD_Dc)/12))*(q/(math.pi*((ID_Casing**2)-(OD_Dc**2))/144)))**n_HB)
                 Turbolant_Delta_P3_HB = 0
           
              elif Re_HB3 > Re_HBc :
                   Fa3 = Y * (Ca_3 * Re_HB3)**(-Z)
                   Turbolant_Delta_P3_HB = (L3_list[i] * Fa3 * MW_HB * q**2)/(1421.22*((ID_Casing-OD_Dc)/12)*(((ID_Casing**2)-(OD_Dc**2))/144)**2)
                   Laminar_Delta_P3_HB = 0
          
              P3_HB = Laminar_Delta_P3_HB + Turbolant_Delta_P3_HB
              P3_HB_list.append(P3_HB)    
 
              if Re_HB4 < Re_HBc:
                 Laminar_Delta_P4_HB = ((4*K_HB*L4_list[i])/(14400*(ID_Casing-OD_Dp)/12))*((YV_HB/K_HB)+(((32*n_HB+16)/(n_HB*Ca_4*(ID_Casing-OD_Dp)/12))*(q/(math.pi*((ID_Casing**2)-(OD_Dp**2))/144)))**n_HB)
                 Turbolant_Delta_P4_HB = 0
           
              elif Re_HB4 > Re_HBc :
                   Fa4 = Y * (Ca_4 * Re_HB4)**(-Z)
                   Turbolant_Delta_P4_HB = (L4_list[i] * Fa4 * MW_HB * q**2)/(1421.22*((ID_Casing-OD_Dp)/12)*(((ID_Casing**2)-(OD_Dp**2))/144)**2)
                   Laminar_Delta_P4_HB = 0
          
              P4_HB = Laminar_Delta_P4_HB + Turbolant_Delta_P4_HB
              P4_HB_list.append(P4_HB)

            # Calculate AF_HB for each iteration
         AF_HB_list = [P1_HB + P2_HB + P3_HB + P4_HB for P1_HB, P2_HB, P3_HB, P4_HB in zip(P1_HB_list, P2_HB_list, P3_HB_list, P4_HB_list)]

    
    # ECD Calculation Based on Power Law Model
    TVD = []
    for section, depth in zip(Section, Drilled_depth):
        if section == 'Vertical':
           tvd = depth
        elif section == 'Deviated':
             tvd = Deviation_start_depth + (depth - Deviation_start_depth) * math.cos(DTR)
        elif section == 'Horizontal':
             tvd = Deviation_end_depth
        TVD.append(tvd)
       
    ECD_PL = []
    for depth, af_pl in zip(TVD, AF_PL_list):
        ecd_pl = ((af_pl / (0.052 * depth)) + MW)
        ECD_PL.append(ecd_pl)

    # ECD Calculation Based on HB Model
    #TVD = Drilled_depth
    ECD_HB = []
    for depth, af_hb in zip(TVD, AF_HB_list):
        ecd_hb = ((af_hb / (0.052 * depth)) + MW)
        ECD_HB.append(ecd_hb)

    # BHP Calculation Based on Power Law Model
    BHP_PL = []
    for depth, ecd_pl in zip(TVD, ECD_PL):
        bhp_pl = ((0.052 * depth) * ecd_pl)
        BHP_PL.append(bhp_pl)
        
    # BHP Calculation Based on HB Model
    BHP_HB = []
    for depth, ecd_hb in zip(TVD, ECD_HB):
        bhp_hb = ((0.052 * depth) * ecd_hb)
        BHP_HB.append(bhp_hb)      
        
    #Mud Window
    TVD_MW = Mud_Window ['TVD(ft)'].tolist()
    PP_ppg = Mud_Window ['PP(ppg)'].tolist()
    FP_ppg = Mud_Window ['FP(ppg)'].tolist()
 
    # Linear interpolation for FP
    FP_interpolated = interp1d(TVD_MW, FP_ppg, kind='linear', fill_value="extrapolate")
    FP_calculated = FP_interpolated(TVD)

    # Linear interpolation for PP
    PP_interpolated = interp1d(TVD_MW, PP_ppg, kind='linear', fill_value="extrapolate")
    PP_calculated = PP_interpolated(TVD) 
    
    FP_ppg =  FP_calculated
    PP_ppg =  PP_calculated       
    # Fracture Pressure Calculation       
    FP=[]
    for depth, fp_ppg in zip(TVD, FP_ppg):
        fp = ((0.052 * depth) * fp_ppg)
        FP.append(fp) 
  
    # Pore Pressure Calculation       
    PP=[]
    for depth, pp_ppg in zip(TVD, PP_ppg):
        pp = ((0.052 * depth) * pp_ppg)
        PP.append(pp) 
        
        
    # Back Pressure Calculation
    Adjusted_BHP_PL = []
    Adjusted_BHP_HB = []  
    BP_PL = []
    BP_HB = []  
    
    for pp, fp, bhp_pl, bhp_hb in zip(PP, FP, BHP_PL, BHP_HB):
        
       # Power Law Model
        if pp < bhp_pl < fp:
           bp_pl = pp + 0.15 * (fp-pp) - bhp_pl
           bhp_pl = pp + 0.15 * (fp-pp)
           Adjusted_BHP_PL.append(bhp_pl)  
           BP_PL.append(bp_pl)       

        elif bhp_pl <= pp:
            bp_pl = pp - bhp_pl + 0.15 * (fp-pp)
            adjusted_bhp_pl = bhp_pl + bp_pl
            Adjusted_BHP_PL.append(adjusted_bhp_pl)
            BP_PL.append(bp_pl)        

        elif bhp_pl >= fp:
            bp_pl = fp - bhp_pl - 0.85 * (fp-pp)
            adjusted_bhp_pl = bhp_pl + bp_pl
            Adjusted_BHP_PL.append(adjusted_bhp_pl)
            BP_PL.append(bp_pl) 

       # HB Model           
        if pp < bhp_hb < fp:
           bp_hb = pp + 0.15 * (fp-pp) - bhp_hb
           bhp_hb = pp + 0.15 * (fp-pp)
           Adjusted_BHP_HB.append(bhp_hb)
           BP_HB.append(bp_hb)  
          
        elif bhp_hb <= pp:
             bp_hb = pp - bhp_hb + 0.15 * (fp-pp)
             adjusted_bhp_hb = bhp_hb + bp_hb
             Adjusted_BHP_HB.append(adjusted_bhp_hb)
             BP_HB.append(bp_hb) 
               
        elif bhp_hb >= fp:
             bp_hb = fp - bhp_hb - 0.85 * (fp-pp)
             adjusted_bhp_hb = bhp_hb + bp_hb
             Adjusted_BHP_HB.append(adjusted_bhp_hb)
             BP_HB.append(bp_hb) 
            

    # Assuming TVD, PP, FP, BHP_PL, BHP_HB, and Adjusted_BHP_HB are defined lists

    plt.figure(figsize=(10, 6),dpi=300)


    # Plot PP
    plt.plot(PP, TVD, label='Pore Pressure', color='blue')

    # Plot FP
    plt.plot(FP, TVD, label='Fracture Pressure', color='red')

    # Plot BHP_PL
    plt.plot(BHP_PL, TVD, label='BHP (Power Law Model)', color='green')

    # Plot BHP_HB
    plt.plot(BHP_HB, TVD, label='BHP (Herschel-Bulkley Model)', color='orange')

    # Plot Adjusted_BHP_HB
    plt.plot(Adjusted_BHP_HB, TVD, label='Adjusted BHP (Herschel-Bulkley Model)', color='purple')
    
    # Plot Adjusted_BHP_PL
    plt.plot(Adjusted_BHP_PL, TVD, label='Adjusted BHP (Power Law Model)', color='yellow')

    # Add labels and title
    plt.xlabel('Pressure (psi)')
    plt.ylabel('TVD (ft)')
    plt.title('Pressure Profiles')
    plt.legend()
    plt.grid(True)

    # Reverse TVD axis
    plt.gca().invert_yaxis()

    # Show plot
    plt.show()
    
    
    Target_BHP_PL = Target_BHP_HB = Adjusted_BHP_HB
    Q_prime_values = range(Q-51, Q+51)

    best_Q_prime_HB = None
    best_Q_prime_PL = None
    closest_adjusted_BHP = float('inf')
    
    pressure_tolerance = 0.001  # Adjust as needed

    closest_pressure_difference_HB = float('inf')
    closest_pressure_difference_PL = float('inf') 
     
    Best_Q_Prime_HB = []
    Best_Q_Prime_PL =[]
    BP_HB_Prime= []
    BP_PL_Prime= []
    
    
    for i in range(len(L1_list)):
        
        P1_prime_list = []
        P2_prime_list = []
        P3_prime_list = []
        P4_prime_list = []
        AF_PL_prime_list = []
        
        P1_HB_prime_list = []
        P2_HB_prime_list = []
        P3_HB_prime_list = []
        P4_HB_prime_list = []
        AF_HB_prime_list = []
        
        
        
        
        for Q_prime in Q_prime_values:
            
            
            # Calculate V_Annulus1
            V_Annulus1_prime = (1.28342246 * Q_prime) / (math.pi * ((Hole_Diameter) ** 2 - (OD_Dc) ** 2))
            
            # Calculate Re_PL1
            Re_PL1_prime = 109000 * ((MW * V_Annulus1_prime ** (2 - n)) / (K)) * ((0.0208 * (Hole_Diameter - OD_Dc)) / (2 + (1 / n))) ** n
            
            # Calculate V_Annulus2
            V_Annulus2_prime = (1.28342246 * Q_prime) / (math.pi * ((Hole_Diameter) ** 2 - (OD_Dp) ** 2))
            
            # Calculate Re_PL2
            Re_PL2_prime = 109000 * ((MW * V_Annulus2_prime ** (2 - n)) / (K)) * ((0.0208 * (Hole_Diameter - OD_Dp)) / (2 + (1 / n))) ** n
        
            # Calculate V_Annulus3
            V_Annulus3_prime = (1.28342246 * Q_prime) / (math.pi * ((ID_Casing) ** 2 - (OD_Dc) ** 2))
            
            # Calculate Re_PL3
            Re_PL3_prime = 109000 * ((MW * V_Annulus3_prime ** (2 - n)) / (K)) * ((0.0208 * (ID_Casing - OD_Dc)) / (2 + (1 / n))) ** n
            
            # Calculate V_Annulus4
            V_Annulus4_prime = (1.28342246 * Q_prime) / (math.pi * ((ID_Casing) ** 2 - (OD_Dp) ** 2))
            
            # Calculate Re_PL4
            Re_PL4_prime = 109000 * ((MW * V_Annulus4_prime ** (2 - n)) / (K)) * ((0.0208 * (ID_Casing - OD_Dp)) / (2 + (1 / n))) ** n



            # Calculate Annular Friction loss (Power Law Model)
            if Re_PL1_prime <= 3470 - 1370 * n:
                Laminar_Delta_P1_prime = (((144 * V_Annulus1_prime * ((2 * n) + 1)) / (((Hole_Diameter - OD_Dc) * 3 * n))) ** n) * (
                        K * L1_list[i] / (300 * (Hole_Diameter - OD_Dc)))
                Turbolant_Delta_P1_prime = 0
              
            elif Re_PL1_prime >= 4270 - 1370 * n:
                f1_prime = (0.0791 / (Re_PL1_prime ** 0.25))
                Turbolant_Delta_P1_prime = (f1_prime * MW * (V_Annulus1_prime ** 2) * L1_list[i]) / (21.1 * (Hole_Diameter - OD_Dc))
                Laminar_Delta_P1_prime = 0
            
           # Append P1 to the list
            P1_prime = Laminar_Delta_P1_prime + Turbolant_Delta_P1_prime
            P1_prime_list.append(P1_prime)

            if Re_PL2_prime <= 3470 - 1370 * n:
               Laminar_Delta_P2_prime = (((144 * V_Annulus2_prime * ((2 * n) + 1)) / (((Hole_Diameter - OD_Dp) * 3 * n))) ** n) * (
                         K * L2_list[i] / (300 * (Hole_Diameter - OD_Dp)))
               Turbolant_Delta_P2_prime = 0
              
            elif Re_PL2_prime >= 4270 - 1370 * n:
                 f2_prime = (0.0791 / (Re_PL2_prime ** 0.25))
                 Turbolant_Delta_P2_prime = (f2_prime * MW * (V_Annulus2_prime ** 2) * L2_list[i]) / (21.1 * (Hole_Diameter - OD_Dp))
                 Laminar_Delta_P2_prime = 0
            
           # Append P2 to the list
            P2_prime = Laminar_Delta_P2_prime + Turbolant_Delta_P2_prime
            P2_prime_list.append(P2_prime)

            if Re_PL3_prime <= 3470 - 1370 * n:
               Laminar_Delta_P3_prime = (((144 * V_Annulus3_prime * ((2 * n) + 1)) / (((ID_Casing - OD_Dc) * 3 * n))) ** n) * (
                         K * L3_list[i] / (300 * (ID_Casing - OD_Dc)))
               Turbolant_Delta_P3_prime = 0
              
            elif Re_PL3_prime >= 4270 - 1370 * n:
                 f3_prime = (0.0791 / (Re_PL3_prime ** 0.25))
                 Turbolant_Delta_P3_prime = (f3_prime * MW * (V_Annulus3_prime ** 2) * L3_list[i]) / (21.1 * (ID_Casing - OD_Dc))
                 Laminar_Delta_P3_prime = 0
            
           # Append P3 to the list
            P3_prime = Laminar_Delta_P3_prime + Turbolant_Delta_P3_prime
            P3_prime_list.append(P3_prime)

            if Re_PL4_prime <= 3470 - 1370 * n:
               Laminar_Delta_P4_prime = (((144 * V_Annulus4_prime * ((2 * n) + 1)) / (((ID_Casing - OD_Dp) * 3 * n))) ** n) * (
                         K * L4_list[i] / (300 * (ID_Casing - OD_Dp)))
               Turbolant_Delta_P4_prime = 0
              
            elif Re_PL4_prime >= 4270 - 1370 * n:
                 f4_prime = (0.0791 / (Re_PL4_prime ** 0.25))
                 Turbolant_Delta_P4_prime = (f4_prime * MW * (V_Annulus4_prime ** 2) * L4_list[i]) / (21.1 * (ID_Casing - OD_Dp))
                 Laminar_Delta_P4_prime = 0
            
           # Append P4 to the list
            P4_prime = Laminar_Delta_P4_prime + Turbolant_Delta_P4_prime
            P4_prime_list.append(P4_prime)
            
            q_prime = 0.0022280093 * Q_prime
            
            Ca_1_prime = 1 - ((1/(1 + n_HB))*(YV_HB/(YV_HB+ K_HB*((2*q_prime*(2*n+1))/((n_HB * math.pi)*((Hole_Diameter-OD_Dc)/24)*(((Hole_Diameter**2)-(OD_Dc**2))/48)))**n_HB)))
            Re_HB1_prime = (4*(2*n_HB+1)/(n_HB))*((MW_HB * (V_Annulus1_prime **(2-n_HB))*((Hole_Diameter-OD_Dc)/24) ** (n_HB))/(YV_HB*((Hole_Diameter-OD_Dc)/(24* V_Annulus1_prime))**n_HB + K_HB *((4*n_HB + 2)/(n_HB * Ca_1_prime))**(n_HB)))

            Ca_2_prime = 1 - ((1/(1 + n_HB))*(YV_HB/(YV_HB+ K_HB*((2*q_prime*(2*n+1))/((n_HB * math.pi)*((Hole_Diameter-OD_Dp)/24)*(((Hole_Diameter**2)-(OD_Dp**2))/48)))**n_HB)))
            Re_HB2_prime = (4*(2*n_HB+1)/(n_HB))*((MW_HB * (V_Annulus2_prime **(2-n_HB))*((Hole_Diameter-OD_Dp)/24) ** (n_HB))/(YV_HB*((Hole_Diameter-OD_Dp)/(24* V_Annulus2_prime))**n_HB + K_HB *((4*n_HB + 2)/(n_HB * Ca_2_prime))**(n_HB)))

            Ca_3_prime = 1 - ((1/(1 + n_HB))*(YV_HB/(YV_HB+ K_HB*((2*q_prime*(2*n+1))/((n_HB * math.pi)*((ID_Casing-OD_Dc)/24)*(((ID_Casing**2)-(OD_Dc**2))/48)))**n_HB)))
            Re_HB3_prime = (4*(2*n_HB+1)/(n_HB))*((MW_HB * (V_Annulus3_prime **(2-n_HB))*((ID_Casing-OD_Dc)/24) ** (n_HB))/(YV_HB*((ID_Casing-OD_Dc)/(24* V_Annulus3_prime))**n_HB + K_HB *((4*n_HB + 2)/(n_HB * Ca_3_prime))**(n_HB)))

            Ca_4_prime = 1 - ((1/(1 + n_HB))*(YV_HB/(YV_HB+ K_HB*((2*q_prime*(2*n+1))/((n_HB * math.pi)*((ID_Casing-OD_Dp)/24)*(((ID_Casing**2)-(OD_Dp**2))/48)))**n_HB)))
            Re_HB4_prime = (4*(2*n_HB+1)/(n_HB))*((MW_HB * (V_Annulus4_prime **(2-n_HB))*((ID_Casing-OD_Dp)/24) ** (n_HB))/(YV_HB*((ID_Casing-OD_Dp)/(24* V_Annulus4_prime))**n_HB + K_HB *((4*n_HB + 2)/(n_HB * Ca_4_prime))**(n_HB)))


            if Re_HB1_prime < Re_HBc:
                Laminar_Delta_P1_HB_prime = ((4*K_HB*L1_list[i])/(14400*(Hole_Diameter-OD_Dc)/12))*((YV_HB/K_HB)+(((32*n_HB+16)/(n_HB*Ca_1_prime*(Hole_Diameter-OD_Dc)/12))*(q_prime/(math.pi*((Hole_Diameter**2)-(OD_Dc**2))/144)))**n_HB)
                Turbolant_Delta_P1_HB_prime = 0
           
            elif Re_HB1_prime > Re_HBc:
                  Fa1_prime = Y * (Ca_1_prime * Re_HB1_prime)**(-Z)
                  Turbolant_Delta_P1_HB_prime = (L1_list[i] * Fa1_prime * MW_HB * q_prime**2)/(1421.22*((Hole_Diameter-OD_Dc)/12)*(((Hole_Diameter**2)-(OD_Dc**2))/144)**2)
                  Laminar_Delta_P1_HB_prime = 0
          
            P1_HB_prime = Laminar_Delta_P1_HB_prime + Turbolant_Delta_P1_HB_prime
            P1_HB_prime_list.append(P1_HB_prime)

            if Re_HB2_prime < Re_HBc:
                Laminar_Delta_P2_HB_prime = ((4*K_HB*L2_list[i])/(14400*(Hole_Diameter-OD_Dp)/12))*((YV_HB/K_HB)+(((32*n_HB+16)/(n_HB*Ca_2_prime*(Hole_Diameter-OD_Dp)/12))*(q_prime/(math.pi*((Hole_Diameter**2)-(OD_Dp**2))/144)))**n_HB)
                Turbolant_Delta_P2_HB_prime = 0
              
            elif Re_HB2_prime > Re_HBc :
                 Fa2_prime = Y * (Ca_2_prime * Re_HB2_prime)**(-Z)
                 Turbolant_Delta_P2_HB_prime = (L2_list[i] * Fa2_prime * MW_HB * q_prime**2)/(1421.22*((Hole_Diameter-OD_Dp)/12)*(((Hole_Diameter**2)-(OD_Dp**2))/144)**2)
                 Laminar_Delta_P2_HB_prime = 0
         
            P2_HB_prime = Laminar_Delta_P2_HB_prime + Turbolant_Delta_P2_HB_prime
            P2_HB_prime_list.append(P2_HB_prime) 
              
            if Re_HB3_prime < Re_HBc:
               Laminar_Delta_P3_HB_prime = ((4*K_HB*L3_list[i])/(14400*(ID_Casing-OD_Dc)/12))*((YV_HB/K_HB)+(((32*n_HB+16)/(n_HB*Ca_3_prime*(ID_Casing-OD_Dc)/12))*(q_prime/(math.pi*((ID_Casing**2)-(OD_Dc**2))/144)))**n_HB)
               Turbolant_Delta_P3_HB_prime = 0
           
            elif Re_HB3_prime > Re_HBc :
                 Fa3_prime = Y * (Ca_3_prime * Re_HB3_prime)**(-Z)
                 Turbolant_Delta_P3_HB_prime = (L3_list[i] * Fa3_prime * MW_HB * q_prime**2)/(1421.22*((ID_Casing-OD_Dc)/12)*(((ID_Casing**2)-(OD_Dc**2))/144)**2)
                 Laminar_Delta_P3_HB_prime = 0
          
            P3_HB_prime = Laminar_Delta_P3_HB_prime + Turbolant_Delta_P3_HB_prime
            P3_HB_prime_list.append(P3_HB_prime)    

            if Re_HB4_prime < Re_HBc:
               Laminar_Delta_P4_HB_prime = ((4*K_HB*L4_list[i])/(14400*(ID_Casing-OD_Dp)/12))*((YV_HB/K_HB)+(((32*n_HB+16)/(n_HB*Ca_4_prime*(ID_Casing-OD_Dp)/12))*(q_prime/(math.pi*((ID_Casing**2)-(OD_Dp**2))/144)))**n_HB)
               Turbolant_Delta_P4_HB_prime = 0
           
            elif Re_HB4_prime > Re_HBc :
                 Fa4_prime = Y * (Ca_4_prime * Re_HB4_prime)**(-Z)
                 Turbolant_Delta_P4_HB_prime = (L4_list[i] * Fa4_prime * MW_HB * q_prime**2)/(1421.22*((ID_Casing-OD_Dp)/12)*(((ID_Casing**2)-(OD_Dp**2))/144)**2)
                 Laminar_Delta_P4_HB_prime = 0
          
            P4_HB_prime = Laminar_Delta_P4_HB_prime + Turbolant_Delta_P4_HB_prime
            P4_HB_prime_list.append(P4_HB_prime)

            # Calculate AF_HB for each iteration
            AF_HB_prime_list = [P1_HB_prime + P2_HB_prime + P3_HB_prime + P4_HB_prime for P1_HB_prime, P2_HB_prime, P3_HB_prime, P4_HB_prime in zip(P1_HB_prime_list, P2_HB_prime_list, P3_HB_prime_list, P4_HB_prime_list)]
            AF_PL_prime_list = [P1_prime + P2_prime + P3_prime + P4_prime for P1_prime, P2_prime, P3_prime, P4_prime in zip(P1_prime_list, P2_prime_list, P3_prime_list, P4_prime_list)]
            
            ECD_PL_prime = []
            for af_pl_prime in AF_PL_prime_list:
                ecd_pl_prime = ((af_pl_prime / (0.052 * TVD[i])) + MW)
                ECD_PL_prime.append(ecd_pl_prime)
                
            ECD_HB_prime = []
            for af_hb_prime in AF_HB_prime_list:
                ecd_hb_prime = ((af_hb_prime / (0.052 * TVD[i])) + MW)
                ECD_HB_prime.append(ecd_hb_prime)    
    
            BHP_PL_prime = []
            for ecd_pl_prime in ECD_PL_prime:
                bhp_pl_prime = ((0.052 * TVD[i]) * ecd_pl_prime)
                BHP_PL_prime.append(bhp_pl_prime)
            
            BHP_HB_prime = []
            for ecd_hb_prime in ECD_HB_prime:
                bhp_hb_prime = ((0.052 * TVD[i]) * ecd_hb_prime)
                BHP_HB_prime.append(bhp_hb_prime)
            
            
            Pressure_Difference_HB=[]
            Pressure_Difference_PL=[]
            
            for bhp_hb_prime, bhp_pl_prime in zip(BHP_HB_prime,BHP_PL_prime):       
       
               pressure_difference_HB = abs(Target_BHP_HB[i] - bhp_hb_prime)
               pressure_difference_PL = abs(Target_BHP_PL[i] - bhp_pl_prime)
           
               Pressure_Difference_HB.append(pressure_difference_HB)
               Pressure_Difference_PL.append(pressure_difference_PL)
               
            if pressure_difference_HB < closest_pressure_difference_HB:
                best_Q_prime_HB = Q_prime
                closest_pressure_difference_HB = pressure_difference_HB
                

            if pressure_difference_PL < closest_pressure_difference_PL:
               best_Q_prime_PL = Q_prime
               closest_pressure_difference_PL = pressure_difference_PL
               
           
 
            if pressure_difference_HB < pressure_tolerance and pressure_difference_PL < pressure_tolerance:
 
               break 
           
        Best_Q_Prime_HB.append(best_Q_prime_HB)
        Best_Q_Prime_PL.append(best_Q_prime_PL)
        BP_PL_Prime.append(pressure_difference_PL)
        BP_HB_Prime.append(pressure_difference_HB)
        

except Exception as e:
    print("An error occurred:", e)

while True:
    try:
        depth_input = int(input("Enter Depth (MD): "))
        if depth_input not in Drilled_depth:
            raise ValueError("Depth not found in drilled depths.")
        
        index = Drilled_depth.index(depth_input)
        data = {
            "MD (ft)": [depth_input],
            "Section": [Section[index]],
            "TVD (ft)": [TVD[index]],
            "PP (psi)": [PP[index]],
            "FP (psi)": [FP[index]],
            "Q (gpm)": [Q],
            "BHP (psi)": [BHP_HB[index]],
            "BP (psi)": [BP_HB[index]],
            "Adjusted BHP (psi)": [Adjusted_BHP_HB[index]],
            "Q' (gpm)": [Best_Q_Prime_HB[index]],
            "BP' (psi)": [BP_HB_Prime[index]]
        }

        df = pd.DataFrame(data)
        print("\nData at Depth", depth_input)
        
        plt.figure(figsize=(12, 8), dpi=300)
        plt.title("Drilling Data at Depth " + str(depth_input), fontsize=20, color='black')
        plt.axis('off')
        
        col_colors = ['lightgrey'] * len(df.columns)
        col_colors[df.columns.get_loc("Q (gpm)")] = 'steelblue'
        col_colors[df.columns.get_loc("BP (psi)")] = 'steelblue'
        col_colors[df.columns.get_loc("BHP (psi)")] = 'steelblue'
        col_colors[df.columns.get_loc("Q' (gpm)")] = 'skyblue'
        col_colors[df.columns.get_loc("BP' (psi)")] = 'skyblue'
        
        table = plt.table(cellText=df.values, colLabels=df.columns, loc='center', cellLoc='center', colColours=col_colors)
        table.auto_set_font_size(False)
        table.set_fontsize(16)
        table.scale(3, 3)
        plt.show()
        
        answer = input("Do you want to enter another depth (Y/N)? ").upper()
        if answer != "Y":
            break
            
    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")