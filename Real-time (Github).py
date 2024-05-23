
import pandas as pd
import math
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import time
import schedule


def run_code():
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
    ROP = pd.read_excel("Data.xlsx", sheet_name='Drilling Parameters', usecols="B", nrows=5, skiprows=4).iloc[0, 0]
    Step = 5
    Time_step = 3600 * (Step/ROP)
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
        Step = 10
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
             Re_HB1 = (4*(2*n_HB+1)/(n_HB))*((MW_HB * (V_Annulus1 **(2-n_HB))*((Hole_Diameter-OD_Dc)/24) ** (n_HB))/(YV_HB*((Hole_Diameter-OD_Dc)/(24* V_Annulus1)) + K_HB *((4*n_HB + 2)/(n_HB * Ca_1))**(n_HB)))

             Ca_2 = 1 - ((1/(1 + n_HB))*(YV_HB/(YV_HB+ K_HB*((2*q*(2*n+1))/((n_HB * math.pi)*((Hole_Diameter-OD_Dp)/24)*(((Hole_Diameter**2)-(OD_Dp**2))/48)))**n_HB)))
             Re_HB2 = (4*(2*n_HB+1)/(n_HB))*((MW_HB * (V_Annulus2 **(2-n_HB))*((Hole_Diameter-OD_Dp)/24) ** (n_HB))/(YV_HB*((Hole_Diameter-OD_Dp)/(24* V_Annulus2)) + K_HB *((4*n_HB + 2)/(n_HB * Ca_2))**(n_HB)))

             Ca_3 = 1 - ((1/(1 + n_HB))*(YV_HB/(YV_HB+ K_HB*((2*q*(2*n+1))/((n_HB * math.pi)*((ID_Casing-OD_Dc)/24)*(((ID_Casing**2)-(OD_Dc**2))/48)))**n_HB)))
             Re_HB3 = (4*(2*n_HB+1)/(n_HB))*((MW_HB * (V_Annulus3 **(2-n_HB))*((ID_Casing-OD_Dc)/24) ** (n_HB))/(YV_HB*((ID_Casing-OD_Dc)/(24* V_Annulus3)) + K_HB *((4*n_HB + 2)/(n_HB * Ca_3))**(n_HB)))

             Ca_4 = 1 - ((1/(1 + n_HB))*(YV_HB/(YV_HB+ K_HB*((2*q*(2*n+1))/((n_HB * math.pi)*((ID_Casing-OD_Dp)/24)*(((ID_Casing**2)-(OD_Dp**2))/48)))**n_HB)))
             Re_HB4 = (4*(2*n_HB+1)/(n_HB))*((MW_HB * (V_Annulus4 **(2-n_HB))*((ID_Casing-OD_Dp)/24) ** (n_HB))/(YV_HB*((ID_Casing-OD_Dp)/(24* V_Annulus4)) + K_HB *((4*n_HB + 2)/(n_HB * Ca_4))**(n_HB)))

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
        
        Depth_MD = pd.read_excel("Data.xlsx", sheet_name='Dynamic Sensor Data', 
                                 usecols="B", nrows=1).iloc[0, 0]
        PP_ppg_sensor = pd.read_excel("Data.xlsx", sheet_name='Dynamic Sensor Data',
                                      usecols="B", nrows=2, skiprows=1).iloc[0, 0] 
        FP_ppg_sensor = pd.read_excel("Data.xlsx", sheet_name='Dynamic Sensor Data',
                                      usecols="B", nrows=3, skiprows=2).iloc[0, 0]
        index = Drilled_depth.index(Depth_MD) 
    
        if PP_ppg[index] != PP_ppg_sensor :
           PP_ppg[index] = PP_ppg_sensor
        if FP_ppg[index] != FP_ppg_sensor :
           FP_ppg[index] = FP_ppg_sensor    
             
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
        
    # Plot Adjusted_BHP_HB
        plt.plot(Adjusted_BHP_PL
                 , TVD, label='Adjusted BHP (Power Law Model)', color='yellow')
    # Add current drilling point
        plt.scatter(Adjusted_BHP_HB[index+1], TVD[index+1], label='Current Drilling Point', color='black', marker='o')
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

    except Exception as e:
        print("An error occurred:", e)

# Run the code every 720 seconds indefinitely
while True:
    print("Running the code...")
    run_code()
    print("Code execution completed.")
    print("Waiting for 720 seconds before running again...")
    time.sleep(720)
    
