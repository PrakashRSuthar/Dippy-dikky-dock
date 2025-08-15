
        from pymol import cmd,stored
        
        set depth_cue, 1
        set fog_start, 0.4
        
        set_color b_col, [36,36,85]
        set_color t_col, [10,10,10]
        set bg_rgb_bottom, b_col
        set bg_rgb_top, t_col      
        set bg_gradient
        
        set  spec_power  =  200
        set  spec_refl   =  0
        
        load "data/4lde_20250814_144255.pdb", protein
        create ligands, protein and organic
        select xlig, protein and organic
        delete xlig
        
        hide everything, all
        
        color white, elem c
        color bluewhite, protein
        #show_as cartoon, protein
        show surface, protein
        #set transparency, 0.15
        
        show sticks, ligands
        set stick_color, magenta
        
        
        
        
        # SAS points
 
        load "data/4lde_20250814_144255.pdb_points.pdb.gz", points
        hide nonbonded, points
        show nb_spheres, points
        set sphere_scale, 0.2, points
        cmd.spectrum("b", "green_red", selection="points", minimum=0, maximum=0.7)
        
        
        stored.list=[]
        cmd.iterate("(resn STP)","stored.list.append(resi)")    # read info about residues STP
        lastSTP=stored.list[-1] # get the index of the last residue
        hide lines, resn STP
        
        cmd.select("rest", "resn STP and resi 0")
        
        for my_index in range(1,int(lastSTP)+1): cmd.select("pocket"+str(my_index), "resn STP and resi "+str(my_index))
        for my_index in range(1,int(lastSTP)+1): cmd.show("spheres","pocket"+str(my_index))
        for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_scale","0.4","pocket"+str(my_index))
        for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_transparency","0.1","pocket"+str(my_index))
        
        
        
        set_color pcol1 = [0.361,0.576,0.902]
select surf_pocket1, protein and id [3302,3304,2638,1981,2631,2635,2636,2637,2640,2639,3176,2658,3142,3143,3337,2002,2005,2009,2010,3152,2717,3177,2691,2029,2030,2031,2743,2744,1804,1972,1976,1973,1980,1823,1825,2624,2620,3345,3347,1834,3312,3313,3314,3335,1832,2629,3276,3336,3338,3374,3375,2001,2003] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.341,0.702]
select surf_pocket2, protein and id [4330,4310,4296,4375,4376,4377,3841,4312,4313,3840,3839,1633,1640,1644,1645,1630,1631,1632,1639,4359,2131,4353,2160,2187,3804,2186,2217,2223,2226,2227,2228,2229,2230,2225,4332,2161,2219,3803,3814,2198,3941,3924,3995,3813,3493,4362,1654,1653] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.424,0.361,0.902]
select surf_pocket3, protein and id [3254,3299,3302,3304,3305,2634,2638,2648,2631,2640,2643,192,2540,2541,176,177,2650,2652,194,196,3201,3202,3232,3233,3234,2520,2522,2523,3236,3176,3196,3198,3199,3200,2656,2658,2653,3204,3209,3271,3275,170,180,3274,3276,3277,3278,2627,2630] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.435,0.278,0.702]
select surf_pocket4, protein and id [155,1205,156,308,284,323,326,878,143,151,214,215,1181,1203,239,880,249,883,1173,885,888,240,230,234,233,247,248,283,325] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.698,0.361,0.902]
select surf_pocket5, protein and id [1800,1767,1770,1995,1769,1956,1928,1930,1931,1932,1933,1924,1925,1926,1929,1927,1965,1967,1962,1877,1797,1799,1801,1802,1893] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.651,0.278,0.702]
select surf_pocket6, protein and id [2629,2617,526,527,1853,2616,2630,2606,2607,2612,2615,186,524,525,513,516,518,519,523,2591,2594,2592,2609,185,1860,1825,1826,2538,2536,2559,2567,502,2558,2563,2564,2541] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.361,0.824]
select surf_pocket7, protein and id [3007,2976,2977,2877,3032,2872,2875,2878,3775,2897,3773,2998,2999,4349,3000,4347,2975,3001,3002] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.278,0.533]
select surf_pocket8, protein and id [1432,1439,1463,3343,3344,3346,3349,1412,1413,1441,1833,3315,3316,3345,1834,3313,3321] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.361,0.545]
select surf_pocket9, protein and id [1946,1947,1948,2499,2626,2488,2575,2577,2570,2572,1940,1980,1906,1907,1908,2566] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.278,0.318]
select surf_pocket10, protein and id [1745,1744,1995,1993,1994,2024,2018,2412,2414,2354,2356] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.451,0.361]
select surf_pocket11, protein and id [4368,4385,4382,4367,2994,3023,4333,4334,4335,4336,4337,4341,2970,4320,2969,2995] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.702,0.459,0.278]
select surf_pocket12, protein and id [642,643,77,83,94,95,73,840,808,749,810,819,784] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.902,0.729,0.361]
select surf_pocket13, protein and id [3424,3433,3435,1508,3425,1536,3533,3573,3459,3536,3538,3401,1506,1511] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.702,0.675,0.278]
select surf_pocket14, protein and id [962,856,855,971,997] 
set surface_color,  pcol14, surf_pocket14 
   
        
        deselect
        
        orient
        