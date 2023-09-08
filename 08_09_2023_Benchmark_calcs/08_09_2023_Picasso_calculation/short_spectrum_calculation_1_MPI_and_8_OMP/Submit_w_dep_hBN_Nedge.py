#!/bin/env python3

# import numpy as np
import sys
import os
import subprocess
import shutil


#program = "/Users/apicon/Documents/ProgramasC/Cpp/cbwe/cBWE.out"
program = "sbatch"
path    = os.getcwd()

def print_input(folder,delay,dec, Nk, qTF, dt, eps, cycles, intens, 
	Ncut, t_fin, epsStepAbs, t1, delay_XUV, labelInput,
	start_print_time, end_print_time, step_print):
	# s="%."+str(dec)+"f"
	# delayd = s % (i)
	s="""NEWINPUT

tdse{
tightbinding Haldane_CoreBN_Nedge gap=  7.25 t1=  """ + str(t1)+ """ t2= 0 phi0= 0.0                  
dt  """ + str(dt)+ """ au               
t_fin """ + str(t_fin)+ """ au
dynamical_dt_evolution
epsStepAbs """ + str(epsStepAbs)+ """
Taylor
TaylorOrder 6
}




nkPT{
1 """ + str(int(Nk))+ """ """ + str(int(Nk))+ """
}

Coulomb{
qTF """ + str(qTF)+ """
Ncut """ + str(Ncut)+ """
epsilon_static """ + str(eps)+ """
# G_distance 400
Read_Coulomb_from_files
labelInput """ + labelInput + """

Rytova_Keldysh
}


decoherence{
corehole 0.00396892 au
}



laserpump sin2  {
cycles """ + str((cycles))+ """    
wavelength """ + str((delay))+ """ nm
intensity """ + str(intens)+ """ wcm2
polarization circular 0 1 0 0 0 1
}

laserprobe gaussian {
sigma 0.08 fs
delay  """ + str((delay_XUV))+ """ fs
intensity  10.e+06
frequency 15.06
polarization 1 0 0
}



observables{  
it_resolution 10 
PrintPopulation   """ + str(start_print_time)+ """ """ + str(end_print_time)+ """  """ + str(step_print)+ """                                     
TAbs
Current
}
"""
	outfile = open(str(folder)+"/input.txt",'w')
	outfile.write(s)
	outfile.close()














# main text of program
if (len(sys.argv)!=5):
	print("1st arg: namefolder")
	print("2nd arg: Min")
	print("3rd arg: Max")
	print("4th arg: step")
	quit()

namefolder_init = sys.argv[1]
minT = float(sys.argv[2])
maxT = float(sys.argv[3])
dT   = float(sys.argv[4])

time_au_fs =   0.02418884326505;


Nk = 150
qTF =  0.01
dt = 0.1 # time step


eps = 1.0 # area and epsilon, we divide our results on it. This is in Vacuum Galvani 2019


t_laser = 4; # in fs

t1 = -2.3 # Galvani TB hopping
# t1 = -2.8 # TB hopping crys DFT
# cycles  = 30
intens = 5.e+9 # was 5.e+9 
# intens = 0 #calibrate ATA
epsStepAbs =2.e-2;
labelInput = "G_num502657_Ncut5_"

t_fin = 4 /time_au_fs #1.5 * t_laser /time_au_fs; # in a.u.

start_print_time = 1110.0 /time_au_fs
end_print_time = 2.0/ time_au_fs
step_print = 0.02/ time_au_fs
# t_fin_0 = 0.7*cycles*9.6 / time_au_fs # fs -> au t final - 4 wave pockets, different for different pump pulse
Ncut = 5


delay_XUV = 0#-14.5 #fs
delay_XUV_max = 60.1
delay_XUV_step = 1000.5
while delay_XUV < delay_XUV_max:
	# str_param = str(eps) + "_Nk_" + str(Nk) + "intens"+str(intens/1e+11) + "_Ncut" +str(Ncut) + "_qTF_" + str(qTF) + "_dt_" + str(dt) + "_cycles_" + str(cycles)
	str_param = "UncellArea_" + str(eps) + "_Nk_" + str(Nk) + "epsStepAbs" + str(epsStepAbs) + "intens"+str(intens)+ "_t1_" +str(t1) + "_Ncut" +str(Ncut) + "_laser"+ str(t_laser) + "_qTF_" + str(qTF) + "_dt_" + str(dt)+"_delayXUV_" + str(delay_XUV) + "w_Pump_" 
	# str_param = "NOCOULOMB_Nk_" + str(Nk) + "intens_" +str(intens) + "_t1_" +str(t1) + "_dt_" + str(dt) +"_delayXUV_" + str(delay_XUV) + "w_Pump_" # + "_cycles_" + str(cycles)


	pathCalc    = os.getcwd()
	src_path = pathCalc + "/Coulomb_G_num502657_Ncut5_"
	print(pathCalc)

	namefolder = namefolder_init + str_param
	if (not os.path.isdir(namefolder)):
		os.mkdir(namefolder)
	os.chdir("./" + namefolder)
	pathCalc    = os.getcwd()

	dec  = 2# (int(len(str(dT)))-2) number of symbols after point
	# print("Delays -> [min,max,dT]=["+str(minT)+","+str(maxT)+","+str(dT)+"] fs")

	c_eV = 1240.0



	# we put range in eV 
	i = minT
	while i < maxT: # in range(minT,maxscancel T + 0.1*dT,dT):
		# t_fin = t_fin_0 / i
		W_nm = c_eV/i # from nm to eV 
		s="%s_%."+str(dec)+"f"

		cycles = int(t_laser * i / 4) 

		namefolder_pr = namefolder # + "_cycles_" + str(cycles)
		folder = s % (namefolder_pr,i)
		if (not os.path.isdir(folder)):
			os.mkdir( folder)
			if (not os.path.isdir(folder + "/Output")):
				os.mkdir( folder + "/Output")
				os.mkdir( folder + "/Output/Coulomb")


			path_copy = src_path + "/V_Hartree_" + labelInput +".txt"
			path_paste = pathCalc + "/" + folder + "/Output/Coulomb/V_Hartree_" + labelInput +".txt"

			shutil.copy(path_copy, path_paste)
			path_copy = src_path + "/A_coeff" + labelInput +".txt"
			path_paste = pathCalc + "/" + folder + "/Output/Coulomb/A_coeff" + labelInput +".txt"
			shutil.copy(path_copy, path_paste)

			shutil.copy(src_path + "/B_coeff" + labelInput +".txt", pathCalc + "/" + folder + "/Output/Coulomb/B_coeff" + labelInput +".txt")
			shutil.copy(src_path + "/C_coeff" + labelInput +".txt", pathCalc + "/" + folder + "/Output/Coulomb/C_coeff" + labelInput +".txt")
			shutil.copy(src_path + "/D_coeff" + labelInput +".txt", pathCalc + "/" + folder + "/Output/Coulomb/D_coeff" + labelInput +".txt")
			shutil.copy(src_path + "/Screen_const_" + labelInput +".txt", pathCalc + "/" + folder + "/Output/Coulomb/Screen_const_" + labelInput +".txt")
			print_input(folder,W_nm, dec, Nk, qTF, dt, eps, 
				cycles, intens, Ncut, t_fin, epsStepAbs, t1, delay_XUV, labelInput,
				start_print_time, end_print_time, step_print)
			os.chdir(pathCalc + "/" + folder)
			print(os.getcwd())
			#proc = subprocess.Popen([program,"input.txt"])
			proc = subprocess.Popen([program,"../../Submission_picasso.do"])
			os.chdir(pathCalc)
			proc.wait()

		i+= dT

	os.chdir(path)
	delay_XUV += delay_XUV_step


	#
	#s2=print_E_molcas.do(fich,outfile,"DO Y= " + geom)
	#g2.append(geom)

	# Coulomb_band_reconstruction


	# Coulomb_diag_basis



# Coulomb{
# qTF """ + str(qTF)+ """
# Ncut """ + str(Ncut)+ """
# epsilon_static """ + str(eps)+ """
# G_distance 400
# Read_Coulomb_from_files
# labelInput """ + labelInput + """

# Rytova_Keldysh
# }


# decoherence{
# corehole 0.00396892 au
# }

# observables{   
# PrintPopulation                                            
# TAbs
# Current
# }

# laserpump sin2  {
# cycles 5    
# wavelength 210.16949152542372 nm
# intensity 100000.0 wcm2
# polarization circular 0 1 0 0 0 1                 
# }