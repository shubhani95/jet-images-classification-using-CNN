import numpy as np
import matplotlib.pyplot as plt
import os

#my_file = 'allcuts.saf'														# path to histogram data file

def plotter(dirname,filename):												# plot histograms from an saf file
	os.system("mkdir "+dirname)
	file = open(filename).read().splitlines()
	for i in range(len(file)):
		line = file[i]
		data = np.array([])
		if "<Description>" in line:
			hist_name = file[i+1].strip()[1:-1]
			[nbins, xmin, xmax] = file[i+3].split()							# extract histogram name, nbins, xmin, xmax
			nbins = int(nbins)
			xmin = float(xmin)
			xmax = float(xmax)
			#xmax = float(xmax) - float(xmax)/float(nbins)					# shift xmax along to align bars properly
			nevents = float(file[i+8][:20].strip())
			sum_weights = float(file[i+9][:20].strip())
			for d in range(nbins):
				data = np.append(data,float(file[i+18+d].strip()[0:12])/100000)	# add data to list
				n_data = [x*nevents/sum_weights for x in data]				# normalise data to number of events
			x = np.linspace(xmin, xmax, nbins)
			fig = plt.figure()
			plt.step(x,n_data,where='mid')
			plt.ylim(bottom=0)
			plt.title(hist_name)
			#plt.show()
			plt.savefig(dirname+'/'+hist_name, dpi=fig.dpi)					# weird formatting issues when saving figure
			plt.close()
	return True

def plotter_compare(dirname,filename1,filename2,filename3,filename4):							# plots histograms from two saf files (comparing the two)
	os.system("mkdir "+dirname)
	file1 = open(filename1).read().splitlines()
	file2 = open(filename2).read().splitlines()
	for i in range(len(file1)):
		line1 = file1[i]
		line2 = file2[i]
		data1 = np.array([])
		data2 = np.array([])
		if "<Description>" in line1:
			hist_name = file1[i+1].strip()[1:-1]
			[nbins, xmin, xmax] = file1[i+3].split()										# extract histogram name, nbins, xmin, xmax
			nbins = int(nbins)
			xmin = float(xmin)
			#xmax = float(xmax) - float(xmax)/float(nbins)
			xmax = float(xmax)																# shift xmax along to align bars properly
			bin_size = float(nbins)/(xmax-xmin)
			nevents1 = float(file1[i+8][:20].strip())
			nevents2 = float(file2[i+8][:20].strip())
			nentries1 = float(file1[i+10][:20].strip())
			nentries2 = float(file2[i+10][:20].strip())
			sum_weights1 = float(file1[i+9][:20].strip())
			sum_weights2 = float(file2[i+9][:20].strip())
			sum_entries1 = float(file1[i+11][:20].strip())
			sum_entries2 = float(file2[i+11][:20].strip())

			for d in range(nbins):
				data1 = np.append(data1,float(file1[i+18+d].strip()[0:12]))					# add data to list
				data2 = np.append(data2,float(file2[i+18+d].strip()[0:12]))
				if "Nbjets" in hist_name:
					n_data1 = [x/(sum_entries1*10000/nevents1) for x in data1]				# normalise data to unity
					n_data2 = [y/(sum_entries2*10000/nevents2) for y in data2]
				else:
					n_data1 = [x/sum_entries1 for x in data1]									# normalise data to unity
					n_data2 = [y/sum_entries2 for y in data2]
			x = np.linspace(xmin, xmax, nbins)
			fig = plt.figure()
			if ("withcuts" in dirname) or ("lowptcut" in  dirname) or ("ptcut30" in dirname):
				if "Nbjets" in hist_name:
					plt.step(x,n_data1,color='r', label='$m_h = 40$GeV\nZero $b$-jets: '+str(1-(nevents1/10000)), where='mid')
					plt.step(x,n_data2,color='b', label='$m_h = 60$GeV\nZero $b$-jets: '+str(1-(nevents2/10000)), where='mid')
				else:
					plt.step(x,n_data1,color='r', label='$m_h = 40$GeV', where='mid')
					plt.step(x,n_data2,color='b', label='$m_h = 60$GeV', where='mid')
			if "variableR_withcuts" in dirname:
				plt.legend(loc='upper right', title="Variable R, Anti-$k_T$\nReduced cuts")
			if "variableR_nocuts" in dirname:
				plt.legend(loc='best', title="Variable R, Anti-$k_T$, no cuts")
			if "lowptcut_40vs60gev_AK4" in dirname:
				plt.legend(loc='best', title="Anti-$k_T$, $R=0.4$\nReduced cuts")
			if "lowptcut_40vs60gev_AK8" in dirname:
				plt.legend(loc='best', title="Anti-$k_T$, $R=0.8$\nReduced cuts")
			if "nocuts_40vs60gev_AK4" in dirname:
				plt.legend(loc='best', title="Anti-$k_T$, $R=0.4$, no cuts")
			if "nocuts_40vs60gev_AK8" in dirname:
				plt.legend(loc='best', title="Anti-$k_T$, $R=0.8$, no cuts")
			if "ptcut30_40vs60gev_AK4" in dirname:
				plt.legend(loc='best', title="Anti-$k_T$, $R=0.4$\nFull cuts ($p_T > 30$GeV)")
			if "ptcut30_40vs60gev_AK8" in dirname:
				plt.legend(loc='best', title="Anti-$k_T$, $R=0.8$\nFull cuts ($p_T > 30$GeV)")
			plt.ylim(bottom=0)
			plt.title(hist_name)
			if "eta" in hist_name:
				plt.xlim(left=xmin)
				plt.xlim(right=xmax)
			else:
				plt.xlim(left=0)
			#plt.legend(loc="upper right")
			plt.ylabel("Events (Normalised to 1)")
			if "bb mass" in hist_name:
				plt.xlabel("$bb$ Mass (GeV)")
				plt.xlim(right=100)
				plt.title("Selected $bb$-jet Pair Mass (GeV)")
			if "b dijet mass" in hist_name:
				plt.xlabel("$m_{bb}$ (GeV)")
				plt.xlim(right=100)
				plt.title("Selected $bb$-jet Pair Mass (GeV)")
			if "bbbb mass" in hist_name:
				plt.xlabel("$m_{bbbb}$ (GeV)")
				plt.xlim(right=250)
				plt.title("$bbbb$ Four-jet Mass (GeV)")
			if "bjet mass1234" in hist_name:
				plt.xlabel("$m_{bbbb}$ (GeV)")
				plt.xlim(right=250)
				plt.title("$bbbb$-jet Mass (GeV)")
			if "jet pt" in hist_name:
				plt.xlabel("$b$-jet $p_T$ (GeV)")
				plt.xlim(right=100)
				plt.title("$b$-jet $p_T$ (GeV)")
			if "Nbjets" in hist_name:
				plt.xlabel("Number of $b$-jets in event")
				plt.title("Event $b$-jet Multiplicity")
			if "Reff" in hist_name:
				plt.title("$b$-jet $R_{eff}$")
				plt.xlabel("$R_{eff}$")
			if "Reff1" in hist_name:
				plt.title("Leading $b$-jet $R_{eff}$")
			if "Reff4" in hist_name:
				plt.title("Lowest $p_T$ $b$-jet $R_{eff}$ (Four $b$-jet Events Only)")
			if "LHiggs pt" in hist_name:
				plt.ylabel("Events (Normalised to 1)")
				plt.title("Light Higgs ($h$) $p_T$")
				plt.xlabel("$h$ $p_T$ (GeV)")
			if "bb DelR" in hist_name:
				plt.ylabel("Events (Normalised to 1)")
				plt.title("$b$-quark Pair $R$ Separation")
				plt.xlabel("$R$")
			if "LHiggs DelR" in hist_name:
				plt.ylabel("Events (Normalised to 1)")
				plt.title("Light Higgs ($h$) $R$ Separation")
				plt.xlabel("$R$")

			plt.savefig(dirname+'/'+hist_name, dpi=fig.dpi)									# weird formatting issues when saving figure
			plt.close()
	return True

def plotter_compare2(dirname,filename1,filename2,filename3,filename4):							# plots histograms from two saf files (comparing the two)
	os.system("mkdir "+dirname)
	file1 = open(filename1).read().splitlines()
	file2 = open(filename2).read().splitlines()
	file3 = open(filename3).read().splitlines()
	file4 = open(filename4).read().splitlines()
	for i in range(len(file1)):
		line1 = file1[i]
		line2 = file2[i]
		line3 = file3[i]
		line4 = file4[i]
		data1 = np.array([])
		data2 = np.array([])
		data3 = np.array([])
		data4 = np.array([])
		if "<Description>" in line1:
			hist_name = file1[i+1].strip()[1:-1]
			[nbins, xmin, xmax] = file1[i+3].split()										# extract histogram name, nbins, xmin, xmax
			nbins = int(nbins)
			xmin = float(xmin)
			#xmax = float(xmax) - float(xmax)/float(nbins)
			xmax = float(xmax)																# shift xmax along to align bars properly
			bin_size = float(nbins)/(xmax-xmin)
			nevents1 = float(file1[i+8][:20].strip())
			nevents2 = float(file2[i+8][:20].strip())
			nevents3 = float(file3[i+8][:20].strip())
			nevents4 = float(file4[i+8][:20].strip())
			nentries1 = float(file1[i+10][:20].strip())
			nentries2 = float(file2[i+10][:20].strip())
			nentries3 = float(file3[i+10][:20].strip())
			nentries4 = float(file4[i+10][:20].strip())
			sum_weights1 = float(file1[i+9][:20].strip())
			sum_weights2 = float(file2[i+9][:20].strip())
			sum_weights3 = float(file3[i+9][:20].strip())
			sum_weights4 = float(file4[i+9][:20].strip())
			sum_entries1 = float(file1[i+11][:20].strip())
			sum_entries2 = float(file2[i+11][:20].strip())
			sum_entries3 = float(file3[i+11][:20].strip())
			sum_entries4 = float(file4[i+11][:20].strip())
			for d in range(nbins):
				data1 = np.append(data1,float(file1[i+18+d].strip()[0:12]))					# add data to list
				data2 = np.append(data2,float(file2[i+18+d].strip()[0:12]))
				data3 = np.append(data3,float(file3[i+18+d].strip()[0:12]))
				data4 = np.append(data4,float(file4[i+18+d].strip()[0:12]))
				if "Nbjets" in hist_name:
					n_data1 = [x/(sum_entries1*100000/nevents1) for x in data1]				# normalise data to unity
					n_data2 = [y/(sum_entries2*100000/nevents2) for y in data2]
					n_data3 = [z/(sum_entries3*100000/nevents3) for z in data3]
					n_data4 = [r/(sum_entries4*100000/nevents4) for r in data4]
				else:
					n_data1 = [x/sum_entries1 for x in data1]									# normalise data to unity
					n_data2 = [y/sum_entries2 for y in data2]
					n_data3 = [z/sum_entries3 for z in data3]
					n_data4 = [r/sum_entries4 for r in data4]
			x = np.linspace(xmin, xmax, nbins)
			fig = plt.figure()
			if "Nbjets" in hist_name:
				plt.step(x,n_data1,color='r', label=r'signal', where='mid')
				plt.step(x,n_data2,color='b', label=r'bbbb', where='mid')
				plt.step(x,n_data3,color='c', label=r'ttbar', where='mid')
				plt.step(x,n_data4,color='g', label=r'zbb', where='mid')
			else:
				plt.step(x,n_data1,color='r', label=r'signal', where='mid')
				plt.step(x,n_data2,color='b', label=r'bbbb', where='mid')
				plt.step(x,n_data3,color='c', label=r'ttbar', where='mid')
				plt.step(x,n_data4,color='g', label=r'zbb', where='mid')

			plt.legend(loc='best', title="Anti-$k_t$, var-$R$")
			plt.ylim(bottom=0)
			plt.title(hist_name)
			if "eta" in hist_name:
				plt.xlim(left=xmin)
				plt.xlim(right=xmax)
			else:
				plt.xlim(left=0)
			#plt.legend(loc="upper right")
			plt.ylabel("Events (Normalised to 1)")
			if "bb mass" in hist_name:
				plt.xlabel("$m_{bb}$ (GeV)")
				plt.xlim(right=1500)
				plt.axvline(x=700, color='g')
				plt.title("double-$b$ tagged fatjet Pair Mass (GeV)")
			if "bjet mass" in hist_name:
				plt.xlabel("$m_{b}$ (GeV)")
				plt.xlim(right=250)
				plt.axvline(x=125, color='g')
				plt.title("double-$b$ tagged fatjet Mass (GeV))")
			if "jet pt" in hist_name:
				plt.xlabel("$b$-jet $p_T$ (GeV)")
				plt.xlim(right=1100)
				plt.title("$b$-jet $p_T$ (GeV)")
			if "Nbjets" in hist_name:
				plt.title("Event Multiplicity")
				plt.xlim(left=1)
				plt.xlabel("Number of double-$b$ tagged fatjets in an event")
			if "bjet pt all" in hist_name:
				plt.xlabel("double-$b$ tagged fatjets $p_T$ (GeV)")
				plt.title("double-$b$ tagged fatjets $p_T$")
			if "LHiggs pt" in hist_name:
				plt.ylabel("Number of Events (Normalised to 1)")
				plt.title("Light Higgs ($h$) $p_T$")
				plt.xlabel("$h$ $p_T$ (GeV)")
			if "bb DelR" in hist_name:
				plt.ylabel("Number of Events (Normalised to 1)")
				plt.title("$b$-quark Pair $R$ Separation")
				plt.xlabel("$R$")
			if "LHiggs DelR" in hist_name:
				plt.ylabel("Number of Events (Normalised to 1)")
				plt.title("Light Higgs ($h$) $R$ Separation")
				plt.xlabel("$R$")
			if "PassSelection" in hist_name:
				plt.xlim(right=1.0)
			if "bjet mass1" in hist_name:
				plt.xlabel("$m_{Leading\ fat-bjet}$ (GeV)")
				plt.title("double-$b$ tagged leading fatjet mass")
				plt.xlim(right=250)
			if "bjet mass2" in hist_name:
				plt.xlabel("$m_{sub-Leading\ fat-bjet}$ (GeV)")
				plt.title("double-$b$ tagged sub-leading fatjet mass")
				plt.xlim(right=250)
			if "bjet DelR" in hist_name:
				plt.xlabel("$R$")
				plt.title("double-$b$ tagged fatjets pair $R$ Separation")
			if "bquark pt" in hist_name:
				plt.ylabel("Number of Events (Normalised to 1)")
				plt.title("$b$-quarks $p_T$")
				plt.xlabel("$b$-quarks $p_T$ (GeV)")
			if "bjet pt1" in hist_name:
				plt.xlabel("double-$b$ tagged leading fatjet $p_T$ (GeV)")
				plt.title("double-$b$ tagged leading fatjet $p_T$")
			if "bjet pt2" in hist_name:
				plt.xlabel("double-$b$ tagged sub-leading fatjet $p_T$ (GeV)")
				plt.title("double-$b$ tagged sub-leading fatjet $p_T$")



			#plt.show()
			plt.savefig(dirname+'/'+hist_name, dpi=fig.dpi)									# weird formatting issues when saving figure
			plt.close()
	return True

def plotter_compare3(dirname,filename1,filename2,filename3,filename4):							# plots histograms from two saf files (comparing the two)
 	os.system("mkdir "+dirname)
 	file1 = open(filename1).read().splitlines()
 	file2 = open(filename2).read().splitlines()
 	file3 = open(filename3).read().splitlines()
 	file4 = open(filename4).read().splitlines()
 	for i in range(len(file1)):
 		line1 = file1[i]
 		line2 = file2[i]
 		line3 = file3[i]
 		line4 = file4[i]
 		data1 = np.array([])
 		data2 = np.array([])
 		data3 = np.array([])
 		data4 = np.array([])
 		if "<Description>" in line1:
 			hist_name = file1[i+1].strip()[1:-1]
 			[nbins, xmin, xmax] = file1[i+3].split()										# extract histogram name, nbins, xmin, xmax
 			nbins = int(nbins)
 			xmin = float(xmin)
 			#xmax = float(xmax) - float(xmax)/float(nbins)
 			xmax = float(xmax)																# shift xmax along to align bars properly
 			bin_size = float(nbins)/(xmax-xmin)
 			nevents1 = float(file1[i+8][:20].strip())
 			nevents2 = float(file2[i+8][:20].strip())
 			nevents3 = float(file3[i+8][:20].strip())
 			nevents4 = float(file4[i+8][:20].strip())
 			nentries1 = float(file1[i+10][:20].strip())
 			nentries2 = float(file2[i+10][:20].strip())
 			nentries3 = float(file3[i+10][:20].strip())
 			nentries4 = float(file4[i+10][:20].strip())
 			sum_weights1 = float(file1[i+9][:20].strip())
 			sum_weights2 = float(file2[i+9][:20].strip())
 			sum_weights3 = float(file3[i+9][:20].strip())
 			sum_weights4 = float(file4[i+9][:20].strip())
 			sum_entries1 = float(file1[i+11][:20].strip())
 			sum_entries2 = float(file2[i+11][:20].strip())
 			sum_entries3 = float(file3[i+11][:20].strip())
 			sum_entries4 = float(file4[i+11][:20].strip())
 			for d in range(nbins):
 				data1 = np.append(data1,float(file1[i+18+d].strip()[0:12]))					# add data to list
 				data2 = np.append(data2,float(file2[i+18+d].strip()[0:12]))
 				data3 = np.append(data3,float(file3[i+18+d].strip()[0:12]))
 				data4 = np.append(data4,float(file4[i+18+d].strip()[0:12]))
 				if "Nbjets" in hist_name:
 					n_data1 = [x/(sum_entries1*100000/nevents1) for x in data1]				# normalise data to unity
 					n_data2 = [y/(sum_entries2*100000/nevents2) for y in data2]
 					n_data3 = [z/(sum_entries3*100000/nevents3) for z in data3]
 					n_data4 = [r/(sum_entries4*100000/nevents4) for r in data4]
 				else:
 					n_data1 = [x/sum_entries1 for x in data1]									# normalise data to unity
 					n_data2 = [y/sum_entries2 for y in data2]
 					n_data3 = [z/sum_entries3 for z in data3]
 					n_data4 = [r/sum_entries4 for r in data4]
 			x = np.linspace(xmin, xmax, nbins)
 			fig = plt.figure()
 			if "Nbjets" in hist_name:
 				plt.step(x,n_data1,color='r', label=r'Signal', where='mid')
 				plt.step(x,n_data2,color='b', label=r'ttbar', where='mid')
 				plt.step(x,n_data3,color='c', label=r'bbbb', where='mid')
 				plt.step(x,n_data4,color='g', label=r'zbb', where='mid')
 			else:
 				plt.step(x,n_data1,color='r', label=r'Signal', where='mid')
 				plt.step(x,n_data2,color='b', label=r'ttbar', where='mid')
 				plt.step(x,n_data3,color='c', label=r'bbbb', where='mid')
 				plt.step(x,n_data4,color='g', label=r'zbb', where='mid')

 			plt.legend(loc='best', title="Anti-$k_t$, $R=0.8$")
 			plt.ylim(bottom=0)
 			plt.title(hist_name)
 			if "eta" in hist_name:
 				plt.xlim(left=xmin)
 				plt.xlim(right=xmax)
 			else:
 				plt.xlim(left=0)
 			#plt.legend(loc="upper right")
 			plt.ylabel("Events (Normalised to 1)")
 			if "bb mass" in hist_name:
 				plt.xlabel("$m_{bb}$ (GeV)")
 				plt.xlim(right=1100)
 				plt.title("Selected $bb$-jet Pair Mass (GeV)")
 			if "b dijet mass" in hist_name:
 				plt.xlabel("$m_{bb}$ (GeV)")
 				plt.xlim(right=1100)
 				plt.title("Selected $bb$-jet Pair Mass (GeV)")
 			if "bbbb mass" in hist_name:
 				plt.xlabel("$m_{bbbb}$ (GeV)")
 				plt.xlim(right=1100)
 				plt.title("$bbbb$ Four-jet Mass (GeV)")
 			if "Two bjet mass" in hist_name:
 				plt.xlabel("Two bjet mass (GeV)")
 				plt.xlim(right=500)
 				plt.title("Two bjet Mass (GeV)")
 			if "jet pt" in hist_name:
 				plt.xlabel("$b$-jet $p_T$ (GeV)")
 				plt.xlim(right=1100)
 				plt.title("$b$-jet $p_T$ (GeV)")
 			if "Nbjets" in hist_name:
 				plt.xlabel("Number of $b$-jets in event")
 				plt.title("Event $b$-jet Multiplicity")
 			if "Reff" in hist_name:
 				plt.title("$b$-jet $R_{eff}$")
 				plt.xlabel("$R_{eff}$")
 			if "Reff1" in hist_name:
 				plt.title("Leading $b$-jet $R_{eff}$")
 			if "Reff4" in hist_name:
 				plt.title("Lowest $p_T$ $b$-jet $R_{eff}$ (Four $b$-jet Events Only)")
 			if "LHiggs pt" in hist_name:
 				plt.ylabel("Number of Events (Normalised to 1)")
 				plt.title("Light Higgs ($h$) $p_T$")
 				plt.xlabel("$h$ $p_T$ (GeV)")
 			if "bb DelR" in hist_name:
 				plt.ylabel("Number of Events (Normalised to 1)")
 				plt.title("$b$-quark Pair $R$ Separation")
 				plt.xlabel("$R$")
 			if "LHiggs DelR" in hist_name:
 				plt.ylabel("Number of Events (Normalised to 1)")
 				plt.title("Light Higgs ($h$) $R$ Separation")
 				plt.xlabel("$R$")
 			if "PassSelection" in hist_name:
 				plt.xlim(right=1.0)
 			if "bjet1 mass" in hist_name:
 				plt.xlabel("$m_{Leading\ Bjet}$ (GeV)")
 				plt.xlim(right=250)
 			if "bjet2 mass" in hist_name:
 				plt.xlabel("$m_{sub-Leading\ Bjet}$ (GeV)")
 				plt.xlim(right=250)
 			if "bjet1 pt" in hist_name:
 				plt.xlabel("$b$-jet $p_T$ (GeV)")
 				plt.xlim(right=500)
 				plt.xlim(left=390)
 				plt.title("$b$-jet $p_T$ (GeV)")

 			#plt.show()
 			plt.savefig(dirname+'/'+hist_name, dpi=fig.dpi)									# weird formatting issues when saving figure
 			plt.close()
 	return True
#plotter('parton_level_40gev','parton_level_40gev.saf')

#plotter_compare('ptcut30_40vs60gev_AK8','ptcut30_40gev_AK8.saf','ptcut30_60gev_AK8.saf')

#plotter_compare2('ttbar_plots_withisr_mpi_pu_AK_jetpileupid_10_const_ptmin50_10k_ma5','check_ttbar_AK_R10_ptmin50_const.saf')
plotter_compare2('check_double_btag_signal_back_varr_withisr_mpi_nopartcut_100k','ml_signal_100k_withmpi_varr.saf','ml_bbbb_100k_withmpi_isr_varr.saf','ml_ttbar_100k_withmpi_isr_varr.saf','ml_zbb_100k_withmpi_isr_varr.saf')
