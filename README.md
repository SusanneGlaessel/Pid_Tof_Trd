# Particle Identification Framework

## General information
The Pid framework identifies hadrons based on the Tof-*m*<sup>2</sup> and/or the Trd-*dE/dx*. It assigns a probability to be a certain particle species to
every track depending on its *m*<sup>2</sup> and its momentum and/or on its Trd-dEdx and its momentum.

The procedure of particle identification with the  Tof-*m*<sup>2</sup> is described in Sec. 4.4 of this [thesis](https://publikationen.ub.uni-frankfurt.de/opus4/frontdoor/deliver/index/docId/51779/file/main.pdf).
It consists of two steps - determination of particle type hypothesis (fitting) and application of hypothesis to a set of particles (filling).

The procedure of particle identification with the the Trd-*dE/dx* is described in https://cbm-wiki.gsi.de/pub/PWG/Hadron/Pid_with_TRD_dEdx.pdf. It constists of three steps - creating of MC-input, calculating MC-probabilities and assigning probabilities and pid hypothesis to a set of particles (filling). In addition, the separtion of electrons and pions with the RICH detector is applied.

Running of these steps is described below, in section "First Run".

## Installation
ROOT6 is needed for installation.
Follow CERN ROOT [instructions](https://root.cern/install/) to install it.
Version compiled with c++17 flag is preferred, otherwise CMAKE_CXX_STANDARD flag needs to be explicitly specified (see below).

The Pid framework is agnostic to the I/O format.
In the current implementation the [AnalysisTree](https://github.com/HeavyIonAnalysis/AnalysisTree) interface is provided.
However one can use another data format (a proper interface should be developed then).

In order to install the Pid framework together with AnalysisTree interface, one needs either (1) to have pre-installed AnalysisTree or (2) to install the AnalysisTree automatically together with the Pid.
Option (2) is strongly recommended.\
(1) Type before performing the *cmake* command:

    export AnalysisTree_DIR=path-to-analysistree-installation/lib/cmake/AnalysisTree/
OR\
(2) Use following cmake flags when running the *cmake* command:

    -DPID_BUNDLED_AT=ON
    -DPID_BUNDLED_AT_VERSION=v2.3.2

Now perform following commands for the installation:

    git clone https://github.com/HeavyIonAnalysis/Pid.git
    cd Pid
    mkdir build install
    cd build
    source path-to-root-installation/bin/thisroot.sh
    cmake -DCMAKE_INSTALL_PREFIX=../install ../
    make -j install

*path-to-root-installation* and *path-to-analysistree-installation* must be replaced with your actual location of Root and AnalysisTree install directories respectively.

## First run

Before running the Pid framework the environment variables need to be set:

    source path-to-pid-installation/bin/PidConfig.sh

The following steps need to be performed (for more details see below):

- Preparation: creating mc-histograms

   -- for Tof pid: preparation outside of the Pid-framework

   -- for Trd pid run:

		./create_mcinput_trd filelist.sh outfilename (for one analysistree.root)

	and merge output-files for runs into one file:

		hadd mcfilename.root *.mcfilename.root

- Step 1: producing probabilities that are required to identify the tracks.

   -- for Tof pid:

		./run_pid_dcm12

   -- for Trd pid:

		./calculate_mcprob_trd mcfilename
	
- Step 2: assigning probabilities for particle species and a pid-hypothesis to reconstructed tracks from an analysistree and writing them into a new analysistree.

	-- for Tof and Trd pid simoultanously:

		./fill_pid 0 filelist.txt outputfilename pid_file_tof pid_file_trd truncation_mode probability_mode min_hits (purity)

   -- for Tof pid only:

		./fill_pid 1 filelist.txt outputfilename pid_file_tof

   -- for Trd pid only:

		./fill_pid 2 filelist.txt outputfilename pid_file_trd truncation_mode probability_mode min_hits (purity)

The preparation and the first step only need to be performed once.

In the following the above steps are described in more detail.

### Preparation: Creating MC-Input

#### TOF
The MC-input file for a specific system and energy needs to be created by the user before running the Pid framework. The input file needs to contain *m*<sup>2</sup>-*p* histograms in the format like input/m2_pq_vtx.apr20.dcmqgsm.12agev.root. In the folder "reco_info" there is a histogram with all reconstructed tracks matched to TOF hits. In folders "reco_vs_sim_info_from_tof", "reco_vs_sim_info_via_vtx" there are histograms for separate particles species (MC-true) - determined from "TOF hit - simulated particles" and "TOF hit - reconstructed track - simulated particles" matchings respectively.\

Another file which you need contains graphic cuts of certain particles species - it allows to ignore during fitting the particles entries in a non-specific parts of the *m*<sup>2</sup>-*p* histogram, which mainly originate from mismatching between track and TOF hit.
This file is input/cuts.root.\

#### TRD
The MC-input file for a specific system and energy needs to be created by running tasks/create_mcinput_trd.cpp. The MC-input file contains *dE/dx*-*p*  histograms for all reconstructed tracks and hits as well as for reconstructed tracks and hits separated by particle species according to the MC-matching information. Separate histograms are createdf for:

- tracks and hits
- number of Trd-hits (1-4)
- truncation (1-4) (including average)
- all tracks/hits and separated by particles species (MC-matching information)

An example file for 5 million minbias Au+Au events at *p*<sub>beam</sub> = 12 AGeV reconstructed with cbmroot jul25 can be found in input/dEdx_p.jul25.phqmd.auau.12agev.root.

To run the executable type:

	./create_mcinput_trd filelist.sh outfilename

The input for the creation of the MC-histrograms must be an analysistree including mc-information in the format jobId.analysistree.root.

The path and name of the analysistree is read in via filelist.txt. Only one analysistree is read in one run. A separate filelist.sh is required for every run.

The ouput filename will be in the format jobId.mcfilename.root

For reliable MC-histograms a set of run needs to be performed and the resulting files for every run need to get merged into one file with:

    hadd mcfilename.root *.mcfilename.root

### Calculating probabilities

#### TOF: Fitting 
The task which runs the fitting process is tasks/run_pid_dcm12.cpp (look for comments in this file for understanding the configuration example), and it produces a C++ executable install/bin/run_pid_dcm12.
To run this executable just type:

    ./run_pid_dcm12

This run will produce several files, let us digest them.
File pionpos.root contains a set of 1-D histograms named "h_..." where <...> means average value of momentum p, and histograms visualize fitting the pions yield as a function of *m*<sup>2</sup> in a certain momentum range.
Graphs pionpos_0, pionpos_1 and pionpos_2 represent the momentum dependence of fitting parameters of the previous fits: height, mean and sigma of the gaussian.
The height is saved as is, while mean and sigma are fitted.
Finally, the chi2 graph shows the chi2 dependence on momentum when *m*<sup>2</sup> fits were performed.
Similar structre is in files with different particles species (kaonpos, proton, pionneg, kaonneg).\
File allpos.root represents result of fitting all positively charged particles simultaneously.
It contains a set of 1-D histograms named "h_..." where <...> means average value of momentum p, and histograms visualize fitting all particles yields and background as a function of m2 in a certain momentum range.
Graphs pionpos_0, pionpos_1 and pionpos_2, kaonpos_0, ..., bgpos_1, bgpos_2 represent the momentum dependence of fitting parameters of the previous fits for all particles species (height, mean and sigma of the gaussian) and background (three coefficients of 2-nd order polynomial).
Finally, the chi2 graph shows the chi2 dependence on momentum when *m*<sup>2</sup> fits were performed.\
Running the macro/drawFits.C macro on the allpos.root (allneg.root) file allows to build all these distributions on the same canvas.

Another file produced during the fitting procedure is pidtof_getter.root.
There is a Pid::GetterTof which contains results of fitting.
In order to determine bayesian probabilities (purities) of particle species in a certain *m*<sup>2</sup>-*p* point (*p*=3 GeV/*c*, *m*<sup>2</sup>=1 (GeV/*c*<sup>2</sup>)<sup>2</sup>), one runs

    pid_getter_tof->GetBayesianProbability(3, 1)

and obtains result in the following form:

    (std::map<int, double>) { -321 => 0.0000000, -211 => 0.0000000, -1 => 0.0000000, 1 => 0.0011928893, 211 => 8.0296922e-81, 321 => 7.1999527e-46, 2212 => 0.99880711 }
    // Here 1 and -1 stand for background in the right and left side of the histogram respectively.

Running the macro/RunGetterTof.C on the pid_getter_tof.root file produces a set of plots showing purity distribution of different particles species and background.

#### TRD

The task which runs the calculation of mc-probabilities is tasks/calculate_prob_trd.cpp. To run this executable just type:

	./calculate_mcprob_trd mcfilename 

Mc-probabilities for all particle species are calculated for every *dE/dx*-*p*  for all histograms with two methods: Total probability and likelyhood method (see PID_with_Trd_dEdx.pdf for more details). The probabilities are saved in pid_getter_trd.root. The getter is required for the filling procedure. The probabilites for all histograms in mcfilename.root are also saved in *mcfilename*_probabilities.root.

**** Additional information on how to work the pid_getter_trd.root independently of the filling procedure: ****

The Pid::GetterTrd contains the results of the probability calculations. For a reconstructed track with certain *dE/dx*  for each Trd-hit and *p*  value, probabilities can get obtained for the total proability method with:

	pid_getter_trd->GetTrdProbabilities(trdtrack);

and for the likelihood method with:

	pid_getter_trd->GetTrdProbabilitiesMulti(trdtrack);

The trdtrack is an object of the folllowing form:

	TrdContainer trdtrack(mom, pT, charge, nhits_trd, dEdx_hits);

Options for the application of probabilities and pid hypothesis (see below for a more detailed explanation of options) can be set with:

	pid_getter_trd->SetProbabilityMode(mode); // default is total probability method
	pid_getter_trd->SetTruncationMode(mode); // default is average
	pid_getter_trd->SetMinHits(minhits);       // default is 1 hit
	pid_getter_trd->SetPurity(minpurity);       // default is 0.0

For example, for a postive track with an average *dE/dx* = 25 keV/cm and *p* = 5 GeV/*c*, the obtained result with the total probability method is in the form:

	(std::map<int, double>) {(2212, 0.5803), (211, 0.2162), (321, 0.0366), (1000010020, 0.0281), (1000010030, 0.0011), (1000020030, 0.0087), (1000020040 = 0.0009, (-11; 0.1277, (-13, 0.0004), (1, 0.0001) }
	// Here 1 stands for background

Alternatively, probabilites can get applied for the average mode and for the total probability method direcly with:

	getter->GetTrdProbabilitiesTotalAverage(mom, dEdx, charge, nhits_trd);
	// with dEdx = average dE/dx of a track

and for the likelihood method with:
	
	getter->GetTrdProbabilitiesLikelihoodAverage(mom, dEdx_hits, charge);
	// with dEdx_hits = dE/dx vector for hits of a track

The macro/RunGetterTrd.C shows some examples on how to access the pid_getter_trd.root to assign probabilities to a track. I also produces a set of plots showing probability distributions of different particles species and background.

### Filling
Once fitting / probability calculatoin is preformed and the pid_getter_tof.root and/or pid_getter_trd.root files are produced, filling the root-file containing reconstructed tracks can be done - each track will be assigned with probabilities of belonging to different particle species and (optionally) its particle type hypothesis according to the Tof and/or Trd measurement.

This is done in the at_interface/PidFiller.*pp, which are managed by the task tasks/fill_pid.cpp. It produces an executable fill_pid which can be run for either Tof and Trd pid simoultanously or for both separately.

The first argument manages the selection of the detector: =0: Tof and Trd, =1: Tof, =2: Trd.

The following arguments depend on the selection of the detector.

To run Tof and Trd pid simoultanously type:

	./fill_pid 0 filelist.txt outputfilename pid_file_tof pid_file_trd truncation_mode probability_mode min_hits (purity)

To run Tof pid only type:

	./fill_pid 1 filelist.txt outputfilename pid_file_tof

To run Trd pid only type

	./fill_pid 2 filelist.txt outputfilename pid_file_trd truncation_mode probability_mode min_hits (purity)

- filelist.list is a text file containing names of the AnalysisTree root-files to be worked on (an example of AnalysisTree file can be downloaded [here](https://sf.gsi.de/f/3ba5a9e3ff5248edba2c/?dl=1)).

- pid_file_tof / pid_file_tof are the getter-files that were produced in the step before

- for the Trd pid the following options need to get selected:

	-- truncation_mode: modes for calculation of energy loss for up to 4
layers:

	=0: <dEdx> average over all hits (default)

	=1-4: Select hits with lowest dEdx: =1: 1 hit, =2: 2 hits, =3: 3 hits, =4: 4 hits

	-- probability_mode:

	=0: total probability - probability based on particle multiplicites (default)

	=1: likelihood - probability based on dEdx-distribution of particle
species

	-- min_hits: minimum number of required hits per track

	optional:

	-- purity: minimum purity for the pid-hypothesis (most probable particle species) (default purity is 0)

After running this exacutable a pid.analysistree.root file is produced.
It is based on the input file. The branch VtxTracks is replaced with RecParticles branch.
The RecParticles differs from VtxTracks in, firstly, the type of branch (Particles instead of Track), and secondly - in few additional fields:
- for Tof:

	prob_K, prob_d, prob_p, prob_pi, prob_bg - probability of the particle to be a kaon, deuteron, proton, pion or background (undefined type).

	pid - Tof pid hypothesis which is the most probable particle type

- for Trd:

	prob_trd_p, prob_trd_pi, prob_trd_K, prob_trd_d, prob_trd_t, prob_trd_He3, prob_trd_He4, prob_trd_e, prob_trd_mu, prob_trd_bg - probability of the particle to be a proton, pion, kaon, deuteron, triton, He3, He4, electron, myon or background (undefined type).

	pid_trd - Trd pid hypothesis which is the most probable particle type

	In addition, when running with Trd, the RICH electron hyposthesis will be added:

	electron_rich - RICH electron hypothesis: true = electron, false = no electron

### Doxygen documentation
    doxygen docs/Doxyfile
File Doxygen/html/index.html with documentation will be created
